import sys
import numpy as np
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy import stats
from google.cloud import bigquery
import pandas_gbq as gbq
from functools import reduce
import helper_functions
from helper_functions import *

def CoexpressionAnalysis(client, data_resource, input_genes, cor_threshold, p_threshold, adj_method):

    '''
    The gene correlation information is used to detect SL pairs.
    '''

    if data_resource=='PanCancerAtlas':
        table_name='pancancer-atlas.Filtered.EBpp_AdjustPANCAN_IlluminaHiSeq_RNASeqV2_genExp_filtered'
        gene_col_name='Symbol'
        entrez_col_name='Entrez'
        exp_name='normalized_count'
        sample_barcode='SampleBarcode'

    elif data_resource=='CCLE':
        table_name='syntheticlethality.DepMap_public_20Q3.CCLE_gene_expression'
        gene_col_name='Hugo_Symbol'
        exp_name='TPM'
        sample_barcode='DepMap_ID'
        entrez_col_name='Entrez_ID'
       # input_genes = [str(x) for x in input_genes]


    else :
        print("The database name can be either PanCancerAtlas or CCLE")
        return()

    sql_correlation= """ CREATE TEMPORARY FUNCTION tscore_to_p(a FLOAT64, b FLOAT64, c FLOAT64)
     RETURNS FLOAT64
    LANGUAGE js AS
    \"\"\"
    return jStat.ttest(a,b,c); //jStat.ttest( tscore, n, sides)
    \"\"\"
    OPTIONS (
     library="gs://javascript-lib/jstat.min.js"
    );

    WITH
    table1 AS (
    SELECT
    entrez_id,
   (RANK() OVER (PARTITION BY entrez_id ORDER BY data ASC)) + (COUNT(*) OVER ( PARTITION BY entrez_id, CAST(data as STRING)) -  1)/2.0 AS rnkdata,
   ParticipantBarcode
	FROM (
   SELECT
   CAST(__ENTREZ_ID__ AS INT64)  entrez_id,
      AVG( __EXP_NAME__)  AS data,
      __SAMPLE_ID__ AS ParticipantBarcode
   FROM `__TABLE_NAME__`
   WHERE  CAST(__ENTREZ_ID__ AS INT64)  IN (__GENE_LIST__) # labels
         AND __EXP_NAME__ IS NOT NULL
   GROUP BY
      ParticipantBarcode, entrez_id
       )
    )
    ,
    table2 AS (
    SELECT
    entrez_id,
   (RANK() OVER (PARTITION BY entrez_id ORDER BY data ASC)) + (COUNT(*) OVER ( PARTITION BY entrez_id, CAST(data as STRING)) - 1)/2.0 AS rnkdata,
   ParticipantBarcode
    FROM (
   SELECT
      CAST(__ENTREZ_ID__ AS   INT64)  entrez_id,
      AVG(__EXP_NAME__)  AS data,
      __SAMPLE_ID__ AS ParticipantBarcode
   FROM `__TABLE_NAME__`
   WHERE  __ENTREZ_ID__ IS NOT NULL  # labels
         AND __EXP_NAME__ IS NOT NULL
   GROUP BY
      ParticipantBarcode, entrez_id
       )
    )
,
summ_table AS (
SELECT
   n1.entrez_id as entrez_id1,
   n2.entrez_id as entrez_id2,
   COUNT( n1.ParticipantBarcode ) as n,
   CORR(n1.rnkdata , n2.rnkdata) as correlation

FROM
   table1 AS n1
INNER JOIN
   table2 AS n2
ON
   n1.ParticipantBarcode = n2.ParticipantBarcode
   AND n2.entrez_id  NOT IN (__GENE_LIST__)

GROUP BY
   entrez_id1, entrez_id2
UNION ALL
SELECT
   n1.entrez_id as entrez_id1,
   n2.entrez_id as entrez_id2,
   COUNT( n1.ParticipantBarcode ) as n,
   CORR(n1.rnkdata , n2.rnkdata) as correlation

FROM
   table1 AS n1
INNER JOIN
   table1 AS n2
ON
   n1.ParticipantBarcode = n2.ParticipantBarcode
   AND n1.entrez_id <  n2.entrez_id
GROUP BY
   entrez_id1, entrez_id2
)
SELECT *,
   tscore_to_p( ABS(correlation)*SQRT( (n-2)/((1+correlation)*(1-correlation))) ,n-2, 2) as pvalue
   #`cgc-05-0042.Auxiliary.significance_level_ttest2`(n-2, ABS(correlation)*SQRT( (n-2)/((1+correlation)*(1-correlation)))) as alpha
FROM summ_table
WHERE n > 25
#AND correlation > __COR_THRESHOLD__
GROUP BY 1,2,3,4,5
#HAVING pvalue <= __P_THRESHOLD__
ORDER BY entrez_id1 ASC, correlation DESC """

    input_genes = [str(x) for x in input_genes]
    input_genes_for_query= ','.join(input_genes)

    sql_correlation = sql_correlation.replace('__GENE_LIST__', input_genes_for_query)
    sql_correlation = sql_correlation.replace('__P_THRESHOLD__', str(p_threshold))
    sql_correlation = sql_correlation.replace('__COR_THRESHOLD__', str(cor_threshold))
    sql_correlation = sql_correlation.replace('__TABLE_NAME__', table_name)
    sql_correlation = sql_correlation.replace('__GENE_SYMBOL__', gene_col_name)
    sql_correlation = sql_correlation.replace('__ENTREZ_ID__', entrez_col_name)
    sql_correlation = sql_correlation.replace('__EXP_NAME__', exp_name)
    sql_correlation = sql_correlation.replace('__SAMPLE_ID__', sample_barcode)

    results= client.query(sql_correlation).result().to_dataframe()
    report=results[['entrez_id1', 'entrez_id2', 'correlation', 'pvalue']]
    report.columns=['Inactive', 'SL_Candidate', 'Correlation', 'PValue']

    if (adj_method!='none'):
        FDR=multipletests(report['PValue'], alpha=0.05, method= adj_method, is_sorted=False)[1]
        report['PValue']=FDR

    report=report.loc[(report['PValue'] < p_threshold)&(report['Correlation'] > cor_threshold)]
    if report.shape[0]<1:
        return report
    Inactive_Gene_Entrezs=ConvertGene(client,  report.Inactive.unique(), 'EntrezID', ['Gene'])
    SL_Candidate_Entrezs=ConvertGene(client,   report.SL_Candidate.unique(), 'EntrezID', ['Gene'])
    report=pd.merge(Inactive_Gene_Entrezs, report, left_on = 'EntrezID', right_on='Inactive', how = 'inner')
    report=pd.merge(report, SL_Candidate_Entrezs, left_on = 'SL_Candidate', right_on='EntrezID', suffixes=('_Inactive',     '_SL_Candidate'), how = 'inner')
    report=report[['EntrezID_Inactive', 'Gene_Inactive', 'EntrezID_SL_Candidate', 'Gene_SL_Candidate', 'Correlation', 'PValue']]
    #report.drop_duplicates(subset=None, keep='first', inplace=True)
    report=report.groupby('Gene_Inactive').apply(lambda x: x.sort_values('PValue'))
    return report

def SurvivalOfFittest(client, data_source, p_threshold, input_genes, input_mutations, percentile_threshold, cn_threshold,adj_method):
    '''
    Gene expression, Copy Number Alteration, Somatic Mutations are used to decide whether gene is inactive.
    The SL pair detection according to difference in gene effect/dependency score
    given one gene is inactive vs not-inactive
    '''
    if data_source=='PanCancerAtlas':
        gene_exp_table='pancancer-atlas.Filtered.EBpp_AdjustPANCAN_IlluminaHiSeq_RNASeqV2_genExp_filtered'
        mutation_table='pancancer-atlas.Filtered.MC3_MAF_V5_one_per_tumor_sample'
        cn_table='pancancer-atlas.Filtered.all_CNVR_data_by_gene_filtered'
        sample_id='SampleBarcode'
        gene_col_name='Symbol'
        gene_exp='normalized_count'
        cn_gene_name='Gene_Symbol'
        mutation_gene_name='Hugo_Symbol'
        mutation_sample_id='Tumor_SampleBarcode'
        cn_gistic='GISTIC_Calls'
        entrez_id='Entrez'

    elif data_source=='CCLE':
        mutation_table='syntheticlethality.DepMap_public_20Q3.CCLE_mutation'
        gene_exp_table='syntheticlethality.DepMap_public_20Q3.CCLE_gene_expression'
        cn_table='syntheticlethality.DepMap_public_20Q3.CCLE_gene_cn'
        sample_id='DepMap_ID'
        gene_col_name='Hugo_Symbol'
        gene_exp='TPM'
        cn_gene_name='Hugo_Symbol'
        mutation_gene_name='Hugo_Symbol'
        mutation_sample_id='Tumor_Sample_Barcode'
        cn_gistic='CNA'
        cn_threshold=np.log2(2**(cn_threshold)+1)
        entrez_id='Entrez_ID'

    else :
        print("The data source name can be either PanCancerAtlas or CCLE")
        return()


    sql_genes_universe= '''SELECT __SYMBOL_GENE_EXP__ as Hugo_Symbol FROM  __GENE_EXP_TABLE__
    UNION DISTINCT
    SELECT __SYMBOL_CN__ FROM __CN_TABLE__
    UNION DISTINCT
    SELECT __SYMBOL_MUTATION__ FROM  __MUTATION_TABLE__'''

    sql_genes_universe = sql_genes_universe.replace('__SYMBOL_GENE_EXP__', gene_col_name)
    sql_genes_universe = sql_genes_universe.replace('__SYMBOL_CN__', cn_gene_name)
    sql_genes_universe = sql_genes_universe.replace('__SYMBOL_MUTATION__', mutation_gene_name)
    sql_genes_universe = sql_genes_universe.replace('__GENE_EXP_TABLE__', gene_exp_table)
    sql_genes_universe = sql_genes_universe.replace('__CN_TABLE__', cn_table)
    sql_genes_universe = sql_genes_universe.replace('__MUTATION_TABLE__', mutation_table)

    results_all_genes= client.query(sql_genes_universe).result().to_dataframe()

    input_aliases=ConvertGene(client, input_genes, 'EntrezID', ['Gene'])
    input_symbols=pd.merge(results_all_genes, input_aliases, left_on = 'Hugo_Symbol', right_on='Gene', how='inner')
    genes_intermediate_representation = ["'"+str(x)+"'" for x in input_symbols['Hugo_Symbol']]

    input_genes_query= ','.join(genes_intermediate_representation)

    sql_without_mutation= '''
    WITH
    table1 AS (
    (SELECT   entrez_id, Barcode FROM
    (SELECT GE.__EXP_GENE_NAME__ AS entrez_id, GE.__SAMPLE_ID__ AS Barcode ,
    PERCENT_RANK () over (partition by __EXP_GENE_NAME__ order by __GENE_EXPRESSION__ asc) AS Percentile
    FROM  __GENE_EXP_TABLE__ GE
    WHERE GE.__EXP_GENE_NAME__ in (__GENELIST__)) AS NGE
    WHERE NGE.Percentile< __CUTOFFPRC__

    INTERSECT DISTINCT

    SELECT entrez_id ,  Barcode FROM
    (SELECT CN.__CN_GENE_NAME__ AS entrez_id, CN.__SAMPLE_ID__ AS Barcode,
    CN.__CN_GISTIC__ AS NORM_CN
    FROM  __CN_TABLE__ CN
    WHERE CN.__CN_GENE_NAME__ in (__GENELIST__)) AS NC
    WHERE NC.NORM_CN<__CUTOFFSCNA__ )'''


    sql_mutation_part='''

    UNION DISTINCT
    SELECT M.__MUTATION_GENE_NAME__  AS entrez_id , M.__MUTATION_SAMPLE_ID__ AS Barcode
    FROM __MUTATION_TABLE__ M
    WHERE __MUTATION_GENE_NAME__ IN (__GENELIST__) AND
    M.Variant_Classification IN (__MUTATIONLIST__))'''


    rest_of_the_query= '''
     , table2 AS (
    SELECT
        __SAMPLE_ID__ Barcode,  __CN_GENE_NAME__ entrez_id,
        (RANK() OVER (PARTITION BY __CN_GENE_NAME__ ORDER BY __CN_GISTIC__ ASC)) + (COUNT(*) OVER ( PARTITION BY __CN_GENE_NAME__, CAST(__CN_GISTIC__ as STRING)) - 1)/2.0  AS rnkdata
    FROM
       __CN_TABLE__
       where __CN_GENE_NAME__ IS NOT NULL
       ),
summ_table AS (
SELECT
   n1.entrez_id as entrez_id1,
   n2.entrez_id as entrez_id2,
   COUNT( n1.Barcode) as n_1,
   SUM( n2.rnkdata )  as sumx_1,
FROM
   table1 AS n1
INNER JOIN
   table2 AS n2
ON
   n1.Barcode = n2.Barcode
GROUP BY
    entrez_id1, entrez_id2 ),

statistics AS (
SELECT entrez_id1, entrez_id2, n1, n, U1,
      (n1n2/2.0 - U1)/den as zscore

FROM (
   SELECT  entrez_id1, entrez_id2, n_t as n,
       n_1 as n1,
       sumx_1 - n_1 *(n_1 + 1) / 2.0 as U1,
       n_1 * (n_t - n_1 ) as n1n2,
       SQRT( n_1 * (n_t - n_1 )*(n_t + 1) / 12.0 ) as den
   FROM  summ_table as t1
   LEFT JOIN ( SELECT entrez_id, COUNT( Barcode ) as n_t
            FROM table2
            GROUP BY entrez_id)  t2
   ON entrez_id2 = entrez_id
   WHERE n_t > 20 and n_1>3
)
WHERE den > 0
)
SELECT entrez_id1, entrez_id2, n1, n, U1,
    `cgc-05-0042.functions.jstat_normal_cdf`(zscore, 0.0, 1.0 ) as pvalue
FROM statistics
GROUP BY 1,2,3,4,5,6
#HAVING pvalue <= 0.01
ORDER BY pvalue ASC '''


    if input_mutations is None:
        sql_sof=sql_without_mutation +  ')' +' ' +  rest_of_the_query
    else:
        mutations_intermediate_representation = ["'"+x+"'" for x in input_mutations]
        input_mutations_for_query = ','.join(mutations_intermediate_representation)
        sql_sof=sql_without_mutation + ' '+ sql_mutation_part + ' ' +  rest_of_the_query
        sql_sof = sql_sof.replace('__MUTATION_TABLE__', mutation_table)
        sql_sof = sql_sof.replace('__MUTATION_SAMPLE_ID__', mutation_sample_id)
        #sql_sof = sql_sof.replace('__MUTATION_ENTREZ_ID__', mutation_entrez_id)
        sql_sof = sql_sof.replace('__MUTATIONLIST__', input_mutations_for_query)


    sql_sof = sql_sof.replace('__GENELIST__', input_genes_query)
    sql_sof = sql_sof.replace('__CUTOFFPRC__', str(percentile_threshold/100))
    sql_sof = sql_sof.replace('__CUTOFFSCNA__', str(cn_threshold))
    sql_sof = sql_sof.replace('__CN_TABLE__', cn_table)
    sql_sof = sql_sof.replace('__GENE_EXP_TABLE__', gene_exp_table)
    sql_sof = sql_sof.replace('__SAMPLE_ID__', sample_id)
    sql_sof = sql_sof.replace('__ENTREZ_ID__', entrez_id)
    sql_sof = sql_sof.replace('__CN_TABLE__', cn_table)
    sql_sof = sql_sof.replace('__GENE_EXPRESSION__', gene_exp)
    sql_sof = sql_sof.replace('__CN_GISTIC__', cn_gistic)
    sql_sof = sql_sof.replace('__EXP_GENE_NAME__', gene_col_name)
    sql_sof = sql_sof.replace('__CN_GENE_NAME__', cn_gene_name)
    sql_sof = sql_sof.replace('__MUTATION_GENE_NAME__', mutation_gene_name)

    results= client.query(sql_sof).result().to_dataframe()
    report=results[['entrez_id1', 'entrez_id2',  'pvalue']]
    report.columns=['Gene_Inactive', 'Gene_SL_Candidate', 'PValue']

    if (adj_method!='none'):
        FDR=multipletests(report['PValue'], method= adj_method, is_sorted=False)[1]
        report['PValue']=FDR

    report=report.loc[report['PValue'] < p_threshold]
    print(report)
    if report.shape[0]<1:
        return report
    rel_symbols= report.Gene_SL_Candidate.unique()

    candidate_symbols=ConvertGene(client, rel_symbols, 'Gene', ['EntrezID'])
    candidate_symbols.rename(columns={'Gene': 'Gene_SL_Candidate'}, inplace=True)
    inactive_symbols=report.Gene_Inactive.unique()
    inactive_symbols=ConvertGene(client, inactive_symbols, 'Gene', ['EntrezID'])

    report=pd.merge(input_symbols, report, left_on = 'Gene', right_on='Gene_Inactive', how = 'inner')

    report=pd.merge(report[['EntrezID', 'Gene_Inactive', 'Gene_SL_Candidate', 'PValue']], candidate_symbols, left_on = 'Gene_SL_Candidate', right_on='Gene_SL_Candidate',  how = 'inner', suffixes=('_Inactive', '_SL_Candidate'))

    report=report[['EntrezID_Inactive', 'Gene_Inactive', 'EntrezID_SL_Candidate', 'Gene_SL_Candidate', 'PValue']]
   # report.drop_duplicates(subset=None, keep='first', inplace=True)
    report=report.groupby('Gene_Inactive').apply(lambda x: x.sort_values('PValue'))
    return report

def FunctionalExamination(client, database, p_threshold, input_genes, percentile_threshold, cn_threshold, adj_method, input_mutations=None):

    '''
    Gene expression, Copy Number Alteration, Somatic Mutations (optional) are used to decide whether gene is inactive.

    The SL pair detection according to difference in gene effect/dependency score
    given one gene is inactive vs not-inactive
    '''

    if database=='CRISPR':
        mutation_table='syntheticlethality.DepMap_public_20Q3.CCLE_mutation'
        gene_exp_table='syntheticlethality.DepMap_public_20Q3.CCLE_gene_expression'
        cn_table='syntheticlethality.DepMap_public_20Q3.CCLE_gene_cn'
        dep_score_table='syntheticlethality.DepMap_public_20Q3.Achilles_gene_effect'
        sample_info_table='syntheticlethality.DepMap_public_20Q3.sample_info'
        sample_id='DepMap_ID'
        mutation_sample_id='DepMap_ID'
        entrez_id='Entrez_ID'
        mutation_entrez_id='Entrez_Gene_Id'
        cn_threshold=np.log2(2**(cn_threshold)+1)
        gene_exp='TPM'
        effect='Gene_Effect'



    elif database=='siRNA':
        mutation_table='syntheticlethality.DEMETER2_v6.CCLE_mutation'
        gene_exp_table='syntheticlethality.DEMETER2_v6.RNAseq_IRPKM'
        cn_table='syntheticlethality.DEMETER2_v6.WES_snp_cn'
        dep_score_table='syntheticlethality.DEMETER2_v6.D2_combined_gene_dep_score'
        sample_info_table='syntheticlethality.DEMETER2_v6.sample_info'
        sample_id='CCLE_ID'
        mutation_sample_id='Tumor_Sample_Barcode'
        entrez_id='Entrez_ID'
        mutation_entrez_id='Entrez_Gene_Id'
        gene_exp='RPKM'
        effect='Combined_Gene_Dep_Score'


    else :
        print("The database name can be either CRISPR or siRNA")
        return()


    sql_without_mutation= """
    WITH
    table1 AS (
    (SELECT   entrez_id, Barcode FROM
    (SELECT GE.__ENTREZ_ID__ AS entrez_id, GE.__SAMPLE_ID__ AS Barcode ,
    PERCENT_RANK () over (partition by entrez_id order by __GENE_EXPRESSION__ asc) AS Percentile
    FROM  __GENE_EXP_TABLE__ GE
    WHERE GE.__ENTREZ_ID__ in (__GENELIST__)) AS NGE
    WHERE NGE.Percentile< __CUTOFFPRC__

    INTERSECT DISTINCT

    SELECT entrez_id ,  Barcode FROM
    (SELECT CN.__ENTREZ_ID__ AS entrez_id, CN.__SAMPLE_ID__ AS Barcode,
    CN.CNA AS NORM_CN
    FROM  __CN_TABLE__ CN
    WHERE CN.__ENTREZ_ID__ in (__GENELIST__)) AS NC
    WHERE NC.NORM_CN<__CUTOFFSCNA__ )"""


    sql_mutation_part="""

    UNION DISTINCT
    SELECT M.__MUTATION_ENTREZ_ID__  AS entrez_id , M.__MUTATION_SAMPLE_ID__ AS Barcode
    FROM __MUTATION_TABLE__ M
    WHERE __MUTATION_ENTREZ_ID__ IN (__GENELIST__) AND
    M.Variant_Classification IN (__MUTATIONLIST__))"""


    rest_of_the_query= """
     , table2 AS (
    SELECT
        __SAMPLE_ID__ Barcode,  entrez_id,
        (RANK() OVER (PARTITION BY __ENTREZ_ID__ ORDER BY __EFFECT__ ASC)) + (COUNT(*) OVER ( PARTITION BY __ENTREZ_ID__, CAST(__EFFECT__ as STRING)) - 1)/2.0  AS rnkdata
    FROM
       __ACHILLES_TABLE__
       where __ENTREZ_ID__ IS NOT NULL AND __EFFECT__ IS NOT NULL
       ),
summ_table AS (
SELECT
   n1.entrez_id as entrez_id1,
   n2.entrez_id as entrez_id2,
   COUNT( n1.Barcode) as n_1,
   SUM( n2.rnkdata )  as sumx_1,
FROM
   table1 AS n1
INNER JOIN
   table2 AS n2
ON
   n1.Barcode = n2.Barcode
GROUP BY
    entrez_id1, entrez_id2 ),

statistics AS (
SELECT entrez_id1, entrez_id2, n1, n, U1,
       (U1 - n1n2/2.0)/den as zscore
FROM (
   SELECT  entrez_id1, entrez_id2, n_t as n,
       n_1 as n1,
       sumx_1 - n_1 *(n_1 + 1) / 2.0 as U1,
       n_1 * (n_t - n_1 ) as n1n2,
       SQRT( n_1 * (n_t - n_1 )*(n_t + 1) / 12.0 ) as den
   FROM  summ_table as t1
   LEFT JOIN ( SELECT entrez_id, COUNT( Barcode ) as n_t
            FROM table2
            GROUP BY entrez_id)  t2
   ON entrez_id2 = entrez_id
   WHERE n_t > 20 and n_1>3
)
WHERE den > 0
)
SELECT entrez_id1, entrez_id2, n1, n, U1,
    `cgc-05-0042.functions.jstat_normal_cdf`(zscore, 0.0, 1.0 ) as pvalue
FROM statistics
GROUP BY 1,2,3,4,5,6
#HAVING pvalue <= 0.01
ORDER BY pvalue ASC """

    #input_genes=input_genes_for_query.split(',')
    #input_genes=[x.strip('\'') for x in input_genes]

    genes_intermediate_representation = [str(x) for x in input_genes]
    input_genes_query= ','.join(genes_intermediate_representation)

    if input_mutations is None:
        sql_func_ex=sql_without_mutation +  ')' +' ' +  rest_of_the_query
    else:
        mutations_intermediate_representation = ["'"+x+"'" for x in input_mutations]
        input_mutations_for_query = ','.join(mutations_intermediate_representation)
        sql_func_ex=sql_without_mutation + ' '+ sql_mutation_part + ' ' +  rest_of_the_query
        sql_func_ex = sql_func_ex.replace('__MUTATION_TABLE__', mutation_table)
        sql_func_ex = sql_func_ex.replace('__MUTATION_SAMPLE_ID__', mutation_sample_id)
        sql_func_ex = sql_func_ex.replace('__MUTATION_ENTREZ_ID__', mutation_entrez_id)
        sql_func_ex = sql_func_ex.replace('__MUTATIONLIST__', input_mutations_for_query)

    sql_func_ex = sql_func_ex.replace('__GENELIST__', input_genes_query)
    sql_func_ex = sql_func_ex.replace('__CUTOFFPRC__', str(percentile_threshold/100))
    sql_func_ex = sql_func_ex.replace('__CUTOFFSCNA__', str(cn_threshold))
    sql_func_ex = sql_func_ex.replace('__CN_TABLE__', cn_table)
    sql_func_ex = sql_func_ex.replace('__SAMPLE_INFO_TABLE__', sample_info_table)
    sql_func_ex = sql_func_ex.replace('__GENE_EXP_TABLE__', gene_exp_table)
    sql_func_ex = sql_func_ex.replace('__SAMPLE_ID__', sample_id)
    sql_func_ex = sql_func_ex.replace('__ENTREZ_ID__', entrez_id)
    sql_func_ex = sql_func_ex.replace('__ACHILLES_TABLE__', dep_score_table)
    sql_func_ex = sql_func_ex.replace('__GENE_EXPRESSION__', gene_exp)
    sql_func_ex = sql_func_ex.replace('__EFFECT__', effect)

    results= client.query(sql_func_ex).result().to_dataframe()
    report=results[['entrez_id1', 'entrez_id2',  'pvalue']]
    report.columns=['Inactive', 'SL_Candidate', 'PValue']

    if (adj_method!='none'):
        FDR=multipletests(report['PValue'],  method= adj_method, is_sorted=False)[1]
        report['PValue']=FDR

    report=report.loc[(report['PValue'] < p_threshold)]#&
    if report.shape[0]<1:
        return report

    rel_symbols= report.SL_Candidate.unique()

    candidate_symbols=ConvertGene(client, rel_symbols, 'Gene', ['EntrezID'])
    #candidate_symbols.rename(columns={'Gene': 'Gene_SL_Candidate'}, inplace=True)

    if database=='siRNA':
        int_vals=np.array([x.isdecimal() for x in  report['SL_Candidate']])
        report.drop(report.index[np.where(int_vals==False)], inplace=True)
        report=report.astype({'SL_Candidate': 'int64'})

   #report['#Total Samples']=np.repeat(n,report.shape[0])
    Inactive_Gene_Entrezs=ConvertGene(client,  report.Inactive.unique(), 'EntrezID', ['Gene'])
    SL_Candidate_Entrezs=ConvertGene(client,  report.SL_Candidate.unique(), 'EntrezID', ['Gene'])
    report=pd.merge(Inactive_Gene_Entrezs, report, left_on = 'EntrezID', right_on='Inactive', how = 'inner')
    report=pd.merge(report, SL_Candidate_Entrezs, left_on = 'SL_Candidate', right_on='EntrezID', suffixes=('_Inactive',  '_SL_Candidate'),  how = 'inner')
    report=report[['EntrezID_Inactive', 'Gene_Inactive', 'EntrezID_SL_Candidate', 'Gene_SL_Candidate', 'PValue']]
    #report.drop_duplicates(subset=None, keep='first', inplace=True)
    report=report.groupby('Gene_Inactive').apply(lambda x: x.sort_values('PValue'))
    return report

def UnionResults(results):
    '''
    This functions merges results from the same inference procedure applied on different datasets. A combined p-value is returned
    '''

    for i in range(len(results)):
        if results[i].shape[0]<1:
            print("At least one of the dataframes is empty, please run only with nonempty dataframes");
            return()
        results[i].reset_index(inplace=True, drop=True)
        if  'Correlation' in results[i].columns:
                results[i].rename(columns = {'Correlation':'Correlation_'+ str(i)}, inplace = True)
        results[i].rename(columns = {'PValue':'PValue_'+ str(i)}, inplace = True)

    combined_results=results[0]
    for i in range(1,len(results)):
        combined_results =pd.merge(combined_results, results[i], on = ['EntrezID_Inactive','Gene_Inactive', 'EntrezID_SL_Candidate', 'Gene_SL_Candidate'], how = 'outer')

    rel_cols=combined_results.columns[np.array([x.startswith('PValue') for x in combined_results.columns])]
    p_matrix=combined_results[rel_cols]

    agg_p=p_matrix.apply(lambda x:stats.combine_pvalues(x.dropna().tolist(), method='fisher', weights=None)[1], axis=1)
    combined_results['PValue']=agg_p
    return(combined_results)

def MergeResults(results):
    '''
    The results from SoF, Coexpression and Fuctional Screening analysis (each of them is optional) are merged,
    a combined p-value is returned for each Synthetic Lethal pair.
    '''
    for i in range(len(results)):
        if results[i].shape[0]<1:
            print("At least one of the dataframes is empty, please run only with nonempty dataframes");
            return();
        results[i].reset_index(inplace=True, drop=True)
        results[i].rename(columns = {'PValue':'PValue'+ str(i)}, inplace = True)

    combined_results=results[0]
    for i in range(1,len(results)):
        combined_results =pd.merge(combined_results, results[i], on = ['EntrezID_Inactive','Gene_Inactive', 'EntrezID_SL_Candidate', 'Gene_SL_Candidate'], how = 'inner')

    rel_cols=combined_results.columns[np.array([x.startswith('PValue') for x in combined_results.columns])]
    p_matrix=combined_results[rel_cols]

    agg_p=p_matrix.apply(lambda x:stats.combine_pvalues(x.dropna().tolist(), method='fisher', weights=None)[1], axis=1)
    combined_results['PValue']=agg_p
    return(combined_results[['EntrezID_Inactive','Gene_Inactive', 'EntrezID_SL_Candidate', 'Gene_SL_Candidate', 'PValue']])
