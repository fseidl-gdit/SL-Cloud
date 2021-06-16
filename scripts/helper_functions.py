import sys
import numpy as np
import pandas as pd
from google.cloud import bigquery
import pandas_gbq as gbq
from functools import reduce


def GetMutation_Frequency(client, tissue_types_query, selected_variants_query):
    '''
    This function returns the number of mutant and all samples given the tissue and variant types, tissue_types_query,   selected_variants_query must be vectors
    '''

    intermediate_tissue_types = ["'"+str(x)+"'" for x in tissue_types_query]
    tissue_types_query= ','.join(intermediate_tissue_types)
    intermediate_selected_variants = ["'"+str(x)+"'" for x in selected_variants_query]
    selected_variants_query= ','.join(intermediate_selected_variants)

    query='''WITH
TABLE1 AS (SELECT Study, Hugo_Symbol, count(ParticipantBarcode) AS Mutated_Count
FROM `pancancer-atlas.Filtered.MC3_MAF_V5_one_per_tumor_sample` WHERE Study in (__TISSUE_TYPE__) and FILTER='PASS'
AND Variant_Classification in (__SELECTED_VARIANTS__)
GROUP BY Study, Hugo_Symbol
ORDER BY Study, count(Hugo_Symbol) DESC),

TABLE2 AS (SELECT STUDY, COUNT(DISTINCT ParticipantBarcode)  AS ALL_SAMPLES_COUNT FROM
  `pancancer-atlas.Filtered.MC3_MAF_V5_one_per_tumor_sample` WHERE Study in (__TISSUE_TYPE__) AND FILTER='PASS'
 AND Variant_Classification in (__SELECTED_VARIANTS__)
 GROUP BY Study
 ORDER BY Study
)

SELECT TABLE1.STUDY, TABLE1.HUGO_SYMBOL, TABLE1.MUTATED_COUNT, TABLE2.ALL_SAMPLES_COUNT, (TABLE1.MUTATED_COUNT/TABLE2.ALL_SAMPLES_COUNT)*100 PERCENTAGE
FROM TABLE1, TABLE2 WHERE TABLE1.STUDY=TABLE2.STUDY
ORDER BY TABLE1.STUDY, PERCENTAGE DESC'''

    query=query.replace('__TISSUE_TYPE__', tissue_types_query)
    query=query.replace('__SELECTED_VARIANTS__', selected_variants_query)

    results= client.query(query).result().to_dataframe()
    return(results)


def ConvertGene(client, input_vector, input_type, output_type):
    '''
    This function provides conversion between EntrezID, Gene and Alias
    Input type can be one of 'Alias', 'Gene', 'EntrezID'
    output type must a vector like ['Gene', 'EntrezID']
    '''


    sql='''
    SELECT DISTINCT __IN_TYPE__,  __OUT_TYPE__
    FROM  `syntheticlethality.gene_information.gene_info_human`
    where  __IN_TYPE__  in (__IN_VECTOR__)
    '''

    if input_type=='EntrezID':
        intermediate_representation = [str(x) for x in input_vector]
    else:
        intermediate_representation = ["'"+str(x)+"'" for x in input_vector]

    input_vector_query= ','.join(intermediate_representation)

    out_type_intermediate_representation = [str(x) for x in output_type]
    output_type_for_query= ','.join(out_type_intermediate_representation)

    sql=sql.replace('__OUT_TYPE__', output_type_for_query)
    sql=sql.replace('__IN_TYPE__', input_type)
    sql=sql.replace('__IN_VECTOR__', input_vector_query)

    result= client.query(sql).result().to_dataframe()
    return(result)


def WriteToExcel(excel_file, data_to_write, excel_tab_names):
    '''
    This function writes the dataframes whose names are given
    in data_to-write parameter to the excel files whose names
    given in excel_file_names parameter
    '''
    with pd.ExcelWriter(excel_file) as writer:
        for i in range(len(excel_tab_names)):
            data_to_write[i].to_excel(writer, sheet_name=excel_tab_names[i], index=False)
