
from google.cloud import bigquery
import pandas as pd

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

