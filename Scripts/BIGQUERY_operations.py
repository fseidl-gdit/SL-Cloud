
import sys
import numpy as np
import pandas as pd
from google.cloud import bigquery
import pandas_gbq as gbq

def CreateDataSet(client, dataset_name, project_id, dataset_description):

    '''
    This function creates a dataset named dataset_name into the project given
    project_id, with the data_description provided
    '''

    dataset_id = client.dataset(dataset_name, project=project_id)
    try:
        dataset=client.get_dataset(dataset_id)
        print('Dataset {} already exists.'.format(dataset.dataset_id))
    except:
        dataset = bigquery.Dataset(dataset_id)
        dataset = client.create_dataset(dataset)
        dataset.description =dataset_description
        dataset = client.update_dataset(dataset, ["description"])
        print('Dataset {} created.'.format(dataset.dataset_id))



def CreateTable(client, data, dataset_name, table_name, project_id, table_desc, table_annotation=None):
    '''
     This function creates a dataset named dataset_name into the project given
     project_id, with the data_description provided
    '''

    dataset_id = client.dataset(dataset_name, project=project_id)
    try:
        dataset=client.get_dataset(dataset_id)
        if table_annotation is None:
            gbq.to_gbq(data, dataset.dataset_id +'.'+ table_name, project_id=project_id, if_exists='replace')

        else:
            gbq.to_gbq(data, dataset.dataset_id +'.'+ table_name, project_id=project_id, table_schema = table_annotation , if_exists='replace')
        print("Table created successfully")

    except:
        print('Table could not be created')
    try:
        table=client.get_table(dataset.dataset_id +'.'+ table_name)
        table.description =table_desc
        table = client.update_table(table, ["description"])
       # print("Table description added successfully")
    except:
        print('Table description could not be updated')
