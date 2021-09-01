# Synthetic Lethality Cloud (SL-Cloud)

This project provides a cloud-based data access platform coupled with software and well documented computational notebooks that re-implement published synthetic lethality (SL) inference algorithms to facilitate novel investigation into synthetic lethality. In addition  we provide general purpose functions that support these prediction workflows e.g. saving data in bigquery tables. We anticipate that computationally savvy users can leverage the resources provided in this project to conduct highly customizable analysis based on their cancer type of interest and particular context. 

## Resource Overview
<img src="https://github.com/IlyaLab/SL-Cloud/blob/main/figures/slhub_overview.png" height="300" width="600">

### Scripts
- [SL library](https://github.com/IlyaLab/SL-Cloud/tree/main/scripts/)


### Notebooks
We provide representative synthetic lethal inference workflow based on highly cited published workflows

#### Sythetic Lethality Inference Workflows 

- [DAISY Pipeline](https://github.com/IlyaLab/SL-Cloud/blob/main/DAISY_pipeline/DAISY_from_library.ipynb) 
- [Mutation-based Conditional SL pipeline](https://github.com/IlyaLab/SL-Cloud/blob/main/mutation_dependent_SL_pipeline/mutation_dependent_SL_pipeline.ipynb)
- [Conservation-based Inference from Yeast Genetic Interactions](https://github.com/bhrtrcn/SyntheticLethality/blob/c7bf444b2eece46777dd545b52f18cd4150d0153/Notebooks/leveraging_conservation_pipeline/YeastOrtholog_SL_pairs.ipynb) (need to transfer here, links to bhrtrcn github)

#### Data Wrangling and Cleaning Associated Procedures 
- Dataset Creation
- Table Creation
- Writing into Excel File
- Gene conversion among gene symbol, EntrezID and alias 

### Synthetic-Lethality Inference Data Resources
This resource provides access to publicly available cancer genomics datasets relevant for SL inference. These data have been pre-processed, cleaned and stored in cloud-based query-able tables leveraging [Google BigQuery](https://cloud.google.com/bigquery)  technology. In addition we leverage relevant datasets available through the Institute for Systems Biology Cancer Genomics Cloud ([ISB-CGC](https://isb-cgc.appspot.com/)) to make inferences of potential synthetic lethal interactions. 
The following represent project-specific datasets with relevance for SL inference:

- **DEPMAP**: DEPMAP shRNA (DEMETER2 V6) and CRISPR (DepMap Public 20Q3) gene expression, sample information, mutation and copy number alterations  for CRISPR experiments and and gene dependency scores for shRNA and gene effect scores.

- **CellMap**: Yeast interaction dataset based on fitness scores after single and double knockouts from SGA experiements.

- **Gene Information**: Tables with relevant gene annotation information such as yeast and human ortholog information, gene-alias-Entrez ID mapping, gene Ensembl-id mapping, gene-Refseq mapping.


### Account Creation
To be able to use our platform, researchers first need to have a Google registered email, a Google Cloud account and have created a Google Project. They can connect the ISB-CGC to their project to be able to use the tables available in ISB-CGC database (Getting isb-cgc-gq dataset is a requirement to run DAISY pipeline). How to use ISB-CGC is explained in detail by videos and notebooks on  [https://isb-cgc.appspot.com]([https://isb-cgc.appspot.com).


### Accessing SL Resource
To  add the syntheticlethality dataset, users need to pin the syntheticlethality project by first clicking "ADD DATA" and after selecting "Pin a project" and "Search for project", you will see the window as in the Figure below. After clicking the project name  syntheticlethality, please click on OPEN. 

<img src="https://github.com/IlyaLab/SL-Cloud/blob/main/figures/add_sl_dataset.png" height="300" width="600">

## Getting Started
Please run this notebook [First notebook] (https://github.com/IlyaLab/SL-Cloud/blob/first_notebook.ipynb) to start using our bigquery tables. 

