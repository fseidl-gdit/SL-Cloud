# Synthetic Lethality Cloud (SL-Cloud)

This project provides a cloud-based data access platform coupled with software and well documented computational notebooks that re-implement published synthetic lethality (SL) inference algorithms to facilitate novel investigation into synthetic lethality. In addition  we provide general purpose functions that support these prediction workflows e.g. saving data in bigquery tables. We anticipate that computationally savvy users can leverage the resources provided in this project to conduct highly customizable analysis based on their cancer type of interest and particular context. 

## Resource Overview
<img src="https://github.com/IlyaLab/SL-Cloud/blob/main/figures/slhub_overview.png" height="300" width="600">

### Scripts
- [SL library](https://github.com/IlyaLab/SL-Cloud/tree/main/scripts/)


### Notebooks
We provide representative synthetic lethal inference workflows. </br>
Firstly, we reimplemented the published workflow DAISY (Jerby-Arnon et al., 2014) using up-to-date large scale data resources. </br>
Secondly, we implemented a mutation-dependent synthetic lethality prediction (MDSLP) workflow based on the rationale that for tumors with mutations that have an impact on protein expression or protein structure (functional mutation), the knockout effects or inhibition of a partner target gene show conditional dependence for the mutated molecular entities.</br>
Thirdly, we present a workflow that leverages cross-species conservation to infer experimentally-derived synthetic lethal interactions in yeast to predict relevant SL pairs in humans. We implemented the Conserved Genetic Interaction (CGI) workflow based, in part, on methods described in (Srivas et al., 2016) summarized here briefly.

#### Sythetic Lethality Inference Workflows 
Example notebooks can be found in the Example_pipelines directory, which including the following notebooks:
- [DAISY Pipeline](https://github.com/IlyaLab/SL-Cloud/blob/main/Example_pipelines/DAISY_example.ipynb) 
- [Mutation Dependent SL pipeline](https://github.com/IlyaLab/SL-Cloud/blob/main/Example_pipelines/MDSLP_example.ipynb)
- [Conservation-based Inference from Yeast Genetic Interactions](https://github.com/IlyaLab/SL-Cloud/blob/main/Example_pipelines/CGI_example.ipynb)

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
To be able to use our platform, researchers first need to have a Google registered email, a Google Cloud account and have created a Google Project. They can connect the ISB-CGC to their project to be able to use the tables available in ISB-CGC database (Getting isb-cgc-bq dataset is a requirement to run DAISY pipeline). How to use ISB-CGC is explained in detail by videos and notebooks on  [https://isb-cgc.appspot.com]([https://isb-cgc.appspot.com).


### Accessing SL Resource
To  add the syntheticlethality dataset, users need to pin the syntheticlethality project by first clicking "ADD DATA" and after selecting "Pin a project" and "Search for project", you will see the window as in the Figure below. After clicking the project name  syntheticlethality, please click on OPEN. 

<img src="https://github.com/IlyaLab/SL-Cloud/blob/main/figures/add_sldataset.png" height="300" width="600">

## Getting Started
Please run the [first notebook](https://github.com/IlyaLab/SL-Cloud/blob/main/first_notebook.ipynb) to start using our bigquery tables. 

