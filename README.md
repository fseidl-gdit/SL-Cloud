# Synthetic Lethality Cloud (SL-Cloud)

This project provides a cloud-based data access platform coupled with software and well documented computational notebooks that re-implement published synthetic lethality (SL) inference algorithms to facilitate novel investigation into synthetic lethality. In addition  we provide general purpose functions that support these prediction workflows e.g. saving data in bigquery tables. We anticipate that computationally savvy users can leverage the resources provided in this project to conduct highly customizable analysis based on their cancer type of interest and particular context. 

## Getting Started

### Account Creation
To be able to use our platform, researchers first need to have a Google registered email, a Google Cloud account and have created a Google Project. How to set up the accounts can be found at  [ISB-CGC Quick Start Guide](https://isb-cancer-genomics-cloud.readthedocs.io/en/latest/sections/HowToGetStartedonISB-CGC.html).


### Accessing SL Resource
To  add the syntheticlethality dataset, users need to pin the syntheticlethality project by first clicking "ADD DATA" and after selecting "Pin a project" and "Enter project name", you will see the window as in the Figure below. After writing syntheticlethality into Projectname box, please click on PIN. 

<img src="https://github.com/IlyaLab/SL-Cloud/blob/main/figures/add_sldataset.png" >

### Accessing ISB-CGC Resources
To add ISB-CGC tables, users need to follow the same steps with Accessing SL resource, only difference is, they need to write isb-cgc-bq into Projectname box.
Adding isb-cgc-bq dataset is a requirement to run DAISY pipeline.

### First Notebook

Please run the [first notebook](https://github.com/IlyaLab/SL-Cloud/blob/main/first_notebook.ipynb) to start using our bigquery tables. 

## What is there in the Project
### Scripts
- [SL library](https://github.com/IlyaLab/SL-Cloud/tree/main/scripts/)



### Notebooks
(Guangrong, I guess you updated this part. It looks like copy paste from manuscript, could you please go over it.)

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



