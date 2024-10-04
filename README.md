# PAC2023Marandi_BSI_prediction
Code repository for the study of machine learning approaches to predict bloodstream infection based on biochemical data by Zargari Marandi et al., 2024.

## Software requirements
data cleaning and preparation in RStudio (R v4.1.2)

Jupyter Notebook for machine leaning in Python 3.10.4

## Description
Firstly, data cleaning and preparation for machine learning are performed using a script called data_cleaning.R. This script handles basic data cleaning tasks to ensure that only the necessary variables and rows are retained for analysis. Following this, another script called data_clean_07022024_pub.ipynb is run, which conducts more detailed data cleaning tasks including merging three datasets. Additionally, this step involves splitting the data by patients.

Secondly, the cleaned datasets (devset and testset) are utilized for machine learning analyses using a script named MAIT_BSI_27092024pub.ipynb. These scripts are extensively commented to ensure readability, and they also provide extra analysis options beyond those reported in the paper.

## Citation
TBA

## Correspondance
Queries shall be forwarded to:
ramtin [DOT] zargari [DOT] marandi [AT] regionh [DOT] dk
