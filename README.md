# Long-read Pan-Cancer Transcriptomics Unveil Distinct Trends in Gene and Isoform Expression Alterations in Tumorigenesis

This repository contains the Jupyter notebooks and associated R scripts used to generate figures for this paper. Below is a description of each notebook and the figures they correspond to.

## Repository Structure

### Data

#### bambu\_result

Outputs of bambu, including quantification of genes and transcripts, and the gtf file.

#### jaffal\_result

Summrised result of JAFFAL.

### Code

#### fig2d\_revised.ipynb

This notebook generates Figure 2d.

#### fig3\&5\_revised.ipynb

This notebook is used to filter out outliers and create the most subfigures of Figures 3 and 5.

#### fig3d\_revised.ipynb

This notebook specifically focuses on generating Figure 3d.

#### fig4\_revised.ipynb

This notebook is responsible for Figure 4, containing filterings of fusion transcripts and code of figure4b, 4c, 4e\&f.

#### fig6\&hierarchy\_system\_revised.ipynb

This notebook generates Figure 6 and includes the hierarchical system analysis.

#### new\_group\_gene\_deseq\_fdr\_revised.R

Script for differential gene expression analysis using DESeq with FDR adjustment of ten tissues.

#### new\_all\_revised.R

Performs global analysis for the whole dataset.

#### transplotr\_revised.R

Visualizes structures of isoforms.

#### fig\_5f\_revised.R

Generates Figure 5f.

#### KEGG\&GOBP\_plot\_revised.R

Produces dot plots and bar plots for KEGG and GO\_bp pathway analysis.

## Prerequisites

To run these notebooks, you need the following:

Python version: >= 3.11.0

Dependencies: Install the required Python packages by running:

pip install -r requirements.txt

Jupyter Notebook: Ensure you have Jupyter Notebook or JupyterLab installed to run .ipynb files.



## Notes

Ensure that the datasets required for the analyses are available and accessible. Update file paths in the notebooks as necessary.
For any questions or clarifications regarding the code or figures, please contact dingyy7@ihcams.ac.cn

