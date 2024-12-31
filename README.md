# Long-read Pan-Cancer Transcriptomics Unveil Distinct Trends in Gene and Isoform Expression Alterations in Tumorigenesis
This repository contains the Jupyter notebooks and associated R scripts used to generate figures for this paper. Below is a description of each notebook and the figures they correspond to.

## Repository Structure

### Data
#### bambu_result
Outputs of bambu, including quantification of genes and transcripts, and the gtf file.

#### sange ab1 files
AB1 files from Sanger sequencing used to validate fusion transcripts.

### Code
#### fig2d.ipynb
This notebook generates Figure 2d.

#### fig3&5.ipynb
This notebook is used to filter out outliers and create the most subfigures of Figures 3 and 5.

#### fig3d.ipynb
This notebook specifically focuses on generating Figure 3d.

#### fig4.ipynb
This notebook is responsible for Figure 4, containing filterings of fusion transcripts and code of figure4b, 4c, 4e&f. Besides, this notebook generates supplementary tables 12-14.

#### fig6&hierarchy_system.ipynb
This notebook generates Figure 6 and includes the hierarchical system analysis.

#### fig5d_5e&5f.ipynb
This notebook analyzes significantly different isoform ratios and generates related figures, including panels D, E, and F of Figure 5.

#### supplementart_tables.ipynb
This notebook generates most supplementary tables, including supplementary tables 4, 6, 8, 9, 10 and the data for ploting figure 3f.

#### figs6.ipynb
This notebook is used for ploting supplementary figure 6.

#### new_group_gene_deseq_fdr.R
Script for differential gene expression analysis using DESeq with FDR adjustment of ten tissues.

#### new_all.R
Performs global analysis for the whole dataset.

#### transplotr.R
Visualizes structures of isoforms.

#### fig_5f.R
Generates Figure 5f.

#### KEGG&GOBP_plot.R
Produces dot plots and bar plots for KEGG and GO_bp pathway analysis.

## Prerequisites
To run these notebooks, you need the following:

Python version: >= 3.11.0

Dependencies: Install the required Python packages by running:

pip install -r requirements.txt

Jupyter Notebook: Ensure you have Jupyter Notebook or JupyterLab installed to run .ipynb files.



## Notes
For any questions or clarifications regarding the code or figures, please contact dingyy7@mail2.sysu.edu.cn
