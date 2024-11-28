# Long-read Pan-Cancer Transcriptomics Unveil Distinct Trends in Gene and Isoform Expression Alterations in Tumorigenesis
This repository contains the Jupyter notebooks and associated R scripts used to generate figures for this paper. Below is a description of each notebook and the figures they correspond to.

## Repository Structure

### Data
#### bambu_result
Outputs of bambu, including quantification of genes and transcripts, and the gtf file.

#### jaffal_result
Summrised result of JAFFAL.
### Code
#### fig2d.ipynb
This notebook generates Figure 2d.

#### fig3&5.ipynb
This notebook is used to create the most subfigures of Figures 3 and 5, and the data processing of these figures.

#### fig3d.ipynb
This notebook specifically focuses on generating Figure 3d.

#### fig4.ipynb
This notebook is responsible for Figure 4, containing filterings of fusion transcripts and code of figure4b, 4c, 4e&f.

#### fig6&hierarchy_system.ipynb
This notebook generates Figure 6 and includes the hierarchical system analysis.

#### new_group_gene_deseq_fdr.R
Script for differential gene expression analysis using DESeq with FDR adjustment of ten tissues.

#### new_all.R
Performs global analysis for the dataset.

#### transplotr.R
Visualizes structures of isoforms.

#### fig_5f.R
Generates Figure 5f.

#### bar_kegg.R
Produces bar plots for KEGG pathway analysis.

## Prerequisites
To run these notebooks, you need the following:

Python version: >= 3.11.0

Dependencies: Install the required Python packages by running:

pip install -r requirements.txt

Jupyter Notebook: Ensure you have Jupyter Notebook or JupyterLab installed to run .ipynb files.



## Notes
Ensure that the datasets required for the analyses are available and accessible. Update file paths in the notebooks as necessary.
For any questions or clarifications regarding the code or figures, please contact dingyy7@mail2.sysu.edu.cn
