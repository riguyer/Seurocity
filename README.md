# Seurocity v1.0.0
*Easy velocity analysis using scVelo with data preprocessed in Seurat*

### Background
Velocyto (https://velocyto.org/) is a powerful tool for studying dynamic gene expression changes within single cell datasets. The Velocyto program is a command line tool that processes raw single cell sequencing data. Velocyto can accept data output by 10X Genetics Cell Ranger pipeline, as well as other common single cell analysis tools, such as SmartSeq2. Advanced users can also provide the required inputs files manually. The program outputs a Loom file (https://linnarssonlab.org/loompy/) containing spliced and unspliced gene counts for every cell. Multiple tools exist for analyzing and processing the data output by Velocyto, including velocyto.R and scVelo (https://scvelo.readthedocs.io/index.html). 

Seurat (https://satijalab.org/seurat/) is a powerful and widely-used R-based backaged for processing, analyzing, and visualizing scRNA seq data. In principle, interoperability exists between Seurat and Loom files. Seurocity also has powerful tools for integrating multiple datasets (https://satijalab.org/seurat/v3.0/integration.html). Further, there is an R-based Seurat wrapper for RNA velocity estimation using Velocyto outpute (https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/velocity.html). 

It is my opinion that Python-based scVelo is more efficient and user-friendly for velocity analysis, while I prefer Seurat's integration tools. **Seurocity** was created to permit easy initial processing, filtering, analysis, clustering, etc in Seurocty, followed by porting of dimensionality reductions and clustering to scVelo for velocity analysis.

### Usage information
