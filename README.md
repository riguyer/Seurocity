# Seurocity v1.0.0
*Easy velocity analysis using scVelo with data preprocessed in Seurat*

### Background
Velocyto (https://velocyto.org/) is a powerful tool for studying dynamic gene expression changes within single cell datasets. The Velocyto program is a command line tool that processes raw single cell sequencing data. Velocyto can accept data output by 10X Genetics Cell Ranger pipeline, as well as other common single cell analysis tools, such as SmartSeq2. Advanced users can also provide the required inputs files manually. The program outputs a Loom file (https://linnarssonlab.org/loompy/) containing spliced and unspliced gene counts for every cell. Multiple tools exist for analyzing and processing the data output by Velocyto, including velocyto.R and scVelo (https://scvelo.readthedocs.io/index.html). 

Seurat (https://satijalab.org/seurat/) is a powerful and widely-used R-based backaged for processing, analyzing, and visualizing scRNA seq data. In principle, interoperability exists between Seurat and Loom files. Seurocity also has powerful tools for integrating multiple datasets (https://satijalab.org/seurat/v3.0/integration.html). Further, there is an R-based Seurat wrapper for RNA velocity estimation using Velocyto outpute (https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/velocity.html). 

It is my opinion that Python-based scVelo is more efficient and user-friendly for velocity analysis, while I prefer Seurat's integration tools. **Seurocity** was created to permit easy initial processing, filtering, analysis, clustering, etc in Seurocty, followed by porting of dimensionality reductions and clustering to scVelo for velocity analysis.

### Usage information
    Distributed under MIT license without warranty
    
    Usage: python Seurocity.py [Options]
    
    Options: 
        -h display this help message
        -l display license
        -r rscript_dir : str
             path to Rscript instance to use, optional
        -w working_dir : str
            path to working directory, optional, default is ./
        -i input_file : str
            name of input rds file, default is 
    
    Expected input file structure in working directory:
        
    .
    ├─── inputs
    |    ├── idents.txt ; file with sample IDs
    |    ├── append.txt ; file with information for porting cell IDs
    |    ├── reductions.txt ; file with reductions to port from Seurat
    |    ├── input.rds ; rds-formatted file with saved Seurat object
    |    └── SAMPLE-ID.loom ; Velocyto-output loom (one for each sample)
    └─── extractSeurat.R
        
    The program will fail if these input files are not provided with 
        the indicated names and locations
        
    Outputs will be generated with the following structure inside the working directory:
        
    .
    └─── outputs
         ├── seurat_dat 
         |   ├── seur_pca_loadings.csv
         |   ├── seur_pca_var.csv 
         |   └── SAMPLE-ID_seur_REDUCTION_emb.csv ; multiple
         ├── proc_loom 
         |   └── SAMPLE-ID_proc.loom ; one for each sample
         └── comb_loom
             └── combined.loom ; a single loom file 
    
    If these output directories do not exist they will be created. If 
        the directories do exist they will be used. However, be aware 
        existing files will be overwritten if they have names matching 
        the output format.

### Development information
    (c) Richard A. Guyer, MD, PhD April 1, 2020
    MGH Pediatric Surgical Research Laboratories
    PI: Allan Goldstein, MD
    riguyer@gmail.com rguyer@partners.org
    
    Created on Wed Apr 1 20:45:00 2020
    Last modified: Wed Apr 15 12:53:38 2020
    Built and tested in Python 3.7.3 on MacBook Air with:
       8 GB 1600 MHz DDR3
       2.2 GHz Intel Core i7
    Using Spyder 4.0.1 development environment (Python code) and Atom 1.45.0 text editor (R code)
    
    R script built and tested in R 3.6.1
        R packages used:
            stringr 1.4.0 
            dplyr 0.8.3  
            Seurat 3.1.2

