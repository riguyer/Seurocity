# import required modules
try: 
    import os, sys, getopt
    import numpy as np
    import scvelo as scv
    import loompy as lp
except ImportError as imperror:
    print(imperror)
    print("Unable to load required modules, exiting program")
    print("")
    sys.exit(1)

# define functions
def a_notin_b(list1, list2):
    """
    use numpy arrays to quickly check if elements of list1 are in list2

    requires numpy

    Parameters
    ----------
    list1 : list
        list a, to be check against list b.
    list2 : list
        list b.

    Returns
    -------
    list of items in list a but not in list b.
    """
    a = np.asarray(list1)
    b = np.asarray(list2)

    return a[np.vectorize(lambda x :  x not in b)(a)].tolist()

def get_duplicates(array):
    """
    return a list of duplicate items in a numpy array
    code from: https://github.com/theislab/scvelo/issues/36#issuecomment-508027076
    
    Parameters
    ----------
    array : numpy array
        numpy array to check for duplicate items
    
    Returns
    -------
    numpy array
        a new numpy array containing one instance of every item duplicated in the input array.
    
    """
    from collections import Counter
    return np.array([item for (item, count) in Counter(array).items() if count > 1])

def files_exist(files_list):
    """
    check whether a list of files exists in specified locations

    Parameters
    ----------
    files_list : list
        list of files with paths relative to the current directory.

    Returns
    -------
    bool : bool
        True if all exist, False if any do not exist
    """
    for file in files_list:
        if not os.path.exists(file):
            return(False)
    return(True)

def run_rscript(path_to_script, arguments, no_return=True, kill_if_nonzero=True, rscript_path=None):
    """
    run an R script from within Python script

    Parameters
    ----------
    path_to_script : str
        relative or absolute path to R script.
    arguments : list
        list of arguments to be passed to the R script.
    no_return : bool, optional
        run script without returing any output. The default is True.
    kill_if_nonzero : bool, optional
        kill Python script if R script exits with non-zero status. The default is True.
    rscript_path : str, optional
        runs Rscript from a specific path if provided, otherwise runs system's default instance of Rscript. Default is None.

    Returns
    -------
    Output of R script if no_return=False, otherwise none
    """
    import subprocess as sub

    if rscript_path==None:
        command = ['Rscript', path_to_script] + arguments
    else:
        command = [rscript_path + 'Rscript', path_to_script] + arguments

    if no_return:
        try:
            sub.check_call(command)
        except sub.CalledProcessError:
            print("ERROR: " + path_to_script + " returned a non-zero exit status")
            if kill_if_nonzero:
                sys.exit(1)
    else:
        try:
            return(sub.check_output(command, universal_newlines=True))
        except sub.CalledProcessError:
            print("ERROR: " + path_to_script + " returned a non-zero exit status")
            if kill_if_nonzero:
                sys.exit(1)

def load_samples(sample_ids):
    """
    Load Velocyto output loom files and seurat data

    Parameters
    ----------
    sample_ids : list
        list of sample IDs.

    Returns
    -------
    samples : dict
        dictionary of dictionaries, {s : dict} for s in sample_ids
            Each subdictionary contains:
                AnnData object (key: 'main') loaded from loom file
                AnnData object (key: from filename) for csv files with seurat data
    """
    samples = {}
    for s in sample_ids:
        samples[s] = {}

        filename = "./inputs/" + s + ".loom"
        comment = "\nLoading sample " + s + " from " + filename
        print(comment)
        samples[s]['main'] = scv.read(filename, cache = False)

        csv_files = filter(lambda x: x.startswith(s) & x.endswith(".csv"), 
                           os.listdir("./outputs/seurat_dat/"))
        for csvf in csv_files:
            varname = os.path.splitext(csvf)[0].casefold()[len(s)+1:]
            samples[s][varname] = scv.read("./outputs/seurat_dat/" + csvf, cache = False)
    return(samples)

def get_reductions():
    """
    Returns list of dimensionality reductions exported from seurat, based on 
        filenames in ./outputs/seurat_dat. This function requires the expected
        file struture for Seurocity.py and output files in the same format 
        used by extractSeurat.R.

    Ensures Seurocity.py will work with recutions successfully exported
        from seurat by extractSeurat.R, not those requested in reductions.txt.

    Returns
    -------
    reductions : list
        list of dimensionality reductions exported from seurat

    """
    reductions = []
    embs = filter(lambda x: x.endswith("_emb.csv"), os.listdir("./outputs/seurat_dat/"))
    for emb in embs:
        emb = emb.split("_")[2]
        reductions.append(emb)
    reductions = list(set(reductions))
    return(reductions)

def get_ids():
    """
    Returns list of sample IDs exported from seurat, based on filenames in 
        ./outputs/seurat_dat. This function requires the expected file struture 
        for Seurocity.py and output files in the same format used by 
        extractSeurat.R
    
    Ensures Seurocity.py will work with sample IDs successfully exported
        from seurat by extractSeurat.R, not those requested in idents.txt.

    Returns
    -------
    ids : list
        list of sample IDs exported from seurat
    """
    ids = []
    embs = filter(lambda x: x.endswith("_emb.csv"), os.listdir("./outputs/seurat_dat/"))
    for emb in embs:
        emb = emb.split("_")[0]
        ids.append(emb)
    ids = list(set(ids))
    return(ids)   

def load_pca_data(dictionary):
    """
    Loads PCA loadings and variance data output by extractSeurat.R

    Parameters
    ----------
    dictionary : dict
        dictionary

    Returns
    -------
    dictionary : dict
        dictionary, with keys added for PCA loadings and variance 
    """
    try:
        comment = "\nImporting PCA loadings and variance data"
        print(comment)
        dictionary['pca_data'] = {}
        dictionary['pca_data']['loadings'] = scv.read("./outputs/seurat_dat/seur_pca_loadings.csv")
        dictionary['pca_data']['variance'] = scv.read("./outputs/seurat_dat/seur_pca_var.csv")
    except OSError as err:
        print("ERROR: " + err)
        print("Proceeding without PCA variance or loadings")
    return(dictionary)

def import_seur_data(dictionary, sample_ids, reductions):
    """
    Import seurat data into main AnnData object for each sample

    Parameters
    ----------
    dictionary : dict
        dictionary output by sequential use of load_samples() and load_pca_data()
    sample_ids : list
        list of sample IDs from get_ids()
    reductions : list
        list of reductions from get_reductions()

    Returns
    -------
    dictionary : dict
        dictionary, with seurat data loaded into main AnnData object
    """
    for s in sample_ids:
        # filter cells and genes based on Seurat processing
        comment = "\n- Processing main AnnData based on imported Seurat data for sample " + s
        print(comment)

        barcodes = dictionary[s]['seur_umap_emb'].obs.index.values.tolist()
        dictionary[s]['main'] = dictionary[s]['main'][barcodes,:]

        genes = dictionary['pca_data']['loadings'].obs.index.tolist()
        loom_genes = dictionary[s]['main'].var.index.tolist()
        discard_genes = a_notin_b(genes, loom_genes) + get_duplicates(dictionary[s]['main'].var_names).tolist()
        genes = a_notin_b(genes, discard_genes)
        dictionary['pca_data']['loadings'] = dictionary['pca_data']['loadings'][genes,:]
        dictionary[s]['main'] = dictionary[s]['main'][:,genes]

        comment = "--- Loading Seurat clusters into main AnnData"
        print(comment)
        dictionary[s]['seur_clust'].X = dictionary[s]['seur_clust'].X.astype(int).astype(str) # converting to int and then to string prevents decimals on cluster IDs
        dictionary[s]['main'].obs['seurat_clust'] = dictionary[s]['seur_clust'].X.tolist()

        # load reduction embeddings
        for r in reductions:
            comment = "--- Loading " + r + " embeddings into main AnnData"
            print(comment)
            obsm_name = "X_" + r.casefold()
            r_location = "seur_" + r + "_emb"
            dictionary[s]['main'].obsm[obsm_name] = dictionary[s][r_location].X

        # load pca variance data and feature loadings from Seurat into
        #    main AnnData object.
        comment = "--- Loading PCA loadings and variance data into main AnnData"
        print(comment)       
        dictionary[s]['main'].uns['pca'] = {}
        dictionary[s]['main'].uns['pca'][dictionary['pca_data']['variance'].obs.index[0]] = dictionary['pca_data']['variance'].X[0]
        dictionary[s]['main'].uns['pca'][dictionary['pca_data']['variance'].obs.index[1]] = dictionary['pca_data']['variance'].X[1]
        
        dictionary[s]['main'].varm['PCs'] = dictionary['pca_data']['loadings'].X
    
    return(dictionary)

def same_genes(dictionary, sample_ids): 
    """
    update main AnnData object for list of samples so all contain the same list
        of genes

    Parameters
    ----------
    dictionary : dict
        dictionary containing AnnData objects
    sample_ids : list
        list of sample IDs to compare against one another.

    Returns
    -------
    dictionary : dict
    """
    # find genes common to all samples
    genes = []
    for s in sample_ids:
        if len(genes)==0:
            genes = dictionary[s]['main'].var.index.tolist()
            continue
        else:
            sample_genes = dictionary[s]['main'].var.index.tolist()
            genes = np.intersect1d(np.asarray(genes),np.asarray(sample_genes)).tolist()
    
    # subset each sample to contain only common genes
    for s in sample_ids:
        dictionary[s]['main'] = dictionary[s]['main'][:,genes]
    
    # return updated dictionary
    return(dictionary)

def arg_error():
    print(""" 
    SEUROCITY v1.0.0
    
    Program for porting datasets integrated and processed by Seurat to scVelo. Allows
        intergration, clustering, and filtering of datasets in Seurat, followed by
        analysis of Velocyto output with scVelo.
    
    This code is intended to be a command line utility to provide a standardized analysis
       pipeline for use in the Goldstein Lab. Others are free to use this software free
       of charge but without any guarantee or warranty.
        
    (c) Richard A. Guyer, MD, PhD April 1, 2020
    MGH Pediatric Surgical Research Laboratories
    PI: Allan Goldstein, MD
    riguyer@gmail.com rguyer@partners.org
    
    Created on Wed Apr 1 20:45:00 2020
    Last modified: Wed Apr 15 12:53:38 2020
    Built and tested in Python 3.7.3 on MacBook Air with:
       8 GB 1600 MHz DDR3
       2.2 GHz Intel Core i7
    Using Spyder 4.0.1 development environment
    
    R script built and tested in R 3.6.1
        R packages used:
            stringr 1.4.0 
            dplyr 0.8.3  
            Seurat 3.1.2 
    
    Distributed under MIT license without warranty
    
    Usage: python Seurocity.py [ARGS]
    -h for help documentation
    -l to view license
          """)

def display_help():
    print("""
    SEUROCITY v1.0.0
    
    Program for porting datasets integrated and processed by Seurat to 
        scVelo. Allows intergration, clustering, and filtering of 
        datasets in Seurat, followed by analysis of Velocyto output 
        with scVelo.
    
    This code is intended to be a command line utility to provide a 
        standardized analysis pipeline for use in the Goldstein Lab. 
        Others are invited to use this software free of charge, but 
        neither guarantee nor warranty is offered.
        
    (c) Richard A. Guyer, MD, PhD April 1, 2020
    MGH Pediatric Surgical Research Laboratories
    PI: Allan Goldstein, MD
    riguyer@gmail.com rguyer@partners.org
    
    Created on Wed Apr 1 20:45:00 2020
    Last modified: Wed Apr 15 12:53:38 2020
    Built and tested in Python 3.7.3 on MacBook Air with:
       8 GB 1600 MHz DDR3
       2.2 GHz Intel Core i7
    Using Spyder 4.0.1 development environment
    
    R script built and tested in R 3.6.1
        R packages used:
            stringr 1.4.0 
            dplyr 0.8.3  
            Seurat 3.1.2 
    
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
        
    Outputs will be generated with the following structure:
        
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
    """)

def display_license():
    print("""
    MIT License
    
    Copyright (c) 2020 Richard A. Guyer, MD, PhD
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
          """)

# main
def main(argv):
    print("SEUROCITY v1.0.0, (c) 2020 Richard A. Guyer, MD, PhD\n")
    
    # default for rscript_dir to pass to run_rscript function
    rscript_dir = None
    input_file = "input.rds"
    
    # handle arguments to adjust default settings
    try:
        opts, args = getopt.getopt(argv,"hlr:w:")
    except getopt.GetoptError:
        arg_error()
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-h':
            display_help()
            sys.exit(0)
        elif opt in ("-l"):
            display_license()
            sys.exit(0)
        elif opt in ("-r"):
            rscript_dir = arg
        elif opt in ("-w"):
            os.chdir(arg)
        elif opt in ("-i"):
            input_file = arg
    
    working_dir = os.getcwd() + "/"
    input_dir = working_dir + "inputs/"
    output_dir = working_dir + "outputs/"

	# check for files required by extractSeurat.R, run if all are present
    required_files_for_R = [input_dir + input_file,
                            input_dir + "idents.txt",
                            input_dir + "reductions.txt",
                            input_dir + "append.txt", 
                            working_dir + "extractSeurat.R"]
    if not files_exist(required_files_for_R):
        print("ERROR: Critical files not found in expected locations")
        print("Please ensure proper input file structure")
        print("For help: python Seurocity.py -h")
        print("")
        sys.exit(1)
    else:
        run_rscript(working_dir + "extractSeurat.R", [input_file,"idents.txt","reductions.txt",
                                                      "append.txt"], rscript_path=rscript_dir)

	# get sample IDs and reductions output by Rscript
    sample_ids = get_ids()
    reductions = get_reductions()
    
    # check whether expected loom files exist
    expected_looms = [input_dir + ident + ".loom" for ident in sample_ids]
    if not files_exist(expected_looms):
        print("ERROR: Expected loom files not found in ./inputs")
        print("Please ensure proper input file structure")
        print("For help: python Seurocity.py -h")
        print("")
        sys.exit(1)  
    else:
        print("\nLoading files and processing AnnData objects, this may take a few minutes")
        samples = load_samples(sample_ids)
        samples = load_pca_data(samples)
        samples = import_seur_data(samples, sample_ids, reductions)
        
        # ensure every sample has the same list of genes
        if len(sample_ids) > 1:
            samples = same_genes(samples, sample_ids)
    
        # save main AnnData object for each sample as a loom file 
        comment = "\nSaving main AnnData for each sample as loom files"
        print(comment)
        if os.path.exists(output_dir + "proc_loom"):
            comment = "- WARNING: ./outputs/proc_loom/ exists, files may be overwritten"
            print(comment)
        else:
            os.mkdir(output_dir + "proc_loom")
            
        for s in sample_ids:
            savename = output_dir + "proc_loom/" + s + "_proc.loom"
            comment = "- Saving sample " + s + " to: " + savename
            print(comment)
            samples[s]['main'].varm['PCs'] = np.asarray(samples[s]['main'].varm['PCs']) # currenty is an ArrayView, need to make into numpy array
            samples[s]['main'].write_loom(savename, write_obsm_varm = True)
        
        # remove samples to clean up memory
        del(samples)
        
        # generate combined loom file and import pca data
        #   generate combined file
        comment = "\nGenerating combined loom file with PCA data loaded"
        if os.path.exists(output_dir + "comb_loom"):
            comment = "- WARNING: ./outputs/comb_loom/ exists, combined.loom will be overwritten if it already exists"
            print(comment)
        else:
            os.mkdir(output_dir + "comb_loom")       
            
        processed_files = os.listdir(output_dir + "proc_loom/")
        processed_files = [output_dir + "proc_loom/" + p for p in processed_files]
        lp.combine(processed_files, output_dir + "comb_loom/combined.loom")
        
        #   load combined loom and pca data files
        combined = scv.read(output_dir + "comb_loom/combined.loom", cache = False)
        pca_var = scv.read(output_dir + "seurat_dat/seur_pca_var.csv")
        pca_load = scv.read(output_dir + "seurat_dat/seur_pca_loadings.csv")
        
        #   variance data 
        combined.uns['pca'] = {}
        combined.uns['pca'][pca_var.obs.index[0]] = pca_var.X[0]
        combined.uns['pca'][pca_var.obs.index[1]] = pca_var.X[1]
        
        #   pca loadings
        genes = combined.var.index.tolist()
        combined.varm['PCs'] = np.asarray(pca_load[genes,:].X)
        
        #   save combined, now containing pca loadings and variance data
        combined.write_loom(output_dir + "comb_loom/combined.loom", write_obsm_varm = True)

if __name__ == "__main__":
    main(sys.argv[1:])
