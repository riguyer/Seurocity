# Script to extract cell IDs, dimensionality reductions, PCA loadings, and cluster
#		identities from a Seurat object. Information is extracted as .csv files that
#		can be read into AnnData objets by Scanpy and Scvelo. This allows east porting
#		of Seurat processing to analysus using these alternate packages.
#
# Inputs: 	- rds file containing a Seurat objects
#						- idents.txt, a txt file containing entries in the "orig.ident" metadata
#							category to extract iformation (extracts for all if file is empty).
#							Values should be separated by commas with no spaces.
#						- reductions.txt, a txt file containing reductions to export. values
#							should be separated by commas with no spaces.
#						- append.txt, a txt file containing strings to append to unique cell
#							barcodes for each orig.ident group, in same order as orig.idents
#							are in idents.csv. Values should be separated by commas with no spaces.
#								# First row: prefix before barcode ("" or NULL if none)
#								# second row: postfix after barcode ("" or NULL if none)
#								# third row: single value, regex pattern to replace at end of barcode ("" or NULL if none)
#
#	Usage: Rscript extractSeurat.R [rds file] [idents.txt] [append.txt]
#
# This code (c) Richard A. Guyer, MD, PhD April 10, 2020
#
# Last updated: April 17, 2020 by Richard A. Guyer
#
# MGH Pediatric Surgical Research Laboratories
# PI: Allan Goldstein, MD
# riguyer@gmail.com rguyer@partners.org

# load required libraries
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(matrixStats))

# define local functions
#		mutate cell_ID metadata category as needed
MutateCellID <- function(seurat_obj, prefix, postfix, pattern) {
	if (! is.null(prefix) | prefix == "") {
		seurat_obj[["cell_ID"]] <- sapply(seurat_obj[["cell_ID"]], function(x) paste0(prefix,x))
	}
	if (! is.null(postfix) | postfix == "") {
		if (! is.null(pattern) | pattern == "") {
			seurat_obj[["cell_ID"]] <- sapply(seurat_obj[["cell_ID"]], function(x) str_replace(x, pattern,postfix))
		} else {
			seurat_obj[["cell_ID"]] <- sapply(seurat_obj[["cell_ID"]], function(x) paste0(x,postfix))
		}
	}
	return(seurat_obj)
}

#		export embeddings for a given dimensionality reduction
ExportEmbeddings <- function(seurat_obj, ident, reductions, row_names) {
	cat("Exporting dimenstionality reductions\n")
	for (r in reductions) {
		cat(paste0("Exporting reduction: ",r,"\n"))
		varname <- paste0(ident,"_seur_",r,"_emb")
		export_file <- paste0("./outputs/seurat_dat/",varname,".csv")
		assign(varname, seurat_obj[[r]]@cell.embeddings)
		write.csv(get(varname), file = export_file, row.names = row_names)
		cat(paste0("Saved as: ",export_file,"\n"))
		remove(list = c("varname","export_file"))
	}
}

#		export cluster identities
ExportClusters <- function(seurat_obj, ident, row_names) {
	cat("Exporting Seurat cluster identities\n")
	export_file <- paste0("./outputs/seurat_dat/",ident,"_seur_clust.csv")
	seur_clusters <- seurat_obj[["seurat_clusters"]]
	seur_clusters$seurat_clusters <- as.integer(as.vector(seur_clusters$seurat_clusters))
	write.csv(seur_clusters, file = export_file, row.names = row_names)
	cat(paste0("Saved as: ",export_file,"\n"))
}

#		export PCA loadings and variance data
ExportPCA <- function(seurat_obj) {
	cat("\nExporting PCA loadings and variance information\n")
	write.csv(seurat_obj[['pca']]@feature.loadings, file = "./outputs/seurat_dat/seur_pca_loadings.csv")

	# get and export variance and variance ratio data
	variance <- t(as.matrix((seurat_obj[['pca']]@stdev)^2))
	row.names(variance) <- "variance"
	variance_ratio <- as.matrix(variance/(sum(rowVars(GetAssayData(seurat_obj,
														assay = "integrated", slot = "scale.data")))))
	row.names(variance_ratio) <- "variance_ratio"
	var_data <- rbind(variance, variance_ratio)
	columns <- rep('PC',50); colnums <- c(1:50); columns <- paste(columns, colnums, sep="_")
	colnames(var_data) <- columns
	write.csv(var_data, file = "./outputs/seurat_dat/seur_pca_var.csv")

	# clean up
	remove(list = c("variance","variance_ratio","columns","colnums"))
}

###################################
# Execution script below this point
# 	create output directory, if does not exist
if (! dir.exists("./outputs/seurat_dat")) {
	cat("Outputs to be save in new directory: ./outputs/seurat_dat\n")
	cat("\n")
	dir.create("./outputs")
	dir.create("./outputs/seurat_dat")
} else {
	cat("Outputs to be written to existing directory: ./outputs/seurat_dat\n")
	cat("WARNING: Matching filenames will be overwritten\n")
	cat("\n")
}

# 	get arguments, point to inputs directory
args <- commandArgs(trailingOnly=TRUE)

# 	check whether appropraite arguments have been provided
#		this should not be a problem unless this is run as a standalone script
#		if run from importSeurat.py, expected file structure and arguments will be provided
if (length(args) < 4) {
	cat("WARNING: not all expected input files are present\n")
	cat("\n")
} else if (length(args) > 4) {
	cat("WARNING: unexpected arguments passed, will ignore any arguemnts beyond\n")
	cat("	the firsr three.\n")
	cat("\n")
}

if (! str_ends(args[1], ".rds")) {
	cat("ERROR: arg[1] must be the name of an .rds file containing a Seurat object\n")
	quit(status=1)
} else if (! (args[2] == "idents.txt") & (args[3] == "append.txt")) {
	cat("ERROR: arg[2] and arg[3] must be names of .txt files, as described in the README file\n")
	quit(status=1)
}

# 	redirect arguments to inputs subdirectory
args <- as.vector(sapply(args, function(x) paste0("./inputs/",x)))

rds_file <- args[1]
idents <- scan(args[2], what=character(), sep=",", quiet=TRUE)
reductions <- scan(args[3], what=character(), sep=",", quiet = TRUE)
append <- scan(args[4], what=character(), sep=",", quiet=TRUE)

NoI <- length(idents)
prefix <- append[1:NoI]
postfix <- append[(NoI+1):(NoI*2)]
pattern <- append[1+(NoI*2)]

# 	load Seurat data
cat(paste0("Loading seurat object from: ",rds_file,"\n"))
cat("\n")
cells <- readRDS(rds_file)

# ensure rds file contained a Seurat object
if (! class(cells)[1] == "Seurat") {
	cat("ERROR: input rds file does not contain a Seurat object\n")
	quit(status=1)
}

# 	set cell identities based on orig.ident metadata category
Idents(cells) <- "orig.ident"

curr_index <- 1
for (i in idents) {
	# subset seurat object by identity, mutate cell IDs as needed for use with
	#		Velocyto output, and extract reduced dimension embeddings and Seurat cluster
	#		identities for cells
	cat(paste0("**Subsetting out sample ",i," for export**\n"))
	curr_subset <- subset(cells, ident = i)
	curr_subset[["cell_ID"]] <- Cells(curr_subset)
	curr_subset <- MutateCellID(curr_subset, prefix[curr_index], postfix[curr_index], pattern)
	curr_rownames <- curr_subset[["cell_ID"]][,1]
	ExportEmbeddings(curr_subset, i, reductions, curr_rownames)
	ExportClusters(curr_subset, i, curr_rownames)

	# extract PCA loadings, variance, and variance ratio. These will be the same
	#		for all samples, so this only needs to be done once.
	if (curr_index == 1) {
		# get and export loadings
		ExportPCA(curr_subset)
	}
	cat("\n")

	curr_index <- curr_index + 1
}
