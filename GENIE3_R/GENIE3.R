# R implementation of GEne Network Inference with Ensemble of Trees (GENIE3)
#
# See the accompanying README.txt for usage and references
#
# Modification history:
#    19 July 2010:   - first version
#    13 August 2010: - fixed formatting of link list table
#                    - fixed warnings when K="all"
#                    - added choice of importance measure
#                    - added possibility to fix random number generator seed
#
# Modifications by Van Anh:
#    13 July 2016:   - permutation importances are no longer automatically computed
#                    - changed the data normalization:
#                        * For impurity importances: only the output gene is normalized (by dividing by the standard deviation)
#                        * For permutation importances: no pre-normalization
#    8 August 2016:  - changed the default value of the nodesize parameter (parameter of the randomForest function). Now, nodesize=1 by default, i.e. trees are fully developed.
#
# Authors: Alexandre Irrthum, Van Anh Huynh-Thu (vahuynh@ulg.ac.be)



# read.expr.matrix: read expression matrix from file.
#
# This function will detect gene and/or sample names automatically
# in first row and/or first column if they are provided. Gene and
# sample names must contain alphabetic characters if provided.
#
# Parameters (required):
#    -- filename: name of or full path to expression matrix file
#    -- form: format of the expression matrix on file. Must be one of
#             "rows.are.genes" for a (genes x samples) matrix
#             "rows.are.samples" for a (samples x genes) matrix
#
# Parameters (optional):
#    -- sep: field separator in expression matrix file
#            By default, fields will be splitted on tabs or spaces
#    -- default.gene.label: string to prepend to gene indices to form
#                           default gene names (e.g. "gene_" -> "gene_1"...)
#    -- default.sample.label: string to prepend to sample indices to form
#                             default sample names (e.g. "sample_" -> "sample_1"...)
#
# Returns:
#    expression matrix in format (genes x samples), i.e. every row is a gene
#
read.expr.matrix <- function(filename, form="", sep="", default.gene.label="gene_", default.sample.label="sample_") {
    # Report when form is not correctly set
    if (form != "rows.are.genes" && form != "rows.are.samples") {
        stop("Parameter form must be set to \"rows.are.genes\" or \"rows.are.samples\"")
    }
    # read data
    m <- read.table(filename, sep=sep, as.is=TRUE)
    has.row.names <- FALSE
    has.col.names <- FALSE
    # have row and column names been recognized (cell (1,1) is empty in file) ?
    if (colnames(m)[1] != "V1" && rownames(m)[1] != "1") {
        has.row.names <- TRUE
        has.col.names <- TRUE
    }
    # is first column alphanumeric ?
    if (all(grepl("[[:alpha:]]", m[,1]))) {
        # check duplicated names
        if (any(duplicated(m[,1]))) {
            stop("Duplicated names in first column\n")
        }
        rownames(m) <- m[,1]
        m <- m[,-1]
        has.row.names <- TRUE
    }
    # is first row alphanumeric ?
    if (all(grepl("[[:alpha:]]", m[1,]))) {
        # check duplicated names
        if (any(duplicated(m[1,]))) {
            stop("Duplicated names in first row\n")
        }
        colnames(m) <- m[1,]
        m <- m[-1,]
        has.col.names <- TRUE
    }
    # coerce matrix data to numeric
    col.names <- colnames(m)
    row.names <- rownames(m)
    m <- as.matrix(m)
    m <- apply(m, 2, function(x) { as.numeric(x) })
    colnames(m) <- col.names
    rownames(m) <- row.names
    num.rows <- dim(m)[1]
    num.cols <- dim(m)[2]
    # fill in default gene names in rows if needed
    if (!has.row.names && form=="rows.are.genes") {
        rownames(m) <- paste(default.gene.label, seq(from=1, to=num.rows), sep="")
    }
    # fill in default sample names in rows if needed
    if (!has.row.names && form=="rows.are.samples") {
        rownames(m) <- paste(default.sample.label, seq(from=1, to=num.rows), sep="")
    }
    # fill in default sample names in columns if needed
    if (!has.col.names && form=="rows.are.genes") {
        colnames(m) <- paste(default.sample.label, seq(from=1, to=num.cols), sep="")
    }
    # fill in default gene names in columns if needed
    if (!has.col.names && form=="rows.are.samples") {
        colnames(m) <- paste(default.gene.label, seq(from=1, to=num.cols), sep="")
    }
    # transpose matrix to (genes x samples) if needed
    if (form == "rows.are.samples") m <- t(m)
    return(m)
}

# GENIE3: compute weighted adjacency matrix of inferred network
# 
# Parameters (required):
#    -- expr.matrix: gene expression matrix as returned by read.expr.matrix
#
# Parameters (optional):
#    -- K: choice of number of input genes randomly selected as candidates at each node
#          Must be one of
#          - "sqrt" for square root of total number of input genes (default)
#          - "all" for total number of input genes (minus 1)
#          -  an integer 
#    -- nb.trees: number of trees in ensemble for each target gene (default 1000)
#    -- input.idx: subset of genes used as input genes (default all genes)
#                  Must be either
#                  - a vector of indices, e.g. c(1,5,6,7), or
#                  - a vector of gene names, e.g. c("at_12377", "at_10912")
#    -- importance.measure: Type of variable importance measure
#          Must be one of
#          - "IncNodePurity" for importance measure based on decrease of residual
#             sum of squares (default)
#          - "%IncMSE" for importance measure obtained by permutation of OOB data
#    -- seed: random number generator seed for replication of analyses
#          (default NULL means the seed is not reset)
#    -- trace: index of currently computed gene is reported (default TRUE)
#    -- All additional parameters are passed to the randomForest function
#       (see randomForest manual for more info)
#
# Returns:
#    weighted adjacency matrix of inferred network.
#    element w_ij (row i, column j) gives the importance of the link
#    from regulatory gene i to target gene j
#
GENIE3 <- function(expr.matrix, K="sqrt", nb.trees=1000, input.idx=NULL, importance.measure="IncNodePurity", seed=NULL, trace=TRUE, ...) {
    # set random number generator seed if seed is given
    if (!is.null(seed)) {
        set.seed(seed)
    }
    # to be nice, report when parameter importance.measure is not correctly spelled
    if (importance.measure != "IncNodePurity" && importance.measure != "%IncMSE") {
        stop("Parameter importance.measure must be \"IncNodePurity\" or \"%IncMSE\"")
    }
	# Check if nodesize parameter is in the input arguments
	args <- list(...)
	nodesize.in.args <- "nodesize" %in% names(args)
    # transpose expression matrix to (samples x genes)
    expr.matrix <- t(expr.matrix)
    # setup weight matrix
    num.samples <- dim(expr.matrix)[1]
    num.genes <- dim(expr.matrix)[2]
    gene.names <- colnames(expr.matrix)
    weight.matrix <- matrix(0.0, nrow=num.genes, ncol=num.genes)
    rownames(weight.matrix) <- gene.names
    colnames(weight.matrix) <- gene.names
    # get number of input genes, names of input genes
    if (is.null(input.idx)) {
        input.gene.names <- gene.names
    } else {
        # input gene indices given as integers
        if (is.numeric(input.idx)) {
            input.gene.names <- gene.names[input.idx]
        # input gene indices given as names
        } else {
            input.gene.names <- input.idx
            # for security, abort if some input gene name is not in gene names
            missing.gene.names <- setdiff(input.gene.names, gene.names)
            if (length(missing.gene.names) != 0) {
                for (missing.gene.name in missing.gene.names) {
                    cat(paste("Gene ", missing.gene.name,
                              " was not in the expression matrix\n", sep=""))
                }
                stop("Aborting computation")
            }
        }
    }
    # compute importances for every target gene
    for (target.gene.idx in seq(from=1, to=num.genes)) {
        if (trace) {
            cat(paste("Computing gene ", target.gene.idx, "/", num.genes, "\n", sep=""))
            flush.console()
        }
        target.gene.name <- gene.names[target.gene.idx]
        # remove target gene from input genes
        these.input.gene.names <- setdiff(input.gene.names, target.gene.name)
		num.input.genes <- length(these.input.gene.names)
        x <- expr.matrix[,these.input.gene.names]
		y <- expr.matrix[,target.gene.name]
	    # set mtry
	    if (class(K) == "numeric") {
	        mtry <- K
	    } else if (K == "sqrt") {
	        mtry <- round(sqrt(num.input.genes))
	    } else if (K == "all") {
	        mtry <- num.input.genes
	    } else {
	        stop("Parameter K must be \"sqrt\", or \"all\", or an integer")
	    }
	    if (trace) {
	        cat(paste("K = ", mtry,", ", nb.trees, " trees\n\n",
	                  sep=""))
	        flush.console()
	    }
        if (importance.measure == "%IncMSE") {
			if (nodesize.in.args) {
				rf <- randomForest(x, y, mtry=mtry, ntree=nb.trees, importance=TRUE, ...)
			} else {
				# By default, grow fully developed trees
				rf <- randomForest(x, y, mtry=mtry, ntree=nb.trees, importance=TRUE, nodesize=1, ...)
			}
            
        } else {
			# Normalize output
			y <- y / sd(y)
			if (nodesize.in.args) {
				rf <- randomForest(x, y, mtry=mtry, ntree=nb.trees, importance=FALSE, ...)
			} else {
				# By default, grow fully developed trees
				rf <- randomForest(x, y, mtry=mtry, ntree=nb.trees, importance=FALSE, nodesize=1, ...)
			}
        }
        im <- importance(rf)[,importance.measure]
        im.names <- names(im)
        weight.matrix[im.names, target.gene.name] <- im
    }
    return(weight.matrix / num.samples)
}       

# get.link.list: get sorted list of regulatory links (most likely link first)
#
# Parameters (required):
#    -- weight.matrix: weighted adjacency matrix as returned by GENIE3
#
# Parameters (optional):
#    -- report.max: maximum number of links to report (default all links)
#
# Returns:
#    list of links in data frame. Each line of the data frame has format
#    regulatory_gene target_gene importance_score
#
get.link.list <- function(weight.matrix, report.max=NULL) {
    # set negative weights (for permutation of OOB importance) to 0.0 
    # weight.matrix[weight.matrix < 0.0] <- 0.0
    num.genes <- dim(weight.matrix)[1]
    genes <- colnames(weight.matrix)
    matrix.length <- length(weight.matrix)
    list.length <- num.genes * (num.genes - 1)
    if (!is.null(report.max) && report.max < list.length) {
        list.length <- report.max
    }
    # setup link list
    link.list <- data.frame(from.gene=rep("", list.length),
                            to.gene=rep("", list.length),
                            im=rep(0.0, list.length),
                            stringsAsFactors=FALSE)
    sorted.indices <- order(weight.matrix, decreasing=TRUE)
    # fill in link list
    index.number <- 1
    link.number <- 1
    while (index.number <= matrix.length && link.number <= list.length) {
        i <- sorted.indices[index.number]
        im <- weight.matrix[i]
        row.col <- lin.to.square(i, num.genes)
        row <- row.col[1]
        col <- row.col[2]
        # Only process weights off-diagonal
        if (row != col) {
            from.gene <- genes[row]
            to.gene <- genes[col]
            link.list[link.number,] <- list(from.gene, to.gene, im)
            link.number <- link.number + 1
        }
        index.number <-index.number + 1
    }
    return(link.list)
}

# load required packages
tryCatch( suppressMessages(library(randomForest)),
          error=function(e) { cat("Error: package randomForest must be installed\n");
                                cat("Use install.packages(\"randomForest\")\n") })

# utility function to convert linear index to (row,col) index for matrix
lin.to.square <- function(i, nrow) {
    col <- ((i - 1) %/% nrow) + 1
    row <- ((i - 1) %% nrow) + 1
    return(c(row, col))
}
