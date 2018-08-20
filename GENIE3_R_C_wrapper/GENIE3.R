#' @title GENIE3
#' 
#' @description \code{GENIE3} Infers a gene regulatory network (in the form of a weighted adjacency matrix) from expression data, using ensembles of regression trees.
#'
#' @param expr.matrix Expression matrix (genes x samples). Every row is a gene, every column is a sample. 
#' @param tree.method Tree-based method used. Must be either "RF" for Random Forests (default) or "ET" for Extra-Trees.
#' @param K Number of candidate regulators randomly selected at each tree node (for the determination of the best split). Must be either "sqrt" for the square root of the total number of candidate regulators (default), "all" for the total number of candidate regulators, or a stricly positive integer.
#' @param ntrees Number of trees in an ensemble for each target gene. Default: 1000.
#' @param regulators Subset of genes used as candidate regulators. Must be either a vector of indices, e.g. \code{c(1,5,6,7)}, or a vector of gene names, e.g. \code{c("at_12377", "at_10912")}. The default value NULL means that all the genes are used as candidate regulators.
#' @param ncores Number of cores to use for parallel computing. Default: 1.
#' @param verbose If set to TRUE, a feedback on the progress of the calculations is given. Default: FALSE.
#' @param seed Random number generator seed for replication of analyses. The default value NULL means that the seed is not reset.
#'
#' @return Weighted adjacency matrix of inferred network. Element w_ij (row i, column j) gives the importance of the link from regulatory gene i to target gene j. 
#' 
#' @examples
#' ## Generate fake expression matrix
#' expr.matrix <- matrix(sample(1:10, 100, replace=TRUE), nrow=20)
#' rownames(expr.matrix) <- paste("Gene", 1:20, sep="")
#' colnames(expr.matrix) <- paste("Sample", 1:5, sep="")
#'
#' ## Run GENIE3
#' weight.matrix <- GENIE3(expr.matrix, regulators=paste("Gene", 1:5, sep=""))
#' 
#' ## Get ranking of edges 
#' link.list <- get.link.list(weight.matrix)
#' head(link.list)
#' @export
GENIE3 <- function(expr.matrix, tree.method="RF", K="sqrt", ntrees=1000, regulators=NULL, ncores=1, verbose=FALSE, seed=NULL) {

	dyn.load("GENIE3.so")

	# check input arguments
	if (!is.matrix(expr.matrix) && !is.array(expr.matrix)) {
		stop("Parameter expr.matrix must be a two-dimensional matrix where each row corresponds to a gene and each column corresponds to a condition/sample.")
	}
	
	if (length(dim(expr.matrix)) != 2) {
		stop("Parameter expr.matrix must be a two-dimensional matrix where each row corresponds to a gene and each column corresponds to a condition/sample.")
	}
	
	if (is.null(rownames(expr.matrix))) {
		stop("expr.matrix must specify the names of the genes in rownames(expr.matrix).")
	}
	
	if (tree.method != "RF" && tree.method != "ET") {
		stop("Parameter tree.method must be \"RF\" (Random Forests) or \"ET\" (Extra-Trees).")
	}
	
	if (K != "sqrt" && K != "all" && !is.numeric(K)) {
		stop("Parameter K must be \"sqrt\", or \"all\", or a strictly positive integer.")
	}
	
	if (is.numeric(K) && K<1) {
		stop("Parameter K must be \"sqrt\", or \"all\", or a strictly positive integer.")
	}
	
	if (!is.numeric(ntrees) || ntrees<1) {
		stop("Parameter ntrees should be a stricly positive integer.")
	}
	
	if (!is.null(regulators)) {
		if (!is.vector(regulators)) {
			stop("Parameter regulators must be either a vector of indices or a vector of gene names.")
		}
		
		if (is.character(regulators) && length(intersect(regulators,rownames(expr.matrix))) == 0) {
			stop("The genes must contain at least one candidate regulator.")
		}
		
		if (is.numeric(regulators) && max(regulators) > dim(expr.matrix)[1]) {
			stop("At least one index in regulators exceeds the number of genes.")
		}
	}
	
	if (!is.numeric(ncores) || ncores<1) {
		stop("Parameter ncores should be a stricly positive integer.")
	}
	
	
	# set random number generator seed if seed is given
    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    # transpose expression matrix to (samples x genes)
    expr.matrix <- t(expr.matrix)
	
    # setup weight matrix
    num.samples <- dim(expr.matrix)[1]
    num.genes <- dim(expr.matrix)[2]
    gene.names <- colnames(expr.matrix)
    weight.matrix <- matrix(0.0, nrow=num.genes, ncol=num.genes)
    rownames(weight.matrix) <- gene.names
    colnames(weight.matrix) <- gene.names
	
    # get names of input genes
    if (is.null(regulators)) {
        input.gene.names <- gene.names
    } else {
        # input gene indices given as integers
        if (is.numeric(regulators)) {
            input.gene.names <- gene.names[regulators]
        # input gene indices given as names
        } else {
            input.gene.names <- regulators
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
	
	# tree method
	if (tree.method == 'RF') {
		RF_randomisation <- 1
		ET_randomisation <- 0
		bootstrap_sampling <- 1
	} else {
		RF_randomisation <- 0
		ET_randomisation <- 1
		bootstrap_sampling <- 0
	} 
	
	if (verbose) {
        cat(paste("Tree method: ", tree.method, "\nK: ", K,
	              "\nNumber of trees: ", ntrees, "\n\n",
                  sep=""))
        flush.console()
	}
    
    # compute importances for every target gene
   
	if (ncores==1) {
		# serial computing
		if (verbose) {
		    cat("Using 1 core.\n\n")
		    flush.console()
		}
		
	    for (target.gene.idx in seq(from=1, to=num.genes)) {

            if (verbose) {	
                cat(paste("Computing gene ", target.gene.idx, "/", num.genes, "\n", sep=""))
                flush.console()
			 }

	        target.gene.name <- gene.names[target.gene.idx]
	        # remove target gene from input genes
	        these.input.gene.names <- setdiff(input.gene.names, target.gene.name)
			num.input.genes <- length(these.input.gene.names)
		
	        x <- expr.matrix[,these.input.gene.names]
			y <- expr.matrix[,target.gene.name]
			
			# normalize output data
			y <- y / sd(y)

		    # set mtry
		    if (class(K) == "numeric") {
		        mtry <- K
		    } else if (K == "sqrt") {
		        mtry <- round(sqrt(num.input.genes))
		    } else {
		        mtry <- num.input.genes
		    } 
		
			# some default parameters 
			nmin <- 1
			permutation_importance <- 0
		
	        im <- .C("BuildTreeEns",as.integer(num.samples),as.integer(num.input.genes),
			          as.single(c(x)),as.single(c(y)),as.integer(nmin),
					  as.integer(ET_randomisation),as.integer(RF_randomisation),
					  as.integer(mtry),as.integer(ntrees),
					  as.integer(bootstrap_sampling),as.integer(permutation_importance),
					  as.double(vector("double",num.input.genes)))[[12]]
			
			# some variable importances might be slighly negative due to some rounding error
			im[im<0] <- 0
					  
	        weight.matrix[these.input.gene.names, target.gene.name] <- im
	    }
	} else {
		# parallel computing
	    library(doRNG); library(doParallel); registerDoParallel(); options(cores=ncores)
		
		if (verbose) {
		    message(paste("\nUsing", getDoParWorkers(), "cores."))
		}
		
	    weight.matrix.reg <- foreach(target.gene.name=gene.names, .combine=cbind) %dorng% 
	    {
	        # remove target gene from input genes
	        these.input.gene.names <- setdiff(input.gene.names, target.gene.name)
			num.input.genes <- length(these.input.gene.names)
		
	        x <- expr.matrix[,these.input.gene.names]
			y <- expr.matrix[,target.gene.name]
			
			# normalize output data
			y <- y / sd(y)

		    # set mtry
		    if (class(K) == "numeric") {
		        mtry <- K
		    } else if (K == "sqrt") {
		        mtry <- round(sqrt(num.input.genes))
		    } else {
		        mtry <- num.input.genes
		    } 
			
			# some default parameters 
			nmin <- 1
			permutation_importance <- 0
		
	        im <- .C("BuildTreeEns",as.integer(num.samples),as.integer(num.input.genes),
			          as.single(c(x)),as.single(c(y)),as.integer(nmin),
					  as.integer(ET_randomisation),as.integer(RF_randomisation),
					  as.integer(mtry),as.integer(ntrees),
					  as.integer(bootstrap_sampling),as.integer(permutation_importance),
					  as.double(vector("double",num.input.genes)))[[12]]
					  
		  	# some variable importances might be slighly negative due to some rounding error
		  	im[im<0] <- 0
					  			  
			c(setNames(0, target.gene.name), setNames(im, these.input.gene.names))[input.gene.names]
	    }
	    attr(weight.matrix.reg, "rng") <- NULL
	    weight.matrix[input.gene.names,] <- weight.matrix.reg
	}
    return(weight.matrix / num.samples)
}       



#' @title get.link.list
#' 
#' @description \code{get.link.list} Converts the weight matrix, as returned by \code{\link{GENIE3}}, to a sorted list of regulatory links (most likely links first).
#' 
#' @param weight.matrix Weighted adjacency matrix as returned by \code{\link{GENIE3}}.
#' @param report.max Maximum number of links to report. The default value NULL means that all the links are reported.
#' @param threshold Only links with a weight above the threshold are reported. Default: threshold = 0, i.e. all the links are reported.
#' 
#' @return List of regulatory links in a data frame. Each line of the data frame corresponds to a link. The first column is the regulatory gene, the second column is the target gene, and the third column is the weight of the link.
#'
#' @seealso \code{\link{GENIE3}}
#'
#' @examples
#' ## Generate fake expression matrix
#' expr.matrix <- matrix(sample(1:10, 100, replace=TRUE), nrow=20)
#' rownames(expr.matrix) <- paste("Gene", 1:20, sep="")
#' colnames(expr.matrix) <- paste("Sample", 1:5, sep="")
#'
#' ## Run GENIE3
#' weight.matrix <- GENIE3(expr.matrix, regulators=paste("Gene", 1:5, sep=""))
#' 
#' ## Get ranking of edges 
#' link.list <- get.link.list(weight.matrix)
#' head(link.list)
#' @export
get.link.list <- function(weight.matrix, report.max=NULL, threshold=0) {
    if(!is.numeric(threshold)) {
    	stop("threshold must be a number.")
    } 
	
	library(reshape2)
	
	# Only process weights off-diagonal
	diag(weight.matrix) <- NA
    link.list <- melt(weight.matrix, na.rm=TRUE)
    colnames(link.list) <- c("regulatory.gene", "target.gene", "weight")
    link.list <- link.list[link.list$weight>=threshold,]
    link.list <- link.list[order(link.list$weight, decreasing=TRUE),]
  
    if(!is.null(report.max)) {
    	link.list <- link.list[1:min(nrow(link.list), report.max),]
    } 
  
    rownames(link.list) <- NULL
  
    return(link.list)
}




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
