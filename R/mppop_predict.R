#' Predict genetic variance and genetic correlations in multi-parent populations
#' using a deterministic equation.
#'
#' @description
#' Predicts the genotypic mean, genetic variance, and
#' usefulness criterion (superior progeny mean) in a set of multi-parent populations
#' using marker effects and a genetic map. If more than two traits are specified,
#' the function will also return predictions of the genetic correlation in the population
#' and the correlated response to selection.
#'
#' @param G.in See \code{G.in} in \code{\link[PopVar]{pop.predict}}.
#' @param y.in See \code{y.in} in \code{\link[PopVar]{pop.predict}}.
#' @param map.in See \code{map.in} in \code{\link[PopVar]{pop.predict}}.
#' @param crossing.table A \code{data.frame} with 2 columns (for bi-parental crosses) or 4 columns (for four-way crosses), each of which contains the names of parents to use in a potential cross. Rows contain individual crosses. See Details.
#' @param parents See \code{parents} in \code{\link[PopVar]{pop.predict}}.
#' @param n.parents Integer number of parents per cross. May be 2 or 4. If \code{crossing.table} is passed,
#'  this argument is ignored.
#' @param tail.p See \code{tail.p} in \code{\link[PopVar]{pop.predict}}.
#' @param self.gen The number of selfing generations in the potential cross. Can be an integer or \code{Inf} for
#' recombinant inbreds. Note: \code{self.gen = 1} corresponds to an F2 population.
#' @param DH Indicator if doubled-haploids are to be induced after the number of selfing generations indicated by
#' \code{self.gen}. For example, if \code{self.gen = 0} and \code{DH = TRUE}, then doubled-haploids are assumed
#' to be induced using gametes from F1 plants.
#' @param models See \code{models} in \code{\link[PopVar]{pop.predict}}.
#' @param n.core Number of cores for parallelization. Parallelization is supported
#' only on a Linux or Mac OS operating system; if working on a Windows system, the function
#' is executed on a single core.
#' @param ... Additional arguments to pass depending on the choice of \code{model}.
#'
#' @details
#'
#' Predictions are based on the deterministic equations specified by Allier et al. (2019).
#' 
#' In the case of four-way crosses (i.e. 4 parents), the function assumes that the first two parents are mated, 
#' producing a \eqn{F_1} offspring; then, the next two parents are mated, producing another \eqn{F_1} offspring. 
#' The two \eqn{F_1} offspring are then mated and inbreeding or doubled haploid induction (if specified) proceeds 
#' from there. For example, say cross \emph{i} uses parents P1, P2, P3, and P4. P1 and P2 are first mated, 
#' producing O1; then, P3 and P4 are mated, producing O2; then, O1 and O2 are mated, producing a segregating family.
#'
#' The \code{mppop.predict} function takes similarly formatted arguments as the \code{\link[PopVar]{pop.predict}} function
#' in the \code{PopVar} package. For the sake of simplicity, we also include the \code{mppop_predict2} function, which
#' takes arguments in a format more consistent with other genomewide prediction packages/functions.
#' 
#' If you select a \code{model} other than "rrBLUP", you must specify the following additional arguments:
#' \itemize{
#'   \item{\code{nIter}: See \code{\link[PopVar]{pop.predict}}}.
#'   \item{\code{burnIn}: See \code{\link[PopVar]{pop.predict}}}.
#' }
#' 
#' @return
#' A \code{data.frame} containing predictions of \eqn{\mu}, \eqn{V_G}, and \eqn{\mu_{sp}} for 
#' each trait for each potential multi-parent cross. When multiple traits are provided, the correlated 
#' responses and correlation between all pairs of traits is also returned.
#' 
#'
#' @references
#'
#' Allier, A., L. Moreau, A. Charcosset, S. Teyssèdre, and C. Lehermeier, 2019 Usefulness Criterion and Post-selection Parental
#' Contributions in Multi-parental Crosses: Application to Polygenic Trait Introgression. G3 (Bethesda) 9: 1469–1479.
#' https://doi.org/https://doi.org/10.1534/g3.119.400129
#' 
#' 
#' @examples
#' \donttest{
#' 
#' # Load data
#' data("think_barley")
#' 
#' # Vector with 8 parents
#' parents <- sample(y.in_ex$Entry, 8)
#' 
#' # Create a crossing table with four parents per cross
#' cross_tab <- as.data.frame(t(combn(x = parents, m = 4)))
#' names(cross_tab) <- c("parent1", "parent2", "parent3", "parent4")
#' 
#' out <- mppop_predict2(M = G.in_ex_mat, y.in = y.in_ex, map.in = map.in_ex, 
#'                       crossing.table = cross_tab, models = "rrBLUP")
#' 
#' }
#' 
#'
#'
#'
#' @importFrom qtl mf.h
#' @importFrom rrBLUP mixed.solve
#' @importFrom parallel mclapply
#' @importFrom utils combn
#'
#' @export
#'
mppop.predict <- function(G.in, y.in, map.in, crossing.table, parents, n.parents = 4, tail.p = 0.1,
                          self.gen = 10, DH = FALSE, models = c("rrBLUP", "BayesA", "BayesB","BayesC", "BL", "BRR"), 
                          n.core = 1, ...) {

  ###################
  # Error handling

  ## Check classes
  stopifnot(is.data.frame(G.in))
  stopifnot(is.data.frame(y.in))
  stopifnot(is.data.frame(map.in))

  # Check the number of markers in the map and the G.in objects
  if (ncol(G.in) - 1 != nrow(map.in))
    stop("The number of markers in G.in does not match the markers in map.in")

  ## Make sure there is no missing marker data
  if (any(is.na(G.in))) stop ("There must be no missing marker data.")

  # Self gen cannot be negative
  stopifnot(self.gen >= 0)
  # DH must be logical
  stopifnot(is.logical(DH))

  ## Make sure markers in G.in are in map.in
  markers_Gin <- as.character(unlist(G.in[1,-1, drop = TRUE]))
  markers_mapin <- map.in[,1, drop = TRUE]

  if (any(! markers_Gin %in% markers_mapin)) stop("The marker names in G.in are not all in map.in.")

  # Match arguments
  models <- match.arg(models)


  ###################


  ## Filter map.in to remove markers with unknown position
  map.in[[1]] <- as.character(map.in[[1]])
  map.in[[3]] <- as.numeric(map.in[[3]])
  map.in_use <- subset(map.in, !is.na(map.in[[3]]))
  # Reorder based on chrom and then pos
  map.in_use <- map.in_use[order(map.in_use[[2]], map.in_use[[3]]), , drop = FALSE]

  # Get the names of the markers in the new map
  markers_mapin <- as.character(map.in_use[[1]])

  ## Subset G.in for markers in map.in_use
  G.in_use <- G.in[, c(1, which(markers_Gin %in% markers_mapin) + 1), drop = FALSE]



  # If the crossing table is not missing, check that the parents are in the G.in input
  if (!missing(crossing.table)) {
    parents <- unique(unlist(crossing.table))

    # Make sure there are two or four parents
    npar <- ncol(crossing.table)

    ## If two parents, replicate the parents for the first "cross"
    if (npar == 2) {
      crossing.table <- crossing.table[,rep(seq_len(npar), each = 2), drop = FALSE]

    } else if (npar == 4) {
      # Good!

    } else {
      stop("Specified crosses must be 2-way or 4-way.")

    }

    # Make sure the parent names are not factors
    crossing.table <- as.data.frame(sapply(X = crossing.table, as.character), stringsAsFactors = FALSE)

  } else {
    if (missing(parents))
      stop("If no crossing.table is provided, a list of parents must be supplied.")

    parents <- sort(parents)
    # Create a crossing table with all possible parent combinations
    crossing.table <- as.data.frame(t(combn(x = parents, m = n.parents)), stringsAsFactors = FALSE)
    names(crossing.table) <- paste0("parent", seq_len(n.parents))

  }

  if (any(!parents %in% G.in_use[,1,drop = T]))
    stop("Parents are not in G.in.")

  ## If self.gen is Inf and DH is T, error
  if (is.infinite(self.gen) & DH) stop("Infinite selfing generations and doubled-haploid production cannot both occur.")

  ## Set the factors of line names in y.in to those in the marker df and those in the y.in
  lines_G.in <- sort(as.character(G.in_use[-1,1, drop = TRUE]))
  # The levels of the y.in_use geno factor should be the entry names in G.in
  y.in[[1]] <- factor(x = y.in[[1]], levels = lines_G.in)
  # Subset out NAs
  y.in_use <- subset(y.in, !is.na(y.in[[1]]))
  # Reorder lines
  y.in_use <- y.in_use[order(y.in_use[[1]]),, drop = FALSE]


  ## Number of traits and trait names
  n_traits <- ncol(y.in_use) - 1
  trait_names <- colnames(y.in_use)[-1]


  ## Convert the marker matrix into a useable matrix form
  G.in_pred <- sapply(X = G.in_use[-1, -1], function(x) as.numeric(as.character(x)))
  row.names(G.in_pred) <- as.character(G.in_use[-1,1])
  colnames(G.in_pred) <- as.character(unlist(G.in_use[1,-1]))
  # Reorder markers and lines
  G.in_pred <- G.in_pred[order(row.names(G.in_pred)), markers_mapin]


  ## Create a model.frame from the phenotypes; extract the column name of genotypes/lines
  geno_colname <- colnames(y.in_use)[1]

  # Subset using factors in y.in
  marker_names <- markers_mapin
  M <- G.in_pred[lines_G.in, marker_names, drop = FALSE]

  ## Pass arguments to mppop_predict2
  cross_predictions2 <- mppop_predict2(M = M, y.in = y.in_use, map.in = map.in_use, crossing.table = crossing.table,
                                       self.gen = self.gen, tail.p = tail.p, DH = DH, models = models, n.core = n.core, 
                                       ... = ...)

  ## Return the predictions
  return(cross_predictions2)

} # Close function







#' Predict genetic variance and genetic correlations in multi-parent populations
#' using a deterministic model
#'
#' @rdname mppop.predict
#'
#' @param M A Matrix of marker genotypes of dimensions \code{nLine} x \code{nMarker}, coded as
#' -1, 0, and 1.
#' @param marker.effects A data frame of marker effects. The first column should include the marker name and
#' subsequent columns should include the marker effects. Supercedes \code{y.in} if passed.
#'
#'
#' @importFrom qtl mf.h
#' @importFrom rrBLUP mixed.solve
#' @importFrom parallel mclapply
#' @importFrom methods formalArgs
#' @importFrom stats dnorm qnorm as.dist dist
#'
#' @export
#'
mppop_predict2 <- function(M, y.in, marker.effects, map.in, crossing.table, parents, n.parents = 4, tail.p = 0.1,
                           self.gen = 10, DH = FALSE, models = c("rrBLUP", "BayesA", "BayesB","BayesC", "BL", "BRR"), 
                           n.core = 1, ...) {

  ###################
  # Error handling

  ## Check classes
  stopifnot(is.matrix(M))
  stopifnot(is.data.frame(map.in))

  # Check the number of markers in the map and the G.in objects
  if (ncol(M) != nrow(map.in)) stop("The number of markers in G.in does not match the markers in map.in")

  ## Make sure there is no missing marker data
  if (any(is.na(M))) stop ("There must be no missing marker data.")

  # Self gen cannot be negative
  stopifnot(self.gen >= 0)
  # DH must be logical
  stopifnot(is.logical(DH))

  ## Make sure markers in G.in are in map.in
  stopifnot(!is.null(colnames(M)))
  stopifnot(!is.null(row.names(M)))


  markers_M <- colnames(M)
  markers_mapin <- map.in[[1]]


  # If the crossing table is not missing, check that the parents are in the G.in input
  if (!missing(crossing.table)) {
    parents <- unique(unlist(crossing.table))

    # Make sure there are two or four parents
    npar <- ncol(crossing.table)

    ## If two parents, replicate the parents for the first "cross"
    if (npar == 2) {
      crossing.table <- crossing.table[,rep(seq_len(npar), each = 2), drop = FALSE]

    } else if (npar == 4) {
      # Good!

    } else {
      stop("Specified crosses must be 2-way or 4-way.")

    }

    # Make sure the parent names are not factors
    crossing.table <- as.data.frame(sapply(X = crossing.table, as.character), stringsAsFactors = FALSE)

  } else {
    if (missing(parents))
      stop("If no crossing.table is provided, a list of parents must be supplied.")

    parents <- sort(parents)

    # Create a crossing table with all possible parent combinations
    crossing.table <- as.data.frame(t(combn(x = parents, m = n.parents)), stringsAsFactors = FALSE)
    names(crossing.table) <- paste0("parent", seq_len(n.parents))

  }



  if (any(! markers_M %in% markers_mapin)) stop("The marker names in M are not all in map.in.")

  # Get the names of genotyped entries
  geno_lines <- sort(row.names(M))


  ## Make sure one of y.in or marker.effects are not missing
  if (missing(y.in) & missing(marker.effects))
    stop("You must pass one of 'y.in' or 'marker.effects.'")

  ## Error check depending on what is passed
  if (!missing(marker.effects)) {

    # Check markers
    markers_maref <- marker.effects[[1]]

    if (! all(markers_M %in% markers_maref) ) stop("The marker names in M are not all in marker.effects")

    ## Number of traits and trait names
    n_traits <- ncol(marker.effects) - 1
    trait_names <- colnames(marker.effects)[-1]

    # Set boolean for later
    calc_marker_eff <- FALSE

    # Else check y.in
  } else {

    stopifnot(is.data.frame(y.in))

    # All phenotyped lines should be genotyped
    if (! all(y.in[[1]] %in% geno_lines) ) stop("All entries in 'y.in' should have marker data in 'M'.")

    ## Set the factors of line names in y.in to those in the marker df and those in the y.in
    y.in[[1]] <- factor(x = y.in[[1]], levels = geno_lines)
    # Subset out NAs
    y.in_use <- subset(y.in, !is.na(y.in[[1]]))
    # Reorder lines
    y.in_use <- y.in_use[order(y.in_use[[1]]),, drop = FALSE]


    ## Number of traits and trait names
    n_traits <- ncol(y.in_use) - 1
    trait_names <- colnames(y.in_use)[-1]

    # Set boolean for later
    calc_marker_eff <- TRUE

  }


  # Match arguments
  models <- match.arg(models)
  
  # Capture other arguments
  other.args <- list(...)


  ###################


  # Reorder map based on chrom and then pos
  map.in_use <- map.in[order(map.in[[2]], map.in[[3]]), , drop = FALSE]

  # Get the names of the markers in the new map
  markers_mapin <- as.character(map.in_use[[1]])


  ## Get list of unique parents from the crossing.table
  parents <- unique(unlist(crossing.table))
  # Make sure there are two or four parents
  npar <- ncol(crossing.table)

  ## If two parents, replicate the parents for the first "cross"
  if (npar == 2) {
    crossing.table <- crossing.table[,rep(seq_len(npar), each = 2), drop = FALSE]

  } else if (npar != 4) {
    stop("Specified crosses must be 2-way or 4-way.")


  }

  # Make sure the parent names are not factors
  crossing.table <- as.data.frame(sapply(X = crossing.table, as.character), stringsAsFactors = FALSE)

  if (!all(parents %in% row.names(M)))
    stop("Parents are not in G.in.")

  ## If self.gen is Inf and DH is T, error
  if (is.infinite(self.gen) & DH) stop("Infinite selfing generations and doubled-haploid production cannot both occur.")

  # Reorder markers
  M1 <- M[geno_lines, markers_mapin, drop = FALSE]

  ## Fit models to calculate marker effects, if necessary
  if (calc_marker_eff) {

    ## Create a model.frame from the phenotypes; extract the column name of genotypes/lines
    geno_colname <- colnames(y.in_use)[1]
    
    # Calculate marker effects
    marker_effect_out <- do.call("calc_marker_effects", 
                                 c(M = quote(M1), y.df = quote(y.in_use[-1]), models = models, 
                                   other.args[names(other.args) %in% formalArgs("calc_marker_effects")]))

    ## Create a complete matrix of marker effects for each trait
    mar_eff_mat <- do.call("cbind", lapply(marker_effect_out, "[[", "effects"))
    mar_eff_mat <- structure(mar_eff_mat, dimnames = list(row.names(mar_eff_mat), names(marker_effect_out)))
    
    mar_beta_mat <- do.call("cbind", lapply(marker_effect_out, "[[", "grand_mean"))
    mar_beta_mat <- structure(mar_beta_mat, dimnames = list(row.names(mar_beta_mat), names(marker_effect_out)))


    # Else create a matrix of marker ordered marker effects
  } else {

    # Create matrix
    mar_eff_mat <- as.matrix(marker.effects[,-1,drop = FALSE])
    row.names(mar_eff_mat) <- marker.effects[[1]]
    mar_eff_mat <- mar_eff_mat[markers_mapin,,drop = FALSE]

    # Set the grand mean to zero
    mar_beta_mat <- matrix(0, nrow = 1, ncol = ncol(mar_eff_mat),
                           dimnames = list(NULL, colnames(mar_eff_mat)))

  }



  ## Create an empty matrix
  marker_names <- markers_mapin

  # Split markers by chromosome
  map.in.chr <- split(map.in_use, map.in_use[,2, drop = FALSE])
  markers_chr <- lapply(map.in.chr, "[[", 1)

  ## Split marker effects by trait then chromosome
  mar_eff_mat_trait_chr <- apply(X = mar_eff_mat, MARGIN = 2, FUN = function(snps) {
    lapply(X = markers_chr, FUN = function(snp_chr) as.matrix(snps[snp_chr]))
  })

  # Calculate separate centimorgan distance matrices per chromosome
  chr_cM <- lapply(X = map.in.chr, FUN = function(x) as.matrix(dist(x[,3,drop = FALSE])))
  # Convert to recombination distance
  chr_c <- lapply(X = chr_cM, FUN = mf.h)


  # Function to calculate recombination covariance between QTL at various generations
  # of selfing
  ck <- function(c, k) if (k == 1) c else ((2*c) / (1 + 2*c)) * ( 1 - (0.5 * (1 - 2*c))^k )

  # Function for calulating disequilibrium between marker effects in progeny
  calc_parD <- function(X, .x, .y) crossprod(X[.x,,drop = FALSE] - X[.y,,drop = FALSE]) / 16



  ## Calculate 2 values for all pairs of markers. Value 1 will be multiplied by
  ## the D_ij values of the initiate cross (P1 x P2)(P3 x P4); and value 2 will be multiplied by
  ## the D_ij values of the second cross (F1(1) x F1(2)).
  # If DH and self.gen = 0, DH's are formed from the F1
  if (DH & self.gen == 0) {
    chr_covar_phi1 <- lapply(X = chr_c, function(c) (1 - (2*c)))
    chr_covar_phi2 <- replicate(n = length(chr_c), 1, simplify = FALSE)

    # DHs formed after k generations of selfing
  } else if (DH & self.gen > 0) {
    chr_covar_phi1 <- lapply(X = chr_c, function(c) (1 - (2*ck(c, self.gen)) + ck(c, self.gen-1) ) * (1 - (2*c)) )
    chr_covar_phi2 <- lapply(X = chr_c, function(c) (1 - (2*ck(c, self.gen))) )

    # RILs after self.gen selfing generations
  } else if (!DH & is.finite(self.gen)) {
    chr_covar_phi1 <- lapply(X = chr_c, function(c) (1 - ck(c, self.gen)) * (1 - (2*c)) )
    chr_covar_phi2 <- lapply(X = chr_c, function(c) ( (1 - (2*ck(c, self.gen))) - ((0.5 * (1 - (2*c)))^self.gen) ) )

  } else if (!DH & is.infinite(self.gen)) {
    stop("RILs after infinite genetations of selfing are not yet supported")

  }





  ## Predicted genotypic value of all genotypes
  ## Add non-zero grand mean if y_in was specified
  pgvs <- (M1 %*% mar_eff_mat) + matrix(mar_beta_mat, ncol = ncol(mar_eff_mat), nrow = nrow(M), byrow = TRUE)

  # Determine the k_sp from tail.p
  k_sp <- dnorm(x = qnorm(p = 1 - tail.p)) / tail.p

  ## Verify that all parents are homozygous for all markers
  parents_M <- M[parents,,drop = FALSE]
  # Make sure the markers are coded correctly
  if (!all(parents_M %in% c(-1, 1)))
    stop("This method assumes fully inbred parents. Therefore, marker genotypes other than -1 and 1 are not supported. Please
    edit the marker data and try again.")

  ## Change the number of cores based on system
  os <- Sys.info()["sysname"]
  n.core <- ifelse(os == "Windows", 1, n.core)

  ## Create a list to store dfs
  # Split up the crossing.table
  crossing_table_split <- split(x = crossing.table, rep(seq_len(n.core), length.out = nrow(crossing.table)))


  # Parallelize over cross predictions
  cross_predictions_out <- mclapply(X = crossing_table_split, FUN = function(ctb) {

    cross_predictions <- vector("list", nrow(ctb))
    for (j in seq_along(cross_predictions)) {

      # Character vector of the four parents - it was always be four
      pars <- as.character(ctb[j, 1:4])

      ## Subset the genotype matrix using the parents
      M_par <- M1[pars,,drop = FALSE]
      # Find segregating markers
      seg_mar <- which(!colMeans(M_par) %in% c(1, -1))
      seg_mar_chr <- lapply(markers_chr, intersect, marker_names[seg_mar])
      seg_mar_chr_num <- mapply(markers_chr, seg_mar_chr, FUN = function(.x, .y) which(.x %in% .y))

      par_geno <- lapply(seg_mar_chr, FUN = function(m) M_par[pars, m, drop = FALSE])

      # Calculate the LD covariance

      ## This depends on the number of selfing generations
      ## Formulae are taken from Allier et al 2019

      ## First calculate the LD between parental alleles in the
      ## first and second crosses

      # Disequilibrium in first cross
      D12 <- lapply(par_geno, calc_parD, 1, 2)  # first two parents
      D34 <- lapply(par_geno, calc_parD, 3, 4)  # last two parents

      phi_1 <- mapply(D12, D34, FUN = `+`, SIMPLIFY = FALSE)

      # Disequilibrium in second cross
      D13 <- lapply(par_geno, calc_parD, 1, 3)
      D14 <- lapply(par_geno, calc_parD, 1, 4)
      D23 <- lapply(par_geno, calc_parD, 2, 3)
      D24 <- lapply(par_geno, calc_parD, 2, 4)

      phi_2 <- mapply(D13, D14, D23, D24, FUN = function(.x, .y, .z, .a) .x + .y + .z + .a, SIMPLIFY = FALSE)


      ## Subset phi1 and phi2
      ## Subset the covariance coefficients
      ## Multiply
      if (DH & self.gen == 0) {
        chr_sigma <- mapply(
          phi_1, chr_covar_phi1, phi_2, chr_covar_phi2, seg_mar_chr_num, FUN = function(.a, .b, .c, .d, .e) {
            ((.c * .d[.e, .e]) + (.a * .b[.e, .e])) * .b[.e, .e]
          }, SIMPLIFY = FALSE)

      } else {
        chr_sigma <- mapply(
          phi_1, chr_covar_phi1, phi_2, chr_covar_phi2, seg_mar_chr_num, FUN = function(.a, .b, .c, .d, .e) {
            ((.c * .d[.e, .e]) + (.a * .b[.e, .e]))
          }, SIMPLIFY = FALSE)
      }


      # Eq. 2 Allier et al 2019
      # Calculate variance per trait
      pred_varG_j <- sapply(X = mar_eff_mat_trait_chr, FUN = function(trait_mars) {
        sum(mapply(chr_sigma, trait_mars, seg_mar_chr, FUN = function(.x, .y, .z) crossprod(.y[.z,], .x) %*% .y[.z,]))
      })


      ## Predictions
      # Cross mean
      pred_mu_j <- colMeans(pgvs[pars,,drop = FALSE])


      # Genetic correlations between traits, if more than one trait
      if (n_traits > 1) {

        # Create a pairwise matrix
        pred_corG_j_mat <- matrix(data = 0, nrow = n_traits, ncol = n_traits,
                                  dimnames = list(trait_names, trait_names))
        # Convert to distance object
        pred_corG_j_dist <- as.dist(pred_corG_j_mat)

        ## Create a combination matrix
        trait_combn <- combn(x = seq(n_traits), m = 2)

        # Turn the marker effect list inside-out
        mar_eff_mat_chr_trait <- lapply(
          X = seq_along(seg_mar_chr), function(q) lapply(mar_eff_mat_trait_chr, "[[", q)
        )

        ## Apply over columns
        # Calculate genetic covariance between pairs of traits
        pred_corG_j <- sapply(X = seq_len(ncol(trait_combn)), FUN = function(rho) {
          ij <- trait_combn[,rho]
          # Calculate covariance
          mar_eff_mat_chr_trait_ij <- lapply(X = mar_eff_mat_chr_trait, FUN = "[", ij)
          covar <- sum(mapply(chr_sigma, mar_eff_mat_chr_trait_ij, seg_mar_chr,
                              FUN = function(.x, .y, .z) crossprod(.y[[1]][.z,], .x) %*% .y[[2]][.z,]))

          # Subset predicted variance to calculate correlation
          corG <- covar / prod(sqrt(pred_varG_j[ij]))
        })

        # Add these predictions to the distance object, then convert to matrix
        pred_corG_j_dist[seq_along(pred_corG_j_dist)] <- pred_corG_j

        ## Convert distance object to matrix
        pred_corG_j_mat <- as.matrix(pred_corG_j_dist)
        diag(pred_corG_j_mat) <- NA

        ## Calculate correlated progeny mean
        response_trait_varG <- matrix(pred_varG_j, nrow = length(pred_varG_j),
                                      ncol = length(pred_varG_j), byrow = TRUE)
        correlated_response <- k_sp * pred_corG_j_mat * sqrt(response_trait_varG)
        pred_mu_j_mat <- matrix(pred_mu_j, nrow = length(pred_mu_j), ncol = length(pred_mu_j), byrow = TRUE)
        pred_cor_musp_low <- pred_mu_j_mat - correlated_response
        pred_cor_musp_high <- pred_mu_j_mat + correlated_response

        # Change names
        colnames(pred_cor_musp_low) <- paste0("pred_cor_musp_low_", trait_names)
        colnames(pred_cor_musp_high) <- paste0("pred_cor_musp_high_", trait_names)

      } else {
        pred_corG_mat <- pred_cor_musp_low <- pred_cor_musp_high <- pred_corG_j_mat <- NULL

      }

      cross_predictions[[j]] <- data.frame(
        t(pars),
        trait = trait_names,
        cbind(pred_mu = pred_mu_j, pred_varG = pred_varG_j, pred_corG_j_mat,
              pred_cor_musp_low, pred_cor_musp_high),
        stringsAsFactors = FALSE, row.names = NULL)

    }

    do.call("rbind", cross_predictions)

  }, mc.cores = n.core) # End loop

  ## Bind rows
  cross_predictions1 <- do.call("rbind", cross_predictions_out)
  # Rename col 1:npar
  names(cross_predictions1)[1:4] <- paste0("parent", 1:4)

  ## Calculate response predictions (superior progeny, correlated response, etc.)
  pred_response <- (k_sp * sqrt(cross_predictions1$pred_varG))

  # Superior progeny mean
  cross_predictions1[["pred_musp_low"]] <- cross_predictions1$pred_mu - pred_response
  cross_predictions1[["pred_musp_high"]] <- cross_predictions1$pred_mu + pred_response

  ## Re-order columns
  cor_W_cols <- grep(pattern = paste0(trait_names, collapse = "|"), x = names(cross_predictions1))
  cross_predictions2 <- cross_predictions1[, c(setdiff(seq(ncol(cross_predictions1)), cor_W_cols), cor_W_cols)]
  row.names(cross_predictions2) <- NULL

  ## Return the predictions
  return(cross_predictions2)

} # Close function



