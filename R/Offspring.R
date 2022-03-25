#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Initiate the dataframe with candidate parents
#'
#' @description For all individuals sharing a pedigree parent, all candidate
#'   field parents of at least 1 member are listed, as well as 'INVIS'.
#'
#' @details Id's sharing the same par.ped are considered a sibship cluster, and
#'   all par.field occuring with at least one id in the cluster are considered a
#'   candidate parent for the whole cluster. Resulting Id - par.field pairs not
#'   occuring in ParField are given 'NA' for expl.var, except when
#'   expl.var="Ndays", when they get '0'.
#'
#'   If \code{expl.var='carer'}, there will be a warning about pedigree sibling
#'   clusters with >1 carer.
#'
#'   'INVIS' is added as candidate parent to all sibship clusters. If
#'   \code{expl.var='Ndays'}, the expl.var value for 'INVIS' for each id is
#'   calculated as 11 - the sum of Ndays with known par.field for that id,
#'   otherwise expl.var = NA for INVIS.
#'
#' @param  ParPedigree parent in genetically inferred pedigree, dataframe with
#'   columns id-dam or id-sire, column names ignored and renamed to 'id' and
#'   'par.ped'. May include dummy individuals and mystery samples.
#' @param  ParField candidate parents inferred in the field, dataframe with
#'   columns id - par.field - [expl.var].
#' @param  expl.var  explanatory variable to calculate probability that
#'   par.field is identical to par.ped. Must be a column name in `ParField'.
#'
#'    Defaults to 'prob', and id-par.field
#'   combinations not listed in 'ParField', as well as 'INVIS', will become
#'   'NA'.
#'
#' @return a dataframe with columns id - par.ped - par.field - [expl.var]
#'
#' @seealso \code{\link{parDFtoM}, \link{ParentOfAll}}
#'
#' @examples
#' \dontrun{
#' FieldMums <- setNames(tblLife[!is.na(tblLife$MumCode), c("Code", "MumCode")],
#'                      c("id", "par.field"))
#' FieldMums$carer <- 1
#' cand.dams <- Init.candpar(unique(Pedigree[,c("id", "dam")]), FieldMums,
#'                          expl.var="carer")
#'
#' cand.sires <- Init.candpar(Pedigree[,c("id", "sire")],
#'                           setNames(DaysHeld[, c("Code","StagCode","Ndays")],
#'                                    c("id", "par.field", "Ndays")),
#'                           expl.var="Ndays")
#' }
#'
#' @export

Init.candpar <- function(ParPedigree, ParField, expl.var="prob") {
	candpar.tmp <- merge(setNames(ParPedigree, c("id", "par.ped")),
	                     ParField[, c("id", "par.field", expl.var)], all=TRUE)
	for (x in c("par.ped", "par.field"))  candpar.tmp[,x] <- as.character(candpar.tmp[,x])
	candpar.tmp <- candpar.tmp[!(is.na(candpar.tmp$par.ped) & is.na(candpar.tmp$par.field)), ]
	candpar.tmp$par.field[is.na(candpar.tmp$par.field)] <- "INVIS"

	candpar.L <- dlply(candpar.tmp[!is.na(candpar.tmp$par.ped), ], "par.ped", parDFtoM, expl.var)
	if (expl.var == "Ndays") {  # ! Rum deer - specific
	  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	  CalcInvisDays <- function(M) {
	    inv <- which(colnames(M) == "INVIS")
	    M[,-inv][is.na(M[,-inv])] <- 0
	    M[,inv] <- pmax(0, 11 - rowSums(M[,-inv, drop=FALSE]))
	    return( M )
	  }
	  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	  candpar.L <- llply(candpar.L, function(M) CalcInvisDays(M) )

  } else if (expl.var == "carer") {
	  ncarers <- sapply(candpar.L, function(M) sum(colnames(M)!="INVIS"))
	  if (any(ncarers>1)) {
	    assign("candpar.warn", TRUE, envir = sys.frame(sys.nframe()-1))  # assign value 1 level up
	    warning("The pedigree includes maternal siblings with different field mums!",
	            immediate. = TRUE)
	    for (i in which(ncarers>1)) {
	      cat(names(candpar.L)[i], " : ", colnames(candpar.L[[i]]), "\n")
	    }
	  }
	  hascarer <- function(M)  apply(M, 1, function(x) any(!is.na(x)))
	  candpar.L <- llply(candpar.L, function(M) {M[is.na(M) & hascarer(M)] <- 0; M})
	}
	candpar.df <- ldply(candpar.L, function(M) adply(M, 1:2, .expand=FALSE),
                     .id="par.ped")
  names(candpar.df)[names(candpar.df)=="V1"] <- expl.var
  for (x in c("id", "par.field", "par.ped")) {
   candpar.df[,x] <- as.character(candpar.df[,x])
  }
	return( candpar.df[, c("id", "par.field", "par.ped", expl.var)] )  # fix column order
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title turn long-format dataframe into a matrix
#'
#' @description to be used in llply etc.
#'
#' @param df dataframe subset, with at least columns 'id' and 'par.field'
#' @param expl.var name of column to become matrix contents
#'
#' @return a matrix with 'id' in rows and 'par.field' in columns. If "INVIS" is
#'   not among 'par.field', it will be added.
#'
#' @examples
#' \dontrun{
#' candpar.L <- dlply(candpar.df, "par.ped", parDFtoM, expl.var="prob")
#' }
#' @export

parDFtoM <- function(df, expl.var="prob") {
	M <- adrop(daply(df, c("id", "par.field"), function(x) x[,expl.var], .drop_o=FALSE), drop = 3)
	if (!"INVIS" %in% colnames(M)) {
    M <- abind(M, "INVIS" = matrix(NA,nrow(M),1), along=2)  ## cbind doesn't work if only 1 row
  }
	names(dimnames(M)) <- c("id", "par.field")
	return( M )
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Column products
#'
#' @description probabilities for each candidate parent (columns) to be the
#'   parent of all offspring (rows), given a matrix of per-offspring per-parent
#'   probabilities for a single par.ped.
#'
#' @details  The result for 'INVIS' (invisible candidate parent) includes a
#'   correction factor, as the invisible parent of offspring A may differ from
#'   the invisible parent of B, as a function of the estimated total number of
#'   invisible candidate parents (\code{nInvis}, \eqn{n_I}) and the number of
#'   siblings \eqn{n_O}:
#'   \deqn{P(\text{INVIS}) = P'(\text{INVIS}) \times
#'          \left(\frac{1}{n_I}\right)^{n_O} \times n_I}
#'
#' @param M matrix, with 'id' in rows and 'par.field' in columns, for a single
#'   par.ped
#' @param nInvis single number, estimated (average) number of potential
#'   invisible candidate parents for each individual
#'
#' @return a numeric vector with log10-probabilities, of length equal to the
#'   number of columns in M
#'
#' @examples
#' \dontrun{
#'   candpar.L <- dlply(candpar.df, "par.ped", parDFtoM, expl.var="prob")
#'	 Candpar.L2 <- llply(candpar.L, ParentOfAll, nInvis[p])
#' }
#' @export

ParentOfAll <- function(M, nInvis=10)
{
  PrPar <- c(apply(M, 2, prod))
  if ("INVIS" %in% colnames(M)) {
    PrPar["INVIS"] <- PrPar["INVIS"] * nInvis*(1/nInvis)^nrow(M)  # prob. same invis parent of all
  }
  return( log10(PrPar) )
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Find best-matching field parent for each pedigree sibship
#'
#' @description Match each sibship cluster, sharing the same (dummy or
#'   genotyped) pedigree parent, to a candidate field parent or INVIS, ensuring
#'   that each field parent is matched to at most one pedigree sibship.
#'
#' @param CPM matrix with with log10-probabilities that field parent j (columns)
#'   is the parent of all individuals in sibship i (rows) (Log-Parent-of-All).
#' @param CandPar.OK  dataframe with id's matched in the previous iteration of
#'   \code{\link{MatchIDsPedigreeField}}, with at least columns id - Cand.id -
#'   L.off.
#' @param NonUnique vector with column names that do not represent a unique
#'   field parent and thus may be matched to multiple pedigree sibhips, e.g.
#'   'INVIS'.
#' @param LL.th numeric vector with assignment thresholds for each iteration, on
#'   log10 scale.
#' @param verbose print the update step and new log-likelihood to the console
#'   each iteration.
#'
#' @return a matrix (\code{M.allow}) with the same dimensions as CPM. In case of
#'   a conclusive match between pedigree parent 'i' and field parent
#'   'j', M.allow[i,j]=1 (allowed), M.allow[i,-j]=NA (disallowed) and
#'   M.allow[-i,j]=NA. In absence of a conclusive match, M.allow[i,]=1, except
#'   for those j's where already a conclusive match is made, and vice versa for
#'   M.allow[,j]. Coding as 'NA' instead of '0' enables operations like
#'   \code{which.max(CPM[i,] * M.allow[i,])} and the like.
#'
#' @examples
#' \dontrun{
#' if (r>1) {
#'   CandPar.OK <- cand.ok[[r-1]]  # confirmed and/or name replaced at end prev. iteration
#' } else {
#'   CandPar.OK <- NULL
#' }
#' M.allow <- FilterCand(CPM = candpar.M,
#'                       CandPar.OK = CandPar.OK,
#'                       NonUnique = c("INVIS", "REF"),   #, "SMAL", "YOUN", "MATU"),
#'                       LL.th = c(seq(10,2,-1), seq(1,.5, -.1)),
#'                       verbose = TRUE)
#' candpar.M <- candpar.M * M.allow
#' }
#' @export

FilterCand <- function(CPM = NULL,
                       CandPar.OK = NULL,
                       NonUnique = c("INVIS", "NONE", "REF"), # "SMAL", "YOUN", "MATU"),
                       LL.th = c(seq(10,2,-1), seq(1,.5, -.1)),  # on log10 scale
                       verbose = TRUE)
{
	Pfield <- colnames(CPM)
	NonUnique <- unique(c(NonUnique, "INVIS"))

	M.allow <- CPM  # same dimensions & names
	M.allow[,] <- 1
  if (is.null(CandPar.OK)) {
    GOM.init <- M.allow
    GOM.init[,] <- 0
    GOM.init[,"INVIS"] <- 1

  } else {
		if (!is.data.frame(CandPar.OK)) stop("CandPar.OK must be dataframe or NULL")
		if (!all(c("id", "Cand.id") %in% names(CandPar.OK)))  stop("CandPar.OK must have columns 'id' and 'Cand.id'")

    Par.OK <- CandPar.OK[CandPar.OK$Cand.id %in% rownames(CPM), c("id", "Cand.id", "L.off")]     # confirmed + renamed
    Renamed <- !Par.OK$id %in% rownames(CPM)  # renamed dummy + mystery samples
    Unseen <- !Par.OK$Cand.id %in% colnames(CPM)  # not candidate parent for any indiv
    Sneaky <- with(Par.OK, Cand.id[is.na(L.off) & Cand.id %in% colnames(CPM)])  # not among candidate parents *for its offspring*
    Par.OK$id[Renamed] <- Par.OK$Cand.id[Renamed]
    Par.OK$Cand.id[Unseen] <- "INVIS"
    M.allow[Par.OK$id, ] <- NA
    M.allow[, Par.OK$Cand.id] <- NA
    M.allow[as.matrix(Par.OK[,c("id", "Cand.id")])] <- 1

		CPM[Sneaky, Sneaky] <- CPM[Sneaky, "INVIS"]  # use 'INVIS' prob to calculate total LL

    GOM.init <- M.allow
    GOM.init[is.na(GOM.init)] <- 0
    matched <- rowSums(M.allow[, !Pfield %in% NonUnique], na.rm=TRUE) == 1
    GOM.init[!matched, ] <- 0
    GOM.init[!matched, "INVIS"] <- 1
    if (!all(rowSums(GOM.init) == 1))  stop("Something wrong with GOM.init")
  }

	LogLik <- numeric(length(LL.th))
	GOM <- GOM.init

	#@@@@@@@@@@@@@@@@@@
	dm2 <- function(x) Max(x, m=NA) - Max2(x, m=-Inf)
	#@@@@@@@@@@@@@@@@@@

	for (s in seq_along(LL.th)) {
		d <- apply(CPM * M.allow, 1, dm2)
		II <- order(d, decreasing=TRUE)

		for (i in II) {
			d.i <- dm2(CPM[i,] * M.allow[i,])
			if (!is.na(d.i) & (d.i >= LL.th[s])) {
				j <- which.max(CPM[i,] * M.allow[i,])
				if (!Pfield[j] %in% NonUnique) {
					a.j <- rep(NA, nrow(CPM))
					for (x in which(!is.na(CPM[,j]))) {
						if (!is.na(M.allow[x,j])) {
							a.j[x] <- CPM[x,j] - Max(CPM[x,-j]*M.allow[x,-j], m=-Inf)
						}
					}
  				if (a.j[i]==Inf || (a.j[i] - Max(a.j[-i], m=-Inf)) >= 3*LL.th[s]) {
  				  if (is.na(CPM[i,j]))  stop("Problem in FilterCand!")
  					GOM[i, j] <- 1
  					GOM[i, -j] <- 0
  					M.allow[i, -j] <- NA
  					M.allow[i, Pfield %in% NonUnique] <- 1
  					if (!Pfield[j] %in% NonUnique) {
  						M.allow[-i, j] <- NA
  					}
  				}
				}
			}
		}
		LogLik[s] <- sum(CPM[GOM==1])
		if(verbose)  cat(s, "\t", LL.th[s], "\t: ", sum(GOM[,setdiff(Pfield,"INVIS")]==1),
			"\t", round(LogLik[s],3), "\n")
	}
	if(verbose)  cat("\n")
 	return( M.allow )

}
