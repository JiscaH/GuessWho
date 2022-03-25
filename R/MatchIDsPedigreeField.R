#============================================================================
#============================================================================
#' @title Match IDs between a genetic pedigree and field data
#'
#' @description Match dummy and mystery IDs in the pedigree to field IDs, and
#'   check for consistency between pedigree and field-identified (candidate)
#'   parents.
#'
#' @details All matches i, j with consensus level 'OK' are accepted, as well
#'   as 'LIKELY' in the second and third iteration, and 'CHECKMUM',
#'   'CHECKYEARS', and 'CHECKSEX' in subsequent iterations. Warnings are issued
#'   about 'MISMATCH', 'MERGE?', and 'CHECKMUM', as well as of cases where one
#'   field ID is accepted for two or more pedigree IDs (receiving consensus
#'   level 'MULTI'). For explanation of these levels, see
#'   \code{\link{init.PedFieldMatch}}
#'
#' @param  Pedigree genetically inferred pedigree, dataframe with columns
#'   id-dam-sire. May include dummy individuals and mystery samples.
#' @param  LHextra dataframe with lifehistory data, with columns id - Sex -
#'   BirthYear - FirstRepro - LastRepro. See \code{\link{CalcFirstLastRepro}}.
#'   May include field IDs (genotyped & not-genotyped) as well as mystery
#'   samples.
#' @param  FieldDams  dataframe with id - dam.field - [expl.var], where the 3rd
#'   column is a variable associated with the probability of parentage (e.g.
#'   'carer')
#' @param  FieldSires  dataframe with id - sire.field - [expl.var], as
#'   'FieldDams' (with expl.var e.g. 'Ndays')
#' @param  FieldIDs  character vector with IDs of real individuals (not mystery
#'   samples), both genotyped and non-genotyped. If \code{NULL}, defaults to all
#'   id's in LHextra, and all id's and field parents in FieldDams and
#'   FieldSires.
#' @param MumProbs  numeric vector length 2, with probability that field mum
#'   (carer) is the pedigree mum, and probability that an arbitrary female is
#'   the pedigree mum for an individual without field mum. [Note: having a GLM
#'   model as input as for sires is fairly straightforward to implement.]
#' @param DadMod  A fitted binomial GLM, with pedigree parent/not as response
#'   variable and [expl.var] in FieldSires and optionally field parent age as
#'   explanatory variable.
#' @param AgeFirstRepro vector of length 2, population minimum age of first
#'   reproduction for females and males, respectively.
#' @param AgeLastRepro vector of length 2, population maximum age of last
#'   reproduction for females and males, respectively.
#' @param YearLastRepro vector of length 2, maximum calendar years of last
#'   reproduction after calendar year of death.
#' @param MonthLastRepro vector of length 2, minimum calendar month of last
#'   reproduction in calendar year of death; see
#'   \code{\link{CalcFirstLastRepro}} for details.
#' @param nIter Number of iterations
#' @param nInvis  vector of length 2, number of 'invisible' candidate parents
#'   per pedigree sibship; see \code{\link{ParentOfAll}}.
#' @param Lref  list with 2 numeric vectors, 'dam' and 'sire', with reference
#'   log10-probabilities for candidate parents as a function of sibship size;
#'   the probability of a candidate field parent must at least pass this value,
#'   see \code{\link{MkOffCat}}.
#' @param T.offcat Numeric vector of length 'nIter'; threshold log10-probability
#'   between most-likely and next-most-likely candidate parent to declare a
#'   conclusive match, see \code{\link{MkOffCat}}
#' @param nag  logical, wait for user confirmation to continue when encountering
#'   mismatches (TRUE), or always continue with fingers crossed (FALSE)?
#'
#' @return A list with
#' \describe{
#'   \item{Pedigree}{The pedigree, with id's replaced by field id's where
#'   possible, including in dam and sire columns. Rownames show the original id,
#'   and columns 'type.id', 'type.dam', and 'type.sire' indicate whether the
#'   name in the respective column is a proper, \strong{G}enotyped pedigree +
#'   field ID ('G'); or is a dummy/mystery pedigree ID \strong{R}eplaced by a
#'   field ID ('R'); or is a non-replaced \strong{D}ummy/mystery pedigree ID
#'   ('D')}
#'   \item{Likely.r}{a list of length nIter with dataframes with shortlisted
#'     candidate field ids for each pedigree id}
#'   \item{Match.r}{a list of length nIter with dataframes with accepted
#'     matches (a subset of the corresponding Likely.r element)} }
#'
#' When \code{nag=TRUE} and choosing not to continue when a mismatch is
#' encountered, the object \code{Match.dump} is created in the Global
#' environment, which is a list with a selected subset of intermediate results
#'
#' @examples
#' \dontrun{
#' OUT <- MatchIDsPedigreeField(Pedigree = Pedigree,
#'                              LHextra = LHextra,
#'                              FieldDams = FieldMums,
#'                              FieldSires = FieldCandSires,
#'                              FieldIDs = FieldIDs,
#'                              MumProbs = c("CarerIsMum" = 0.995,
#'                                           "NoCarer_IsMum" = 0.05),
#'                              DadMod = mod.DaysHeld,
#' #                             Repro = list(
#'                                  AgeFirstRepro = c(3, 2),
#'                                  AgeLastRepro = c(25,20),
#'                                  MonthLastRepro = c(9,1),
#'                                  YearLastRepro = c(0, 1), #)
#'                              nInvis = c(3, 7),
#'                              Lref = list(dam = Lref.dam,
#'                                          sire = Lref.sire),
#'                              nIter = 5,
#'                              nag=TRUE)
#' }
#'
#' @importFrom stats setNames na.exclude predict
#' @importFrom plyr dlply ldply llply daply adply
#' @importFrom abind abind adrop
#'
#' @export

MatchIDsPedigreeField <- function(Pedigree = NULL,
                                  LHextra = NULL,
                                  FieldDams = NULL,
                                  FieldSires = NULL,
                                  FieldIDs = NULL,
                                  MumProbs = c("CarerIsMum" = 0.995,
                                               "NoCarer_IsMum" = 0.05),
                                  DadMod = NULL,
                           #       Repro = list(
                                           AgeFirstRepro = c(3, 2),
                                           AgeLastRepro = c(25,20),
                                           YearLastRepro = c(0, 1),
                                           MonthLastRepro = c(9,1),  #),
                                  nIter = 5,
                                  nInvis = c(3, 7),
                                  Lref = NULL,
                                  T.offcat=c(5, 3,2, rep(1, nIter-2)),
                                  nag  = TRUE)
{

  #~~~  Check input & initiate  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (is.null(Pedigree))  stop("Please provide 'Pedigree'")
  if (is.null(LHextra))  stop("Please provide 'LHextra'")
  if (is.null(FieldDams))  stop("Please provide 'FieldDams'")
  if (is.null(FieldSires))  stop("Please provide 'FieldSires'")
  if (nIter <3)  stop("nIter must be at least 3")

  for (x in c("id", "dam", "sire"))  Pedigree[,x] <- as.character(Pedigree[,x])
  for (x in c("id"))  LHextra[,x] <- as.character(LHextra[,x])
  for (x in c("id", "par.field"))  FieldDams[,x] <- as.character(FieldDams[,x])
  for (x in c("id", "par.field"))  FieldSires[,x] <- as.character(FieldSires[,x])

  if (nrow(Pedigree) > length(unique(Pedigree$id)))  stop("Pedigree must not have duplicate ids")
  if (nrow(LHextra) > length(unique(LHextra$id)))  stop("LHextra must not have duplicate ids")

  if (is.null(FieldIDs)) {
    FieldIDs <- unique(na.exclude(c(LHextra$id,
                                    FieldDams$id, FieldDams$par.field,
                                    FieldSires$id, FieldSires$par.field)))
  } else {
    if (!all(c(FieldDams$id, FieldDams$par.field,FieldSires$id, FieldSires$par.field) %in%
        FieldIDs))  stop("All id and par.field in 'FieldDams' and 'FieldSires' must be in FieldIDs")
    FieldIDs <- unique(na.exclude(FieldIDs))
  }
  if (!all(FieldIDs %in% LHextra$id)) {
    LHextra <- rbind(LHextra,
                     merge(LHextra[0,],
                           data.frame(id=setdiff(FieldIDs, LHextra$id)),
                           all=TRUE))
  }

  #~~~~~~~~~
  OK.accept = c(list("1"= "OK"),
                rep(list(c("OK", "LIKELY")), 2),
                rep(list(c("OK", "LIKELY", "MULTI", "CHECKMUM", "CHECKYEARS", "CHECKSEX")), nIter-3))
  Match.r <- vector("list", length=nIter)
  Likely.r <- vector("list", length=nIter)

  rownames(Pedigree) <- Pedigree$id
  Pedigree.init <- Pedigree

  #~~~~~~~~~
  # Calc individual years of potential first & last reproduction
  if (!all(c("FirstRepro","LastRepro") %in% colnames(LHextra))) {
    LHextra <- CalcFirstLastRepro(LHextra,
                                  AgeFirstRepro = AgeFirstRepro,
                                  AgeLastRepro = AgeLastRepro,
                                  YearLastRepro = YearLastRepro,
                                  MonthLastRepro = MonthLastRepro)
  }


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (r in 1:nIter) {
    cat("\n\n", format(Sys.time(), "%a %d %b %Y %X"),
        "\t round ", r, " of ", nIter, "\n\n")

		#~~~  Label, Sex, Age, Mother  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    LH.chk <- CheckBirthYears(Pedigree,    # Pedigree updated at end previous iteration
                              LHextra,
                              AgeFirstRepro = AgeFirstRepro,
                              AgeLastRepro = AgeLastRepro)
    utils::flush.console()
    if (nag & any(!LH.chk$OK)) {
      ANS <- readline(prompt = "\n\n Continue Y/N? ")
      if (!substr(ANS, 1, 1) %in% c("Y", "y", "J", "j", "")) {
        assign("Match.dump", list(LH.chk = LH.chk,
                                  Likely.r = Likely.r,
                                  Match.r = Match.r),
               envir = .GlobalEnv)
				stop("Dumping data in 'Match.dump', and stopping")
			}
    }
    LH.pedigree <- merge(Pedigree[, c("id", "dam")],
                         LH.chk[, c("id", "Sex", "BY.min", "BY.max", "FirstRepro.off", "LastRepro.off")],
                         all.x=TRUE)
    LH.pedigree <- LH.pedigree[match(Pedigree$id, LH.pedigree$id), ]  # fix order
    rownames(LH.pedigree) <- rownames(Pedigree)  # original IDs
#    LH.pedigree <- unique(LH.pedigree)  # in case of 'MULTI' match due to erroneously non-merged sibships
    LH.field <- merge(LHextra[,c("id", "Sex", "BirthYear", "FirstRepro", "LastRepro")],
                      FieldDams, all=TRUE)
    LH.field <- LH.field[match(FieldIDs, LH.field$id), ]

    cat("\n Comparing label, sex, birth years, offspring of mothers ... \n")
    Cat.ID <- MkIDcat(LH.pedigree, LH.field, ids.snpd.field = intersect(rownames(Pedigree), FieldIDs))


    #~~~  Offspring  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # based on field candidate-parents of pedigree offspring

		LL.off <- matrix(NA, dim(Cat.ID)[2], dim(Cat.ID)[3],
                     dimnames = dimnames(Cat.ID)[2:3])

		for (p in 1:2) {
			if (p==1) {
			  cat("\n Comparing field dams of pedigree maternal sibs ... \n")
				candpar.warn <- FALSE
				cand.dams <- Init.candpar(unique(Pedigree[,c("id", "dam")]), FieldDams,
				                          expl.var=colnames(FieldDams)[3])
				if (nag & candpar.warn) {
				  ANS <- readline(prompt = "\n\n Continue Y/N? ")
				  if (!substr(ANS, 1, 1) %in% c("Y", "y", "J", "j", "")) {
				    assign("Match.dump", list(LH.field = LH.field,
				                              LH.pedigree = LH.pedigree,
				                              Cat.ID = Cat.ID,
				                              cand.dams = cand.dams,
				                              Likely.r = Likely.r,
				                              Match.r = Match.r),
				           envir = .GlobalEnv)
				    stop("Dumping data in 'Match.dump', and stopping")
				  }
				}
				cand.dams$prob <- with(cand.dams, ifelse(is.na(carer),
																					ifelse(id %in% FieldIDs & par.ped %in% FieldIDs,
																						 MumProbs["NoCarer_IsMum"],
																						 0.5),   # Dummies & Mystery
																					ifelse(carer == 1,
																							MumProbs["CarerIsMum"],
																							1 - MumProbs["CarerIsMum"])))
				candpar.df <- cand.dams

			} else if (p==2) {
			  cat("\n Comparing field sires of pedigree paternal sibs ... \n")
			  cand.sires <- Init.candpar(unique(Pedigree[,c("id", "sire")]), FieldSires,
			                             expl.var=colnames(FieldSires)[3])
				cand.sires$Age.f <- CalcParAge(cand.sires, LHextra, Factor=TRUE)
				cand.sires$prob <- predict(DadMod, newdata=cand.sires, type="response")
				cand.sires$prob[is.na(cand.sires$prob)] <- 0   # dead / not born yet
				candpar.df <- cand.sires[, c("id", "par.ped", "par.field", "prob")]
			}

			candpar.L <- dlply(candpar.df, "par.ped", parDFtoM, expl.var="prob")
			candpar.L2 <- llply(candpar.L, ParentOfAll, nInvis[p])   # ParentOfAll(): log10(apply(M, 2, prod))

			parfield.ids <- unique(candpar.df$par.field)
			candpar.M <- matrix(NA, length(candpar.L2), length(parfield.ids),
													dimnames = list(names(candpar.L2), parfield.ids))
			for (i in 1:length(candpar.L2)) {
				candpar.M[i, names(candpar.L2[[i]])] <- candpar.L2[[i]]
			}

			if (p == 2) {  # filter out unlikely/impossible candidate sires for each pedigree sibship
			  if (r>1) {
			    CandPar.OK <- Match.r[[r-1]]  # confirmed and/or name replaced at end prev. iteration
			  } else {
			    CandPar.OK <- NULL
			  }
			  M.allow <- FilterCand(CPM = candpar.M,
			                        CandPar.OK = CandPar.OK,
			                        NonUnique = c("INVIS", "REF"),   #, "SMAL", "YOUN", "MATU"),
			                        LL.th = c(seq(10,2,-1), seq(1,.5, -.1)),  # on log10 scale
			                        verbose = TRUE)
				candpar.M <- candpar.M * M.allow
			}

			LL.off[rownames(candpar.M), colnames(candpar.M)] <- candpar.M
		}

		#~~~~~~~~~
		# reference probabilities (minimum)

		Ref.i <- setNames(numeric(nrow(Pedigree)), Pedigree$id)
		nOff.dam <- c(table(Pedigree$dam))
		Ref.i[names(nOff.dam)] <- Lref[["dam"]][nOff.dam]
		nOff.sire <- c(table(Pedigree$sire))
		Ref.i[names(nOff.sire)] <- Lref[["sire"]][nOff.sire]


    #~~~  Match  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # rename highly-likely matches of dummies & mystery ids to build upon in next iteration

		Likely.L <- setNames(vector("list", length=dim(Cat.ID)[2]), rownames(Pedigree))
		# Pedigree$id not necessarily unique after several iterations, if erroneous non-merged sibships
		for (i in 1:dim(Cat.ID)[2]) {
		  these <- which(Cat.ID["Label",i,] == "Match" | Cat.ID["Mother",i, ] == "Match" | !is.na(LL.off[i,]))
		  if (length(these)==0)  next
	    Cat.off <- MkOffCat(LL.off[i,][these, drop=FALSE],  REF=Ref.i[i],
	                        T.off = T.offcat[r])
	    Cat.i <- t(Cat.ID[,i,][,these,drop=FALSE])
	    Consensus <- MatchPedField(Offspring=Cat.off, Cat.i)
	    Likely.L[[i]] <- data.frame("Cand.id" = rownames(Cat.i),
	                                Cat.i,
	                                Off = Cat.off,
	                                L.off = LL.off[i, these],
	                                LR.off = LL.off[i,these] - Ref.i[i],
	                                OK = Consensus,
	                                stringsAsFactors = FALSE)
		}
		Likely.df <- ldply(Likely.L, .id="id.orig")
    Likely.df$id.orig <- as.character(Likely.df$id.orig)
    Likely.df <- Likely.df[Likely.df$Cand.id != "REF", ]
    Likely.df$id <- Pedigree$id[match(Likely.df$id, rownames(Pedigree))]

    #~~~ check & flag issues ~~~
    warn.msg <- c("MISMATCH" = "Some mismatches, please check",
                  "MERGE?" = "Some ids possibly to be merged, please check",
                  "CHECKMUM" = "Some ids have pedigree mum but no field mum, or mismatch",
                  "CHECKYEARS" = "Some birth and/or death years are incompatible")
    for (NotOK in c("MISMATCH", "MERGE?", "CHECKMUM", "CHECKYEARS")) {
      if (any(Likely.df$OK == NotOK)) {
        warning(warn.msg[NotOK], immediate.=TRUE)
        print(Likely.df[which(Likely.df$OK == NotOK), ], row.names = FALSE)
        cat("\n")
      }
    }

    Match.df <- Likely.df[Likely.df$OK %in% OK.accept[[r]] &
                                Likely.df$Cand.id != "INVIS", ]

    # guarantee that each Cand.id matched to at most one id.
    if (any(duplicated(Match.df$Cand.id))) {
      cand.dup <- Match.df$Cand.id[duplicated(Match.df$Cand.id)]
      Likely.df$OK[Likely.df$Cand.id %in% cand.dup & Likely.df$OK %in% c("OK", "LIKELY")] <- "MULTI"
      Match.df <- Likely.df[Likely.df$OK %in% OK.accept[[r]] &
                                Likely.df$Cand.id != "INVIS", ]
      warning("Some candidates match >1 pedigree id", immediate.=TRUE)
      tmp <- Likely.df[Likely.df$OK %in% "MULTI", ]
      print(tmp[order(tmp$Cand.id), ], row.names=FALSE)
    }

    Likely.r[[r]] <- Likely.df
    Match.r[[r]] <- Match.df

    utils::flush.console()

    print(table(ifelse(rownames(Pedigree) %in% FieldIDs, "real",
                       ifelse(Pedigree$id %in% FieldIDs, "renamed", "dummy/mystery")),
                "Matched" = Pedigree$id %in% Match.df$id))

    #~~~
    Replace.OK <- Match.df[!Match.df$id %in% FieldIDs, ]

    Pedigree$id[match(Replace.OK$id, Pedigree$id)] <- Replace.OK$Cand.id
    for (x in c("dam", "sire")) {
      these <- which(Replace.OK$id %in% Pedigree[, x])
      for (i in these) {
        Pedigree[, x] <- gsub(Replace.OK$id[i], Replace.OK$Cand.id[i],
                              Pedigree[, x])
      }
    }

    # conditional continue
    if (nag & any(Likely.df$OK %in% c("MISMATCH", "MERGE?", "MULTI")) & r < nIter) {
      ANS <- readline(prompt = "\n\n Continue with next iteration Y/N? ")
      if (!substr(ANS, 1, 1) %in% c("Y", "y", "J", "j", "")) {
        assign("Match.dump", list(LH.field = LH.field,
                                  LH.pedigree = LH.pedigree,
                                  Cat.ID = Cat.ID,
                                  cand.dams = cand.dams,
                                  cand.sires = cand.sires,
                                  LL.off = LL.off,
                                  Likely.r = Likely.r,
                                  Match.r = Match.r),
               envir = .GlobalEnv)
        stop("Dumping data in 'Match.dump', and stopping")
      }
    }
  }

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  MkType <- function(name) ifelse(is.na(name),
                                  NA,
                                  ifelse(name %in% intersect(rownames(Pedigree), FieldIDs),
                                         "G",  # genotyped, non-replaced ID
                                         ifelse(name %in% FieldIDs,
                                                "R",    # pedigree dummy/mystery ID replaced by field ID
                                                "D")))  # dummy or mystery ID
  for (x in c("id", "dam", "sire")) {
    Pedigree[, paste0("type.", x)] <- MkType(Pedigree[,x])
  }
  return( list(Pedigree = Pedigree,
               Likely.r = Likely.r,
               Match.r = Match.r) )
}



###############################################################################
# global variables keeps giving a NOTE during CHECK, not sure how to fix
if(getRversion() >= "2.15.1") utils::globalVariables(c("PedFieldMatch", "Match.dump"))
if(getRversion() >= "3.1.0") utils::suppressForeignCheck(c("PedFieldMatch", "Match.dump"))
###############################################################################
