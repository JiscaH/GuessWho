#' @title Are pedigree and field ID synonymous?
#'
#' @description Match IDs in the pedigree, including dummy and mystery IDs, to
#'   field IDs based on sample label, sex, birth year and reproductive years,
#'   and field-observed offspring of pedigree mother.
#'
#' @param  LH.pedigree genetically inferred pedigree, and 'pedigree-assisted'
#'   Sex, BY.min (minimum birth year), BY.max (maximum birth year), see
#'   \code{\link{CalcFirstLastRepro}}. May include dummy individuals and mystery
#'   samples.
#' @param  LH.field dataframe with lifehistory data based solely on field
#'   observations, with columns id - dam.field - Sex - BY.min - BY.max. Must
#'   have one row for each 'FieldID' (i.e. non-mystery, non-dummy IDs), even if
#'   there is no associated lifehistory data known.
#' @param ids.snpd.field  character vector with IDs that are both SNP genotyped
#'   and proper field IDs; necessary after some dummy/mystery IDs in pedigree
#'   have been replaced.
#'
#' @return An categorical array with 3 dimensions: Label-Sex-Years-
#'   Mother by no. individuals in LH.pedigree by the no. individuals in
#'   LH.field +1 ('INVIS'). The following categories are distinguished:
#' \item{Label}{Match - MismatchSNPd - MismatchNotSNPd - NotSNPd}
#' \item{Sex}{Match - Mismatch - Unknown (either one unknown sex)}
#' \item{Years}{Match - Mismatch - Unknown}
#' \item{Mother}{Match - Mismatch - NewMum (only pedigree dam) - Unverified
#' (only field dam) - Unknown (dummy or mystery dam in pedigree) - None}
#'
#' @seealso \code{\link{MkOffCat}, \link{MatchPedField}}
#'
#' @examples
#' \dontrun{
#' LH.chk <- CheckBirthYears(Pedigree, LHextra, AgeFirstRepro = AgeFirstRepro,
#'                           AgeLastRepro = AgeLastRepro)
#' LH.pedigree <- LH.chk[match(Pedigree$id, LH.chk$id), ]
#' LH.field <- merge(LHextra, FieldDams, all=TRUE)
#' LH.field <- LH.field[match(FieldIDs, LH.field$id), ]
#' Cat.ID <- MkIDcat(LH.pedigree, LH.field,
#'                   ids.snpd.field = intersect(rownames(pedigree), FieldIDs))
#' }
#'
#' @export

MkIDcat <- function(LH.pedigree = NULL,   # id-dam-sire-Sex-BY.min-BY.max
                    LH.field = NULL,      # id-dam.field-Sex-BY.min-BY.max
                    ids.snpd.field = NULL)
{

  #~~~  Check input  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (is.null(LH.pedigree))  stop("Please provide LH.pedigree")
  if (is.null(LH.field))  stop("Please provide LH.field")
  if (is.null(ids.snpd.field))  ids.snpd.field <- intersect(LH.pedigree$id, LH.field$id)


  #~~~  Initiate  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	LH.field <- rbind(LH.field,
	                  merge(LH.field[0,], data.frame(id="INVIS"), all=TRUE))
  cat.id <- array(dim = c(4, nrow(LH.pedigree), nrow(LH.field)),
                   dimnames = list(c("Label", "Sex", "Years", "Mother"),
                                   LH.pedigree$id, LH.field$id))
  # can become very large -- alternative needed?
  LVLS <- vector("list", length=4)

  #~~~  Label  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # own sample label

  labelmatch <- function(idp, idf, ids.snpd.field = ids.snpd.field) {
    if (length(idp)==1) {
      idp <- rep(idp, length(idf))
    } else if(length(idf)==1) {
      idf <- rep(idf, length(idp))
    } else {
      stop("labelmatch doesn't work with 2 vectors")
    }
    ifelse(idp %in% ids.snpd.field,
           ifelse(idp == idf,
                  1,
                  ifelse(idf %in% ids.snpd.field,
                         2,
                         3)),
           ifelse(idf %in% ids.snpd.field,
                  3,
                  4))
  }
  LVLS[["Label"]] <- c("1"="Match", "2"="MismatchSNPd", "3"="MismatchNotSNPd", "4"="NotSNPd")

	for (i in seq_along(LH.pedigree$id)) {   # quicker than outer()
	  cat.id["Label", i, ] <- labelmatch(LH.pedigree$id[i], LH.field$id, ids.snpd.field)
	}
	# character:  18.96 sec
	# numbers:  3.19 sec

  #~~~  Sex  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Match Sex (0/1/0.5)

  sexmatch <- matrix(c(1,2,3,
                       2,1,3,
                       3,3,3), 3, 3)
  LVLS[["Sex"]] <- c("1"="Match", "2"="Mismatch", "3"="Unknown")

  LH.pedigree$Sex[is.na(LH.pedigree$Sex)] <- 3
  LH.field$Sex[is.na(LH.field$Sex)] <- 3

  for (i in seq_along(LH.pedigree$id)) {
    cat.id["Sex", i, ] <- sexmatch[cbind(LH.pedigree$Sex[i], LH.field$Sex)]
  }
  cat.id["Sex", , "INVIS"] <- 1


  #~~~  Years  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Match birth & reproductive years (based on pedigree parents & offspring for dummies)

  for (i in seq_along(LH.pedigree$id)) {
    cat.id["Years", i, ] <- LH.pedigree$BY.min[i] <= LH.field$BirthYear &
                              LH.pedigree$BY.max[i] >= LH.field$BirthYear &
#                              (is.na(LH.field$FirstRepro) | is.na(LH.pedigree$FirstRepro.off[i]) |
                                  LH.pedigree$FirstRepro.off[i] >= LH.field$FirstRepro &
#															(is.na(LH.field$LastRepro) | is.na(LH.pedigree$LastRepro.off[i]) |
                              LH.pedigree$LastRepro.off[i] <= LH.field$LastRepro
    # TODO? add DeathYear?
  }
  cat.id["Years", , ] <- 2-cat.id["Years", , ]   # 0=FALSE, 1=TRUE --> 1=TRUE, 2=FALSE
  cat.id["Years", , "INVIS"] <- 1
  cat.id["Years", , ][is.na(cat.id["Years", , ])] <- 3      # TODO: check: more NA's? how handled?

  LVLS[["Years"]] <- c("1"="Match", "2"="Mismatch", "3"="Unknown")

  #~~~  Mother  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # based on field offspring of pedigree mother

  mummatch <- function(i, j) {
    ifelse(!is.na(LH.field$par.field[j]),
         ifelse(!is.na(LH.pedigree$dam[i]),
            ifelse(LH.field$par.field[j] == LH.pedigree$dam[i],
               1,  # Match
               ifelse(LH.pedigree$dam[i] %in% LH.field$id,
                  2,   # mismatch
                  ifelse(LH.field$par.field[j] %in% LH.pedigree$id,
                      2,
                      5))),   # dummy or mystery dam in pedigree
            ifelse(LH.field$par.field[j] %in% ids.snpd.field &
                     LH.field$id[j] %in% ids.snpd.field,
               2,
               4)),  # Unverified
         ifelse(!is.na(LH.pedigree$dam[i]),
            3,   # new mum
            6))  # no mum
  }

  cat.id["Mother", ,] <- outer(seq_along(LH.pedigree$id), seq_along(LH.field$id), mummatch)
  cat.id["Mother", , "INVIS"] <- 4

  LVLS[["Mother"]] <- c("Match", "Mismatch", "NewMum", "Unverified", "Unknown", "None")

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  cat.id.num <- cat.id
  cat.id[,,] <- NA
  for (x in c("Label", "Sex", "Years", "Mother")) {
    cat.id[x,,] <-  LVLS[[x]][cat.id.num[x,,]]
  }

  return( cat.id )
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Offspring-based categories
#'
#' @description Turn pedigree vs. field parent probabilities into categorical
#'   variables
#'
#' @param LLO named numerical vector for a single pedigree ID with
#'   log10-probabilities to match a subset of field IDs
#' @param REF reference log10-probability, minimum probability required for
#'   'Conclusive'.
#' @param T.off  threshold log10-probability difference between most-likely and
#'   next-most-likely candidate. see 'Value'.
#'
#' @return A character vector of the same length as 'LLO', with the following
#'   possible values:
#'   \item{None}{No offspring in the pedigree}
#'   \item{Conclusive}{This field id is more likely than the next best
#'     by a margin of at least 'T.off', and is at least 'REF'.}
#'   \item{Inferior}{A different field id is more likely, by a margin of at
#'   least 'T.off'}
#'  \item{Inconclusive}{The likelihoods of the best and next-best candidate
#'  differ by less than 'T.off', or the most-likely is less than 'REF'}
#'  \item{Invisible}{'INVIS' is conclusively the most likely candidate}
#'
#' @seealso \code{\link{MkIDcat}, {MatchPedField}}
#'
#' @examples
#' \dontrun{
#' for (i in 1:dim(Cat.ID)[2]) {
#'   these <- which(Cat.ID["Label",i,] == "Match" |
#'                  Cat.ID["Mother",i, ] == "Match" | !is.na(LL.off[i,]))
#'	 if (length(these)>0) {
#'	   Cat.off <- MkOffCat(LL.off[i,][these, drop=FALSE],  T.off=3)
#'		 ...
#'   }}
#' }
#'
#' @export

MkOffCat <- function(LLO, REF, T.off=3) {
	tol = .Machine$double.eps
	AnyOff <- any(!is.na(LLO))
	Cat.off <- ifelse(rep(!AnyOff, length(LLO)),
                 "None",
                 ifelse(is.na(LLO),
                    "Inconclusive",
                    ifelse(LLO -REF > tol,
                       ifelse((LLO - Max2(LLO, m=-Inf)) >= T.off,
                          ifelse(names(LLO)=="INVIS",
                            "Invisible",
                            "Conclusive"),
                          "Inconclusive"),
                       "Inconclusive")))
#
# 	if (any(is.na(Cat.off))) {
# 	  print(LLO)
# 	  cat('\n')
# 	  print(Cat.off)
# 	  stop('Cat.off NA!')
# 	}

	if (any(Cat.off == "Conclusive")) {
		Cat.off[Cat.off == "Inconclusive"] <- "Inferior"
	} else if (any(Cat.off == "Invisible")) {
		Cat.off[Cat.off == "Inconclusive"] <- "Invisible"
	}
	return( Cat.off )
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Consensus match
#'
#' @description Consensus category if the pedigree id matches any of the
#'   shortlisted candidate field ids, based on label, sex, birth years, pedigree
#'   mother's field offspring and pedigree offspring's field candidate parents.
#'
#' @param Offspring vector with categories describing the match based on
#'   pedigree offspring; output of \code{\link{MkOffCat}}
#' @param Cat.i matrix with 4 columns (Label - Sex - Years - Mother) and one row
#'   per shortlisted candidate field id (optionally including INVIS); subset of
#'   output of \code{\link{MkIDcat}}
#'
#' @return a vector with length equal to the number of field ids (number of rows
#'   in 'DF'), taking values
#'   \item{OK}{match}
#'   \item{X}{not a match}
#'   \item{LIKELY}{Conclusive based on offspring, not SNPd, but undetermined if
#'   dummy mother matches}
#'   \item{MAYBE}{maternal sibling or match: Mother matches, not SNPd,
#'   no/inconclusive based on offspring, Sex & Years match/unknown, }
#'   \item{MERGE?}{one id plausible based on label, and another based on
#'   offspring}
#'   \item{CHECKMUM}{id plausible based on label/offspring, but mum in pedigree
#'   & not in field, or mismatches}
#'   \item{CHECKYEARS}{id plausible based on label/offspring, but inconsistency
#'   in birth/reproductive years}
#'   \item{CHECKSEX}{id plausible based on label/offspring, but sex mismatch}
#'   \item{???}{Not enough info to assess plausibility}
#'
#' @seealso \code{\link{MkIDcat}, \link{MkOffCat}, \link{init.PedFieldMatch}}
#'
#' @examples
#' \dontrun{
#' Likely.L <- setNames(vector("list", length=dim(Cat.ID)[2]), unique(Pedigree$id))
#' for (i in 1:dim(Cat.ID)[2]) {
#'   these <- which(Cat.ID["Label",i,] == "Match" |
#'                  Cat.ID["Mother",i, ] == "Match" |
#'                  !is.na(LL.off[i,]))
#'   if (length(these)==0)  next
#'   Cat.off <- MkOffCat(LL.off[i,][these, drop=FALSE],  REF=Ref.i[i], T.off=3)
#'   Cat.i <- t(Cat.ID[,i,][,these,drop=FALSE])
#'   Consensus <- MatchPedField(Offspring=Cat.off, Cat.i)
#'   Likely.L[[i]] <- data.frame("Cand.id" = rownames(Cat.i),
#'                               Cat.i,
#'                               Off = Cat.off,
#'                               L.off = LL.off[i, these],
#'                               LR.off = LL.off[i,these] - Ref.i[i],
#'                               OK = Consensus,
#'                               stringsAsFactors = FALSE)
#' }
#' }
#'
#' @export

MatchPedField <- function(Offspring, Cat.i) {
	if (!exists("PedFieldMatch"))  init.PedFieldMatch()

  if (nrow(Cat.i)==0) {
    OK <- character(0)
  } else if (nrow(Cat.i) == 1) {
    OK <- ifelse(Cat.i[1,"Label"] == "Match",
       "OK",
       "???")
  } else {

		OK <- PedFieldMatch[cbind(Offspring, Cat.i)]
		names(OK) <- rownames(Cat.i)

    if(sum(OK %in% c("OK", "LIKELY"), na.rm=TRUE) > 1) {
      stop(">1 OK/LIKELY")
    } else if(sum(OK == "OK", na.rm=TRUE) == 1) {
      OK[OK %in% c("MAYBE", "UNLIKELY", "MISMATCH")] <- "X"
    }
    if (sum(OK == "MAYBE", na.rm=TRUE) == 1) {
      OK[OK == "MAYBE"] <- "LIKELY"
    }
		if (any(OK=="MERGE?") && OK["INVIS"]=="MERGE?") { # !!
		  OK[OK=="MERGE?"] <- "OK"
		  OK["INVIS"] <- "X"
		}
  }

	return( OK )
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Initiate array 'PedFieldMatch'
#'
#' @description This array is used to 'translate' the various levels in the
#'   array returned by \code{\link{MkIDcat}}, plus offspring-based likelihoods
#'   being conclusive/inconclusive or none, for no offspring, into a consensus.
#'
#' @details For the Offspring-based category, 'Invisible' denotes that an
#'   unobserved individual is the most likely match. Note that this may be one
#'   of the (other, shortlisted) candidate field parents, if they were not
#'   observed during the relevant breeding period(s).
#'
#' @return a 5D-array is created globally, taking the following values:
#'   \item{OK}{match}
#'   \item{X}{not a match}
#'   \item{LIKELY}{Conclusive based on offspring, not SNPd, but undetermined if
#'   dummy mother matches}
#'   \item{MAYBE}{maternal sibling or match: Mother matches, not SNPd,
#'   no/inconclusive based on offspring, Sex & Years match/unknown, }
#'   \item{MERGE?}{one id plausible based on label, and another based on
#'   offspring}
#'   \item{CHECKMUM}{id plausible based on label/offspring, but mum in pedigree
#'   & not in field, or mismatches}
#'   \item{CHECKYEARS}{id plausible based on label/offspring, but inconsistency
#'   in birth/reproductive years}
#'   \item{CHECKSEX}{id plausible based on label/offspring, but sex mismatch}
#'   \item{???}{Not enough info to assess plausibility}
#'
#' Dimensions 2-5 are as returned by \code{\link{MkIDcat}}:
#' \item{Label}{Match - MismatchSNPd - MismatchNotSNPd - NotSNPd; 'SNPd' and
#'   'NotSNPd' refers to the field ID}
#' \item{Sex}{Match - Mismatch - Unknown}
#' \item{Years}{Match - Mismatch - Unknown}
#' \item{Mother}{Match - Mismatch - NewMum (only pedigree dam) - Unverified
#' (only field dam) - Unknown (dummy or mystery dam in pedigree) - None}
#' \item{Offspring}{Conclusive (is best match) - Inconclusive - Inferior (different candidate is best match) - Invisible ('INVIS' is best match) - None}
#'
#' @seealso \code{\link{MatchPedField}}
#'
#' @export

init.PedFieldMatch <- function() {
		PedFieldMatch <- array(dim=c(5, 4, 3, 3, 6),
													 dimnames = list(Offspring = c("Conclusive", "Inconclusive", "Inferior", "Invisible", "None"),
																					 Label = c("Match", "MismatchSNPd", "MismatchNotSNPd", "NotSNPd"),
																					 Sex = c("Match", "Mismatch", "Unknown"),
																					 Years = c("Match", "Mismatch", "Unknown"),
																					 Mother = c("Match", "Mismatch", "NewMum", "Unverified", "Unknown", "None")
																					 )
													 )
		PedFieldMatch[c("Conclusive","None"), "Match", c("Match", "Unknown"), c("Match","Unknown"), c( "Match", "Unknown", "None")] <- "OK"
		PedFieldMatch[c("Conclusive","None"), "Match", c("Match", "Unknown"), c("Match","Unknown"), "Unverified"] <- "LIKELY"

		PedFieldMatch[c("Inconclusive","Invisible"), "Match", , , ] <- "OK"
		PedFieldMatch["Conclusive", "NotSNPd", , , c("Match", "None")] <- "OK"
		PedFieldMatch["Conclusive", "NotSNPd", , , c("Unknown","Unverified")] <- "LIKELY"  # dummy mum

		PedFieldMatch[c("Inconclusive", "None", "Invisible"),"NotSNPd",c("Match","Unknown"),c("Match","Unknown"),"Match"] <- "MAYBE"  # match or mat sib
		PedFieldMatch[c("Inconclusive", "None", "Invisible"),"NotSNPd",c("Match","Unknown"),c("Match","Unknown"),c("Unknown","None","NewMum", "Unverified")] <- "???"

		PedFieldMatch["Inferior",,,,] <- "X"  # Offspring
		PedFieldMatch[,c("MismatchSNPd", "MismatchNotSNPd"),,,] <- "X"  # label
		PedFieldMatch[,,"Mismatch",,] <- "X"  # Sex
		PedFieldMatch[,,,"Mismatch",] <- "X"  # Years
		PedFieldMatch[,,,,"Mismatch"] <- "X"  # Mother

		PedFieldMatch[,"Match",,"Mismatch",] <- "CHECKYEARS"
		PedFieldMatch["Conclusive",,,"Mismatch",] <- "CHECKYEARS"
		PedFieldMatch[,"Match","Mismatch",,] <- "CHECKSEX"
		PedFieldMatch["Conclusive",,"Mismatch",,] <- "CHECKSEX"
		PedFieldMatch[,"Match",,,c("Mismatch", "NewMum")] <- "CHECKMUM"
		PedFieldMatch["Conclusive",,,,c("Mismatch", "NewMum")] <- "CHECKMUM"
		PedFieldMatch["Inferior","Match",,,] <- "MERGE?"
		PedFieldMatch["Conclusive",c("MismatchSNPd", "MismatchNotSNPd"),,,] <- "MERGE?"

		PedFieldMatch <<- PedFieldMatch   # make global
}
