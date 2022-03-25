#' @title Calculate potential reproductive years
#'
#' @description Calculate minimum year of first reproduction and maximum year of
#'   last reproduction, based on individual birth year, death year and month,
#'   and population age of first and last reproduction.
#'
#' @details If the breeding season is in September to November and offspring are
#'   born the next spring, than a male which died in September of year \emph{t} may
#'   still have sired offspring born in year \emph{t+1}: \code{YearLastRepro = 1} and
#'   \code{MonthLastRepro = 9}. A male which died in August cannot be the sire
#'   of any offspring born in year \emph{t+1}.
#'
#' @param LHextra datAgeFirstReproame with lifehistory data, with columns id - Sex -
#'   BirthYear - DeathYear - DeathMonth. Missing values allowed.
#' @param AgeFirstRepro vector of length 2, population minimum age of first
#'   reproduction for females and males, respectively.
#' @param AgeLastRepro vector of length 2, population maximum age of last
#'   reproduction for females and males, respectively.
#' @param YearLastRepro vector of length 2, maximum calendar years of last
#'   reproduction after calendar year of death.
#' @param MonthLastRepro vector of length 2, minimum calendar month of last
#'   reproduction in calendar year of death. Note for females that foetusses
#'   may be sampled.
#'
#' @return \code{LHextra} with 2 added columns:
#' \describe{
#'    \item{FirstRepro}{The (minimum) year in which the individual could first
#'    have reproduced. NA for individuals that died before the age of first
#'    reproduction.}
#'    \item{LastRepro}{The (maximum) year in which the individual could last
#'    have reproduced}}
#'
#' @seealso \code{\link{CheckBirthYears}}
#'
#' @examples
#' \dontrun{
#' LHextra <- CalcFirstLastRepro(LHextra,
#'                              AgeFirstRepro = c(3, 2),
#'                              AgeLastRepro = c(25, 20),
#'                              MonthLastRepro = c(1, 9),
#'                              YearLastRepro = c(0, 1))
#' LH.chk <- CheckBirthYears(Pedigree,
#'                          LHextra,
#'                          AgeFirstRepro = c(3, 2),
#'                          AgeLastRepro = c(25,20))
#' }
#'
#' @export

CalcFirstLastRepro <- function(LHextra,
                               AgeFirstRepro = c(1, 1),
                               AgeLastRepro = c(100, 100),
                               YearLastRepro = c(0, 0),
                               MonthLastRepro = c(1, 1))
{
  if (!"FirstRepro" %in% names(LHextra))  LHextra$FirstRepro <- NA
  if (!"LastRepro" %in% names(LHextra))   LHextra$LastRepro <- NA

  # values for unknown sex
  LHextra$Sex[!LHextra$Sex %in% 1:3] <- 3
  AgeFirstRepro[3] <- min(AgeFirstRepro)
  AgeLastRepro[3] <- max(AgeLastRepro)
  YearLastRepro[3] <- YearLastRepro[2]
  MonthLastRepro[3] <- MonthLastRepro[2]

  LHextra$FirstRepro <- with(LHextra,
                             ifelse(!is.na(FirstRepro), FirstRepro,
                                ifelse(is.na(BirthYear),
                                       NA,
                                       BirthYear + AgeFirstRepro[Sex])))

  LHextra$LastRepro <- with(LHextra,
                            ifelse(!is.na(LastRepro), LastRepro,
                               ifelse(is.na(DeathYear),
                                  ifelse(is.na(BirthYear),
                                         NA,
                                         BirthYear + AgeLastRepro[Sex]),
                                  ifelse(is.na(DeathMonth) | DeathMonth >= MonthLastRepro[Sex],
                                     DeathYear + YearLastRepro[Sex],  # could have reproduced in year of death
                                     DeathYear + YearLastRepro[Sex] -1))))  # could not have

  LHextra$FirstRepro[which(LHextra$FirstRepro > LHextra$LastRepro)] <- NA  # died before age first repro

  return( LHextra )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Check consistency of birth years and potential reproductive years
#'   throughout pedigree
#'
#' @description Check that offspring birth years are consistent with their
#'   assigned parents' potential reproductive years as calculated by
#'   \code{\link{CalcFirstLastRepro}} based on their birth year, death year and
#'   month, and population age of first and last reproduction.
#'
#' @param Pedigree dataframe with columns id - dam - sire
#' @param LHextra datframe with lifehistory data, with columns id - Sex -
#'   BirthYear - FirstRepro - LastRepro.
#' @param AgeFirstRepro vector of length 2, population minimum age of first
#'   reproduction for females and males, respectively.
#' @param AgeLastRepro vector of length 2, population maximum age of last
#'   reproduction for females and males, respectively.
#'
#' @return a dataframe with the following columns
#'   \item{id}{}
#'   \item{dam}{}
#'   \item{sire}{}
#'   \item{Sex}{1=females,2=males,3=unknown. Pedigree dam/sire takes precedent
#'     over LHextra}
#'   \item{BirthYear.self, FirstRepro.self, LastRepro.self}{based on own data in
#'     LHextra}
#'   \item{BY.min.dam, BY.max.dam}{based on dam FirstRepro & LastRepro}
#'   \item{BY.min.sire, BY.max.sire}{based on sire FirstRepro & LastRepro}
#'   \item{FirstRepro.off, LastRepro.off}{First and last birth year of pedigree
#'     offspring}
#'   \item{BY.min.off, BY.max.off}{based on offspring birth years & AgeFirstRepro}
#'   \item{OK.dam, OK.sire, OK.off}{dam, sire, offspring lifehistory data
#'     consistent with own data (logical)}
#'   \item{OK.dam.off, OK.sire.off}{dam, sire data consistent with offspring
#'     lifehistory data (logical)}
#'   \item{BY.min, BY.max}{Consensus minimum & maximum birth year}
#'   \item{OK}{everything OK/not (logical)}
#'
#' @examples
#' \dontrun{
#' LHextra <- CalcFirstLastRepro(LHextra,
#'                              AgeFirstRepro = c(3, 2),
#'                              AgeLastRepro = c(25, 20),
#'                              MonthLastRepro = c(1, 9),
#'                              YearLastRepro = c(0, 1))
#' LH.chk <- CheckBirthYears(Pedigree,
#'                          LHextra,
#'                          AgeFirstRepro = c(3, 2),
#'                          AgeLastRepro = c(25,20))
#' }
#'
#' @export

CheckBirthYears <- function(Pedigree,
                            LHextra,
                            AgeFirstRepro = c(1, 1),
                            AgeLastRepro = c(100, 100))
{

  YYY <- c("BirthYear", "FirstRepro", "LastRepro")
  for (x in c("id", "Sex", YYY)) {
    if (!x %in% colnames(LHextra))  stop("LHextra must have column ", x)
  }

  Pedigree <- unique(Pedigree)
  LHextra <- unique(LHextra)

  LHchk <- merge(Pedigree, setNames(LHextra[, c("id", "Sex", YYY)],
                                    c("id", "Sex", paste0(YYY, ".self"))),
                 all.x = TRUE)
  LHchk <- merge(LHchk, setNames(LHextra[, c("id", YYY[2:3])],
                                    c("dam", "BY.min.dam", "BY.max.dam")),
                 all.x = TRUE)
  LHchk <- merge(LHchk, setNames(LHextra[, c("id", YYY[2:3])],
                                    c("sire", "BY.min.sire", "BY.max.sire")),
                 all.x = TRUE)

  LHchk$Sex <- ifelse(LHchk$id %in% Pedigree$dam, 1,   # for dummy parents
                    ifelse(LHchk$id %in% Pedigree$sire, 2,
                       ifelse(!is.na(LHchk$Sex), LHchk$Sex, 3)))

  # Offspring birth years
  OffIDs <- c(dlply(Pedigree, "dam", function(df) df$id),
              dlply(Pedigree, "sire", function(df) df$id))
  OffIDs <- OffIDs[names(OffIDs) != "NA"]

  Off.BY <- ldply(OffIDs, function(O) {
    data.frame(FirstRepro.off = Min(LHextra$BirthYear[LHextra$id %in% O]),
               LastRepro.off = Max(LHextra$BirthYear[LHextra$id %in% O])) },
    .id = "id")

  LHchk <- merge(LHchk, Off.BY, all.x=TRUE)
  AgeFirstRepro[3] <- min(AgeFirstRepro)
  LHchk$BY.max.off <- LHchk$FirstRepro.off - AgeFirstRepro[LHchk$Sex]
  LHchk$BY.min.off <- LHchk$LastRepro.off - AgeLastRepro[LHchk$Sex]

  #~~~~~~~~~~~
  LHchk$OK.dam <- with(LHchk, ifelse(BirthYear.self >= BY.min.dam & BirthYear.self <= BY.max.dam, TRUE, FALSE))
  LHchk$OK.sire <- with(LHchk, ifelse(BirthYear.self >= BY.min.sire & BirthYear.self <= BY.max.sire, TRUE, FALSE))
  LHchk$OK.off <- with(LHchk, ifelse(FirstRepro.self <= FirstRepro.off & LastRepro.self >= LastRepro.off, TRUE, FALSE))

  LHchk$OK.dam.off <- with(LHchk, ifelse(!is.na(BirthYear.self), OK.dam & OK.off,
                                      ifelse(BY.max.off >= BY.min.dam, TRUE, FALSE)))
  LHchk$OK.sire.off <- with(LHchk, ifelse(!is.na(BirthYear.self), OK.sire & OK.off,
                                      ifelse(BY.max.off >= BY.min.sire, TRUE, FALSE)))

  LHchk$BY.min <- with(LHchk, ifelse(!is.na(BirthYear.self), BirthYear.self,
                                     pmax(BY.min.dam, BY.min.sire, BY.min.off, na.rm=TRUE)))
  LHchk$BY.max <- with(LHchk, ifelse(!is.na(BirthYear.self), BirthYear.self,
                                     pmin(BY.max.dam, BY.max.sire, BY.max.off, na.rm=TRUE)))
  LHchk$OK <- with(LHchk, OK.dam & OK.sire & OK.off & OK.dam.off & OK.sire.off & BY.min <= BY.max)
  LHchk$OK[is.na(LHchk$OK)] <- TRUE

  if (any(!LHchk$OK))  warning(paste("There are ", sum(!LHchk$OK), " cases of birthyear discrepancies"), immediate. = TRUE)

  return(LHchk)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Calculate parental age
#'
#' @description wrapper to calculate age of candidate field parent
#'
#' @details the default factor levels are Rum red deer specific
#'
#' @param candpar.df dataframe with columns id - par.ped - par.field
#' @param LHextra dataframe with (at least) columns id - BirthYear - LastRepro
#' @param Factor logical, return ages as a factor using the specified
#'   breakpoints (TRUE) or as a numeric vector (FALSE)
#' @param BRKS breaks to pass to \code{\link{cut}}
#' @param LBLS labels to pass to \code{\link{cut}}
#' @param REFLVL  reference level for factor
#' @param Dead value in case the id's birth year is later than the field
#'   parent's last year of reproduction
#'
#' @return a numeric vector or factor of length equal to number of rows in
#'   candpar.df, and in that exact same order.
#'
#' @examples
#' \dontrun{
#' cand.sires <- Init.candpar(unique(Pedigree[,c("id", "sire")]), FieldSires,
#' 			                             expl.var=colnames(FieldSires)[3])
#' cand.sires$Age.f <- CalcParAge(cand.sires, LHextra, Factor=TRUE)
#' cand.sires$prob <- predict(DadMod, newdata=cand.sires, type="response")
#' }
#'
#' @export

CalcParAge <- function(candpar.df, LHextra, Factor=TRUE,
                       BRKS=c(0, 4, 5,6, 11,12, 14, 25),
                       LBLS=c(3.5, 5, 6, 9, 12, 13.5, 15),
                       REFLVL="9",
                       Dead=-42)
{
  for (x in c("id", "par.ped", "par.field"))  candpar.df[,x] <- as.character(candpar.df[,x])
  if (nrow(LHextra) > length(unique(LHextra$id)))  stop("LHextra must not have duplicate ids")
  LHextra[,"id"] <- as.character(LHextra[,"id"])
  candpar.df$order <- 1:nrow(candpar.df)            # merge(,sort=FALSE)  doesn't work.
  candpar.LH <- merge(candpar.df, LHextra[, c("id", "BirthYear")], all.x=TRUE)
  candpar.LH <- merge(candpar.LH, setNames(LHextra[, c("id", "BirthYear", "LastRepro")],
                                           c("par.field", "BY.par.field", "LastRepro.par.field")),
                      all.x=TRUE)
  candpar.LH <- candpar.LH[order(candpar.LH$order),]
  age.par.field <- with(candpar.LH, BirthYear - BY.par.field)
  age.par.field[which(candpar.LH$BirthYear > candpar.LH$LastRepro.par.field)] <- Dead
	if (!Factor) {
		return( age.par.field )
	} else {
		Age.f <- cut(age.par.field, breaks = BRKS, right=TRUE, labels = LBLS)
		Age.f <- stats::relevel(as.factor(Age.f), ref=REFLVL)
		Age.f[is.na(age.par.field)] <- REFLVL
		return( Age.f )
	}
}
