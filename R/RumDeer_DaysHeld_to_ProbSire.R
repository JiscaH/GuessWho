#' @title Days held during conception window to sire probability
#'
#' @description Calculate the probability the male is the offspring's sire from
#'   the number of days its mother was in his harem during the 11-day window
#'   around back-calculated conception (birth date - 235 days).
#'
#' @details The sum of Ndays per offspring may exceed 11, as Some females are
#'   recorded with 2 different males on the same day. Calculating part-days for
#'   these cases was tried, but does not notably affect the relationship between
#'   days held and probability to sire the offspring.
#'
#'   The age levels were determined by fitting a model with each age as a single
#'   level, and iteratively combining ages which effects on paternity
#'   probability did not differ significantly.
#'
#' @param ObsSire dataframe with the following columns:
#' \describe{
#'   \item{id}{}
#'   \item{par.ped}{genetically inferred parent; 1 per id, or NA)}
#'   \item{par.field}{Field observed candidate parent; 0, 1 or several per id}
#'   \item{Ndays}{number of days \code{id}'s mother was in \code{par.field}'s harem during
#'     the 11-day conception window.}
#'   \item{age.par.field}{age in years of \code{par.field} at \code{id}'s birth. Only
#'     used when \code{Model = 'age'}} }
#' @param Model binomial glm to calculate the probability that a male is the
#'   father, using column names in DaysHeld. Implemented are
#'   \describe{
#'     \item{\code{'default'}}{fits Ndays and Ndays^2}
#'     \item{\code{'age'}}{additionally fits age as a factor with levels 2-4, 5,
#'   6, 7-11, 12, 13-14, 15+.} }
#' @param Plot plot observed & predicted values?
#'
#' @return the fitted glm model
#'
#' @section Database query:
#' The in-build query "Sys_BehaviouralPaternities" returns candidate fathers in
#' whose harem the mother was seen at least 6 days during the conception window.
#' To obtain all candidate fathers, including those in whose harem the mother
#' was seen 1-5 days, build a 'Totals' query around "sys_ConceptionBirthDates"
#' and "tblRutcharts", with
#' \itemize{
#'   \item sys_ConceptionBirthDates.Code (Group By)
#'   \item sys_ConceptionBirthDates.MumCode linked to tblRutCharts.HindCode (Group By)
#'   \item tblRutCharts.Date (Where; do not show)
#'   \item Criteria = Abs([sys_ConceptionBirthDates.ConceptPredict]-[tblRutCharts.Date])<6
#'   \item tblRutCharts.StagCode (Group By)
#'   \item Criteria = Not "NONE"
#'   \item Ndays: tblRutCharts.Date (Count) }
#'
#' @examples
#' \dontrun{
#' candsires.tmp <- Init.candpar(Pedigree[,c("id", "sire")], FieldCandSires, expl.var="Ndays")
#' candsires.tmp$age.par.field <- CalcParAge(candsires.tmp, LHextra, Factor=FALSE)
#' candsires.tmp$Age.f <- CalcParAge(candsires.tmp, LHextra, Factor=TRUE, REFLVL="9")
#'
#' # subset those to run model on. don't use Ndays=0: not all present/alive that
#' # year. Suggest to use only those where both id and par.field are SNP genotyped,
#' # such that it is known with high accuracy whether or not par.field == par.ped.
#' # par.ped may be NA, which is treated as par.field != par.ped.
#' candsires.tmp$Use <- with(candsires.tmp, ifelse(Ndays == 0, FALSE,
#'                                      ifelse(is.na(age.par.field), FALSE,
#'                                         ifelse(par.field %in% Pedigree$id
#'                                                & id %in% Pedigree$id, TRUE, FALSE))))
#' candsires.tmp <- candsires.tmp[candsires.tmp$Use, ]
#'
#' mod.DaysHeld <- DaysHeld_to_Prob(ObsSire = candsires.tmp, Model = "age", Plot=TRUE)
#' }
#'
#' @importFrom stats binomial
#'
#' @export

DaysHeld_to_Prob <- function(ObsSire = NULL,
                             Model = "default",
                             Plot = TRUE)
{
  # input check
  if (!Model %in% c("default", "age"))  stop("Model must be 'default' or 'age'")
  if (!all(c("par.ped", "par.field", "Ndays") %in% colnames(ObsSire)))
    stop("ObsSire must have columns 'par.ped', 'par.field', and 'Ndays'")
  if (Model == "age") {
    if (!"Age.f" %in% colnames(ObsSire))  stop("ObsSire must have column 'Age.f'")
  }

  ObsSire$IsFather <- with(ObsSire,
                           ifelse(is.na(par.ped), 0,   # assume all parents are assigned.
                            ifelse(par.field == par.ped, 1, 0)))

  if (Model == "default") {
    mod.prob <- stats::glm(IsFather~ Ndays + I(Ndays^2), family=binomial(link="probit"),
                    data = ObsSire)

  } else if (Model == "age") {
    mod.prob <- stats::glm(IsFather ~ Ndays + I(Ndays^2) + Age.f, family=binomial(link="probit"),
                    data = ObsSire)
  }

  #~~~~~~~~~~~~~~~~~~
  if (Plot) {
    if (Model == "default") {
      pred.prob <- stats::predict(mod.prob, newdata=data.frame(Ndays=c(0:11)), type="response")
      # round(setNames(pred.prob, 0:11), 2)
      #    0    1    2    3    4    5    6    7    8    9   10   11
      # 0.04 0.11 0.21 0.34 0.46 0.57 0.66 0.72 0.76 0.78 0.79 0.78   # age 5+, SNPd male, SNPd offspring

    } else if (Model == "age") {
      df.pred <- data.frame(Ndays=rep(0:11, times=7),
                      Age.f=as.factor(rep(levels(ObsSire$Age.f), each=12)))
      df.pred$prob <- stats::predict(mod.prob, newdata=df.pred, type="response")  # se.fit=TRUE
    }

    prop.sire <- with(ObsSire, tapply(IsFather, round(Ndays), mean, na.rm=T))
    graphics::plot(1:11, prop.sire, pch=16, col="grey", cex=sqrt(table(round(ObsSire$Ndays)))/5,
         ylim=c(0,1), xlim=c(0,12), xlab="Days held", ylab="Probability to be father", las=1,
         main = "Siring probability vs. days held among genotyped males")
    graphics::legend("topleft", legend=c("Observed proportion", "Model prediction"), pch=c(16,NA), lwd=c(NA,2),
             col=c("grey", 1), pt.cex=c(2,NA), inset=.02)

    if (Model == "default") {
      graphics::lines(0:11, pred.prob, lwd=2)
    } else if (Model == "age") {
      COL <- c("3.5"="seagreen1", "5" = "seagreen3", "6" = "forestgreen", "9" = 1,
               "12" = "chocolate4", "13.5" = "darkgoldenrod", "15" = "chocolate")
      MaxDays <- with(ObsSire, tapply(Ndays, Age.f, max))
      for (a in levels(ObsSire$Age.f)) {
        graphics::lines(0:MaxDays[a], df.pred$prob[df.pred$Age.f == a][0:MaxDays[a]+1],
              lwd=ifelse(a=="9", 3, 1), col=COL[a])
      }
      tbl.Age <- table(ObsSire$Age.f)
      graphics::legend("bottomright",
             legend=paste0(c("2-4", 5,6,"7-11", 12, "13-14", "15-17"), " [", tbl.Age[c(2:4,1,5:7)], "]"),
             col=COL, lwd=2, title="Age [n]", inset=.02)
    }
  }
  #~~~~~~~~~~~~~~~~~~

  return( mod.prob )
}


