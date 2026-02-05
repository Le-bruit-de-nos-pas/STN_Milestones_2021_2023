
library(survival)
library(survminer)
library(lubridate)
library(cmprsk)
library(tidyverse)
library(utile.visuals)
library(Hmisc)
library(corrplot)
library(ggcorrplot)
library(pROC)
library(caret)
library(mgcv)



# Global time-to-death ----------------------
GlobalTimeToDeath=read.csv("GlobalTimeToDeath.csv", sep=";")

Surv(GlobalTimeToDeath$time, GlobalTimeToDeath$status)


ggsurvplot(
  fit = survfit(Surv(time, status) ~ 1, data = GlobalTimeToDeath), 
  xlab = "Follow-up (Months)", 
  ylab = "Overall survival probability",
  legend = "none",
  break.x.by = 12,
  risk.table = TRUE,
  tables.col = "strata")

# Check Plot -> Export (...)

# Probability of surviving X's years
summary(survfit(Surv(time, status) ~ 1, data = GlobalTimeToDeath), times = 12)

summary(survfit(Surv(time, status) ~ 1, data = GlobalTimeToDeath), times = 24)


summary(survfit(Surv(time, status) ~ 1, data = GlobalTimeToDeath), times = 48)

summary(survfit(Surv(time, status) ~ 1, data = GlobalTimeToDeath), times = 96)

# By gender 
ggsurvplot(
  fit = survfit(Surv(time, status) ~ sex, data = GlobalTimeToDeath), 
  xlab = "Follow-up (Months)", 
  ylab = "Overall survival probability",
  legend = "none",
  break.x.by = 12,
  risk.table = TRUE,
  tables.col = "strata",
  conf.int = TRUE)

# Check Plot -> Export (...)

coxph(Surv(time, status) ~ sex, data = GlobalTimeToDeath)






#  Cox proportional-hazards regression modeling to measure the association between baseline variables 
#  and survival of Parkinson's disease patients after Deep brain stimulation surgery 

GlobalDeath = read.csv("GlobalDeath.csv", sep = ";", na.strings=c("", " ","NA"))

# Measuring the contribution of each variable individually (i.e. non-adjusted) -------------

res.cox <- coxph(Surv(time, status) ~ AgeSurgery, data = GlobalDeath)
summary(res.cox)


res.cox <- coxph(Surv(time, status) ~ Sex, data = GlobalDeath)
summary(res.cox)

res.cox <- coxph(Surv(time, status) ~ DiseaseDurationSurgery, data = GlobalDeath)
summary(res.cox)

res.cox <- coxph(Surv(time, status) ~ LEDDpreop, data = GlobalDeath)
summary(res.cox)


res.cox <- coxph(Surv(time, status) ~ LCTresp, data = GlobalDeath)
summary(res.cox)

res.cox <- coxph(Surv(time, status) ~ AgeSurgeryMORE, data = GlobalDeath)
summary(res.cox)


res.cox <- coxph(Surv(time, status) ~ DiseaseDurationSurgeryMORE, data = GlobalDeath)
summary(res.cox)

res.cox <- coxph(Surv(time, status) ~ LCTrespMORE, data = GlobalDeath)
summary(res.cox)

res.cox <- coxph(Surv(time, status) ~ GAITOFFmore, data = GlobalDeath)
summary(res.cox)


res.cox <- coxph(Surv(time, status) ~ PostInstOFFmore, data = GlobalDeath)
summary(res.cox)


res.cox <- coxph(Surv(time, status) ~ HYOFFmore, data = GlobalDeath)
summary(res.cox)

res.cox <- coxph(Surv(time, status) ~ factor(Neuropsychologic, levels = c("Normal", "MCI")), data = GlobalDeath)
summary(res.cox)


res.cox <- coxph(Surv(time, status) ~ UPDRS3off, data = GlobalDeath)
summary(res.cox)


res.cox <- coxph(Surv(time, status) ~ UPDRS3on, data = GlobalDeath)
summary(res.cox)


res.cox <- coxph(Surv(time, status) ~ UPDRS2off, data = GlobalDeath)
summary(res.cox)

res.cox <- coxph(Surv(time, status) ~ UPDRS2on, data = GlobalDeath)
summary(res.cox)


res.cox <- coxph(Surv(time, status) ~ factor(Phenotype, levels = c("Tremor", "PIGD", "Indeterminate")), data = GlobalDeath)
summary(res.cox)

# Adjusting for different variables, Multiple models compared --------------

res.coxGLOBAL <- coxph(Surv(time, status) ~ AgeSurgery + factor(Sex, levels = c("F", "M")) + 
                   DiseaseDurationSurgery +  LEDDpreop + 
                   factor(GAITOFFmore, levels = c("ZeroOROne", "EqualORMore2")) +
                   factor(PostInstOFFmore, levels = c("ZeroOne", "EqualMore2")) +
                   factor(Neuropsychologic, levels = c("Normal", "MCI")) + 
                   UPDRS3off +  UPDRS3on +  UPDRS2off +  UPDRS2on + 
                   factor(Phenotype, levels = c("Tremor", "PIGD", "Indeterminate")), data = GlobalDeath)
summary(res.coxGLOBAL)


res.coxGLOBAL2 <- coxph(Surv(time, status) ~ AgeSurgery +  factor(Sex, levels = c("F", "M")) + 
                         DiseaseDurationSurgery +  LEDDpreop + 
                         factor(GAITOFFmore, levels = c("ZeroOROne", "EqualORMore2")) +
                         factor(PostInstOFFmore, levels = c("ZeroOne", "EqualMore2")) +
                         UPDRS3off +  UPDRS3on +  UPDRS2off +  UPDRS2on + 
                         factor(Phenotype, levels = c("Tremor", "PIGD", "Indeterminate")), data = GlobalDeath)
summary(res.coxGLOBAL2)


res.coxGLOBAL3 <- coxph(Surv(time, status) ~ AgeSurgery +factor(Sex, levels = c("F", "M")) + 
                          DiseaseDurationSurgery +  LEDDpreop +   factor(PostInstOFFmore, levels = c("ZeroOne", "EqualMore2")) +
                          UPDRS3off +   UPDRS3on +  UPDRS2off +  UPDRS2on + 
                          factor(Phenotype, levels = c("Tremor", "PIGD", "Indeterminate")), data = GlobalDeath)
summary(res.coxGLOBAL3)


res.coxGLOBAL4 <- coxph(Surv(time, status) ~ AgeSurgery + DiseaseDurationSurgery + LEDDpreop + 
                          factor(PostInstOFFmore, levels = c("ZeroOne", "EqualMore2")) +
                          UPDRS3off + UPDRS3on +  UPDRS2off +  UPDRS2on + 
                          factor(Phenotype, levels = c("Tremor", "PIGD", "Indeterminate")), data = GlobalDeath)
summary(res.coxGLOBAL4)



res.coxGLOBAL5 <- coxph(Surv(time, status) ~ AgeSurgery + LEDDpreop + 
                          factor(PostInstOFFmore, levels = c("ZeroOne", "EqualMore2")) + UPDRS3off + 
                          UPDRS3on +  UPDRS2off +  UPDRS2on + 
                          factor(Phenotype, levels = c("Tremor", "PIGD", "Indeterminate")), data = GlobalDeath)
summary(res.coxGLOBAL5)



res.coxGLOBAL6 <- coxph(Surv(time, status) ~ AgeSurgery +factor(PostInstOFFmore, levels = c("ZeroOne", "EqualMore2")) +
                          UPDRS3off + UPDRS3on +  UPDRS2off +   UPDRS2on + 
                          factor(Phenotype, levels = c("Tremor", "PIGD", "Indeterminate")), data = GlobalDeath)
summary(res.coxGLOBAL6)



res.coxGLOBAL7 <- coxph(Surv(time, status) ~ AgeSurgery + factor(PostInstOFFmore, levels = c("ZeroOne", "EqualMore2")) +
                          UPDRS3off + UPDRS3on +  UPDRS2off + factor(Phenotype, levels = c("Tremor", "PIGD", "Indeterminate")), data = GlobalDeath)
summary(res.coxGLOBAL7)


res.coxGLOBAL8 <- coxph(Surv(time, status) ~ AgeSurgery + UPDRS3off + UPDRS3on +   UPDRS2off + 
                          factor(Phenotype, levels = c("Tremor", "PIGD", "Indeterminate")), data = GlobalDeath)
summary(res.coxGLOBAL8)



res.coxGLOBAL9 <- coxph(Surv(time, status) ~ AgeSurgery + UPDRS3off + UPDRS3on + 
                          factor(Phenotype, levels = c("Tremor", "PIGD", "Indeterminate")), data = GlobalDeath)
summary(res.coxGLOBAL9)


res.coxGLOBAL10 <- coxph(Surv(time, status) ~ AgeSurgery + UPDRS3on + 
                          factor(Phenotype, levels = c("Tremor", "PIGD", "Indeterminate")), data = GlobalDeath)
summary(res.coxGLOBAL10)




res.coxGLOBAL11 <- coxph(Surv(time, status) ~ AgeSurgery + 
                           factor(Phenotype, levels = c("Tremor", "PIGD", "Indeterminate")), data = GlobalDeath)
summary(res.coxGLOBAL11)



# Cumulative Incidence
# After Scrucca L., Santucci A., Aversa F. (2007) Competing risks analysis using R: 
# an easy guide for clinicians. Bone Marrow Transplantation, 40, 381--387.

"CumIncidence" <- function(ftime, fstatus, group, t, strata, rho = 0, 
                           cencode = 0, subset, na.action = na.omit, level,
                           xlab = "Time", ylab = "Probability", 
                           col, lty, lwd, digits = 4)
{
  # check for the required package
  if(!require("cmprsk"))
  { stop("Package `cmprsk' is required and must be installed.\n 
           See help(install.packages) or write the following command at prompt
           and then follow the instructions:\n
           > install.packages(\"cmprsk\")") } 
  
  mf  <- match.call(expand.dots = FALSE)
  mf[[1]] <- as.name("list")
  mf$t <- mf$digits <- mf$col <- mf$lty <- mf$lwd <- mf$level <- 
    mf$xlab <- mf$ylab <- NULL
  mf <- eval(mf, parent.frame())
  g <- max(1, length(unique(mf$group)))
  s <- length(unique(mf$fstatus))
  if(missing(t)) 
  { time <- pretty(c(0, max(mf$ftime)), 6)
  ttime <- time <- time[time < max(mf$ftime)] }
  else { ttime <- time <- t }
  
  fit   <- do.call("cuminc", mf)
  tfit <- timepoints(fit, time)
  
  cat("\n+", paste(rep("-", 67), collapse=""), "+", sep ="")
  cat("\n| Cumulative incidence function estimates from competing risks data |")
  cat("\n+", paste(rep("-", 67), collapse=""), "+\n", sep ="")
  tests <- NULL
  if(g > 1)
  { 
    tests <- data.frame(fit$Tests[,c(1,3,2)], check.names = FALSE)
    colnames(tests) <- c("Statistic", "df", "p-value")
    tests$`p-value` <- format.pval(tests$`p-value`)
    cat("Test equality across groups:\n")
    print(tests, digits = digits) 
  }
  cat("\nEstimates at time points:\n")
  print(tfit$est, digits = digits)
  cat("\nStandard errors:\n")
  print(sqrt(tfit$var), digits = digits)
  
  if(missing(level))
  { 
    if(missing(t))
    { time <- sort(unique(c(ftime, time)))
    x <- timepoints(fit, time) }
    else x <- tfit
    col <- if(missing(col)) rep(1:(s-1), rep(g,(s-1))) else col
    lty <- if(missing(lty)) rep(1:g, s-1) else lty
    lwd <- if(missing(lwd)) rep(1, g*(s-1)) else lwd      
    matplot(time, base::t(x$est), type="s", ylim = c(0,1), 
            xlab = xlab, ylab = ylab, xaxs="i", yaxs="i", 
            col = col, lty = lty, lwd = lwd)
    legend("topleft", legend =  rownames(x$est), x.intersp = 2, 
           bty = "n", xjust = 1, col = col, lty = lty, lwd = lwd)
    out <- list(test = tests, est = tfit$est, se = sqrt(tfit$var))
  }
  else
  { if(level < 0 | level > 1) 
    error("level must be a value in the range [0,1]")
    
    oldpar <- par(ask=TRUE)
    on.exit(par(oldpar))
    if(missing(t))
    { time <- sort(unique(c(ftime, time)))
    x <- timepoints(fit, time) }
    else x <- tfit
    z <- qnorm(1-(1-level)/2)
    lower <- x$est ^ exp(-z*sqrt(x$var)/(x$est*log(x$est)))
    upper <- x$est ^ exp(z*sqrt(x$var)/(x$est*log(x$est)))
    col <- if(missing(col)) rep(1:(s-1), rep(g,(s-1))) 
    else             rep(col, g*(s-1))
    lwd <- if(missing(lwd)) rep(1, g*(s-1)) 
    else             rep(lwd, g*(s-1))      
    
    for(j in 1:nrow(x$est))
    { matplot(time, cbind(x$est[j,], lower[j,], upper[j,]), type="s", 
              xlab = xlab, ylab = ylab, xaxs="i", yaxs="i", 
              ylim = c(0,1), col = col[j], lwd = lwd[j], lty = c(1,3,3))
      legend("topleft", legend =  rownames(x$est)[j], bty = "n", xjust = 1) }
    
    i <- match(ttime, time)
    ci <- array(NA, c(2, length(i), nrow(lower)))
    ci[1,,] <- base::t(lower[,i])
    ci[2,,] <- base::t(upper[,i])
    dimnames(ci) <- list(c("lower", "upper"), ttime, rownames(lower))
    cat(paste("\n", level*100, "% pointwise confidence intervals:\n\n", sep=""))
    print(ci, digits = digits)
    out <- list(test = tests, est = x$est, se = sqrt(tfit$var), ci = ci)
  }
  
  invisible(out)
}



# Falls -----------------
falls=read.csv("falls.csv", sep=";")
str(falls)
attach(falls)
disfalls=factor(disfalls, levels=c(0), labels= c("PD"))
table(disfalls, statusfalls)
fitfalls=CumIncidence (ftimefalls, statusfalls, cencode = 0, xlab="Months", t=c(0, 3, 6, 9, 12, 24, 36, 48,72,96,97), level = 0.95)



# Freezing of Gait -------------
freezing=read.csv("freezing.csv", sep=";")
str(freezing)
attach(freezing)
disfreezing=factor(disfreezing, levels=c(0), labels= c("PD"))
table(disfreezing, statusfreezing)
fitfreezing=CumIncidence (ftimefreezing, statusfreezing, cencode = 0, xlab="Months", t=c(0, 3, 6, 9, 12, 24, 36, 48,72,96,97), level = 0.95)



# Hallucinations ----------
hallucinations=read.csv("hallucinations.csv", sep=";")
str(hallucinations)
attach(hallucinations)
disfreezing=factor(dishallucinations, levels=c(0), labels= c("PD"))
table(dishallucinations, statushallucinations)
fithallucinations=CumIncidence (ftimehallucinations, statushallucinations, cencode = 0, xlab="Months", t=c(0, 3, 6, 9, 12, 24, 36, 48,72,96,97), level = 0.95)


# Dementia ----------------
dementia=read.csv("dementia.csv", sep=";")
str(dementia)
attach(dementia)
disdementia=factor(disdementia, levels=c(0), labels= c("PD"))
table(disdementia, statusdementia)
fitdementia=CumIncidence (ftimedementia, statusdementia, cencode = 0, xlab="Months", t=c(0, 3, 6, 9, 12, 24, 36, 48,72,96,97), level = 0.95)


# Institutionalization ----------------
nursinghome=read.csv("nursinghome.csv", sep=";")
str(nursinghome)
attach(nursinghome)
nursinghome=factor(disnursinghome, levels=c(0), labels= c("PD"))
table(disnursinghome, statusnursinghome)
fitnursinghome=CumIncidence (ftimenursinghome, statusnursinghome, cencode = 0, xlab="Months", t=c(0, 3, 6, 9, 12, 24, 36, 48,72,96,97), level = 0.95)

# Subdistribution Hazards

if(!require(cmprsk))
{ stop("the package 'cmprsk' is required, please install it. \nSee help(install.packages).") }

factor2ind <- function(x, baseline)
{
  #### dummy variable encoding ####
  xname <- deparse(substitute(x))
  n <- length(x)
  x <- as.factor(x)
  if(!missing(baseline)) x <- relevel(x, baseline)
  X <- matrix(0, n, length(levels(x)))
  X[(1:n) + n*(unclass(x)-1)] <- 1
  X[is.na(x),] <- NA
  dimnames(X) <- list(names(x), paste(xname, levels(x), sep = ":"))
  return(X[,-1,drop=FALSE])
}

modsel.crr <- function (object, ..., d = log(object$n)) 
{
  if(class(object) != "crr") 
    stop("object is not of class 'crr'")
  objects <- list(object, ...)
  nmodels <- length(objects)
  modnames <- paste("Model ", format(1:nmodels), ": ", 
                    lapply(objects, function(x) x$call), 
                    sep = "", collapse = "\n")
  
  mod0 <- object
  mod0$loglik <- mod0$loglik.null
  mod0$coef <- mod0$call$cov1 <- mod0$call$cov2 <- NULL
  objects <- c(list(mod0), objects)
  nmodels <- nmodels + 1
  
  modnames <- c("Model 0: Null model", modnames)
  ns <- sapply(objects, function(x) x$n) 
  dfs <- sapply(objects, function(x) length(x$coef)) 
  if(any(ns != ns[1]))
    stop("models were not all fitted to the same dataset")
  out <- matrix(rep(NA, 5 * nmodels), ncol = 5)
  loglik <- sapply(objects, function(x) x$loglik)
  crit <- sapply(objects, function(x) -2*x$loglik + d*length(x$coef))
  out[,1] <- ns
  out[,2] <- loglik
  out[,3] <- dfs
  out[,4] <- crit
  out[,5] <- crit - min(crit)
  if(d==log(object$n)) critname <- "BIC"
  else if(d == 2) critname <- "AIC"
  else critname <- "Criterion"
  colnames(out) <- c("Num.obs", "logLik", "Df.fit", critname, paste(critname, "diff"))
  rownames(out) <- 0:(nmodels-1)
  title <- "Model selection table\n"
  topnote <- modnames
  structure(as.data.frame(out), heading = c(title, topnote), 
            class = c("anova", "data.frame"))
}








# Measuring the contribution of each variable onto the incidence of each milestone event
# in the presence of a competing risk -> Death (event 0 on the column " *variable*status ")


# Falls Regression -----------------
falls_predictions = read.csv("falls_predictions.csv", sep = ";", na.strings=c("", " ","NA"))
falls_predictions$Sex <- as.factor(falls_predictions$Sex)
falls_predictions$PostInstOFFmore <- as.factor(falls_predictions$PostInstOFFmore)
falls_predictions$GAITOFFmore <- as.factor(falls_predictions$GAITOFFmore)
falls_predictions$LCTrespMORE <- as.factor(falls_predictions$LCTrespMORE)
newfalls = cbind(falls_predictions$AgeSurgery, 
                 falls_predictions$DiseaseDurationSurgery, 
                 factor2ind(falls_predictions$Sex, "F"), 
                 factor2ind(falls_predictions$PostInstOFFmore, "ZeroOne"),
                 falls_predictions$UPDRS3off, 
                 falls_predictions$UPDRS2off,
                 factor2ind(falls_predictions$LCTrespMORE, "Less50"))
mod1 = crr(falls_predictions$ftimefalls, falls_predictions$statusfalls, newfalls)
summary(mod1)




# Freezing Regression ------------------
freezing_predictions = read.csv("freezing_predictions.csv", sep = ";", na.strings=c("", " ","NA"))
freezing_predictions$Sex <- as.factor(freezing_predictions$Sex)
freezing_predictions$PostInstOFFmore <- as.factor(freezing_predictions$PostInstOFFmore)
freezing_predictions$GAITOFFmore <- as.factor(freezing_predictions$GAITOFFmore)
freezing_predictions$LCTrespMORE <- as.factor(freezing_predictions$LCTrespMORE)
newfreezing = cbind(freezing_predictions$AgeSurgery, 
                    freezing_predictions$DiseaseDurationSurgery, 
                    factor2ind(freezing_predictions$Sex, "F"), 
                    factor2ind(freezing_predictions$PostInstOFFmore, "ZeroOne"),
                    factor2ind(freezing_predictions$GAITOFFmore, "ZeroOROne"),
                    freezing_predictions$UPDRS3off, 
                    freezing_predictions$UPDRS2off,
                    factor2ind(freezing_predictions$LCTrespMORE, "Less50"))
modfreezing = crr(freezing_predictions$ftimefreezing, freezing_predictions$statusfreezing, newfreezing)
summary(modfreezing)



# Hallucinations Regression ---------------
hallucinations_predictions = read.csv("hallucinations_predictions.csv", sep = ";", na.strings=c("", " ","NA"))
hallucinations_predictions$Sex <- as.factor(hallucinations_predictions$Sex)
hallucinations_predictions$PostInstOFFmore <- as.factor(hallucinations_predictions$PostInstOFFmore)
hallucinations_predictions$GAITOFFmore <- as.factor(hallucinations_predictions$GAITOFFmore)
hallucinations_predictions$LCTrespMORE <- as.factor(hallucinations_predictions$LCTrespMORE)
hallucinations_predictions$Neuropsychologic <- as.factor(hallucinations_predictions$Neuropsychologic)
newhallucinations = cbind(hallucinations_predictions$AgeSurgery, 
                          hallucinations_predictions$DiseaseDurationSurgery, 
                          factor2ind(hallucinations_predictions$Sex, "F"), 
                          hallucinations_predictions$UPDRS3off, 
                          hallucinations_predictions$UPDRS2off,
                          hallucinations_predictions$UPDRS1preop,
                          factor2ind(hallucinations_predictions$LCTrespMORE, "Less50"),
                          factor2ind(hallucinations_predictions$Neuropsychologic, "Normal"),
                          hallucinations_predictions$LEDDpreop)
modnewhallucinations = crr(hallucinations_predictions$ftimehallucinations, hallucinations_predictions$statushallucinations, newhallucinations)
summary(modnewhallucinations)







# Dementia Regression ------------------
dementia_predictions = read.csv("dementia_predictions.csv", sep = ";", na.strings=c("", " ","NA"))
dementia_predictions$Sex <- as.factor(dementia_predictions$Sex)
dementia_predictions$PostInstOFFmore <- as.factor(dementia_predictions$PostInstOFFmore)
dementia_predictions$GAITOFFmore <- as.factor(dementia_predictions$GAITOFFmore)
dementia_predictions$LCTrespMORE <- as.factor(dementia_predictions$LCTrespMORE)
dementia_predictions$Neuropsychologic <- as.factor(dementia_predictions$Neuropsychologic)
newdementia = cbind(dementia_predictions$AgeSurgery, 
                    dementia_predictions$DiseaseDurationSurgery, 
                    factor2ind(dementia_predictions$Sex, "F"), 
                    dementia_predictions$UPDRS3off, 
                    dementia_predictions$UPDRS2off,
                    dementia_predictions$UPDRS1preop,
                    factor2ind(dementia_predictions$LCTrespMORE, "Less50"),
                    factor2ind(dementia_predictions$Neuropsychologic, "Normal"))
modnewdementia = crr(dementia_predictions$ftimedementia, dementia_predictions$statusdementia, newdementia)
summary(modnewdementia)




# Institutionalization Regression ------------------
nursing_predictions = read.csv("nursing_predictions.csv", sep = ";", na.strings=c("", " ","NA"))
nursing_predictions$Sex <- as.factor(nursing_predictions$Sex)
nursing_predictions$PostInstOFFmore <- as.factor(nursing_predictions$PostInstOFFmore)
nursing_predictions$GAITOFFmore <- as.factor(nursing_predictions$GAITOFFmore)
nursing_predictions$LCTrespMORE <- as.factor(nursing_predictions$LCTrespMORE)
nursing_predictions$Neuropsychologic <- as.factor(nursing_predictions$Neuropsychologic)
newnursing = cbind(nursing_predictions$AgeSurgery, 
                   nursing_predictions$DiseaseDurationSurgery, 
                   factor2ind(nursing_predictions$Sex, "F"), 
                   nursing_predictions$UPDRS3off, 
                   nursing_predictions$UPDRS2off,
                   nursing_predictions$UPDRS1preop,
                   factor2ind(nursing_predictions$LCTrespMORE, "Less50"),
                   factor2ind(nursing_predictions$Neuropsychologic, "Normal"))

modnewnursing= crr(nursing_predictions$ftimenursing, nursing_predictions$statusnursing, newnursing)
summary(modnewnursing)





# BACKWARDS STEEPING MODEL , FINAL PREICTORS (P <= 0.2)

# Falls Regression  Backwards Stepping -----------------
falls_predictions = read.csv("falls_predictions.csv", sep = ";", na.strings=c("", " ","NA"))
falls_predictions$Sex <- as.factor(falls_predictions$Sex)
falls_predictions$PostInstOFFmore <- as.factor(falls_predictions$PostInstOFFmore)
falls_predictions$GAITOFFmore <- as.factor(falls_predictions$GAITOFFmore)
falls_predictions$LCTrespMORE <- as.factor(falls_predictions$LCTrespMORE)

newfallsback = cbind(falls_predictions$AgeSurgery, 
                                          falls_predictions$UPDRS2off,
                                         falls_predictions$LEDDpreop)

modnewfallsback = crr(falls_predictions$ftimefalls, falls_predictions$statusfalls, newfallsback)

summary(modnewfallsback)





# Freezing Regression  Backwards Stepping ------------------
freezing_predictions = read.csv("freezing_predictions.csv", sep = ";", na.strings=c("", " ","NA"))
freezing_predictions$Sex <- as.factor(freezing_predictions$Sex)
freezing_predictions$PostInstOFFmore <- as.factor(freezing_predictions$PostInstOFFmore)
freezing_predictions$GAITOFFmore <- as.factor(freezing_predictions$GAITOFFmore)
freezing_predictions$LCTrespMORE <- as.factor(freezing_predictions$LCTrespMORE)

newfreezingback = cbind(freezing_predictions$AgeSurgery, 
                                            factor2ind(freezing_predictions$GAITOFFmore, "ZeroOROne"),
                                            freezing_predictions$UPDRS3on,
                                            freezing_predictions$UPDRS2off,
                                           freezing_predictions$UPDRS2on)


modnewfreezingback = crr(freezing_predictions$ftimefreezing, freezing_predictions$statusfreezing, newfreezingback)


summary(modnewfreezingback)




                 
# Hallucinations Regression Backwards Stepping ---------------
hallucinations_predictions = read.csv("hallucinations_predictions.csv", sep = ";", na.strings=c("", " ","NA"))
hallucinations_predictions$Sex <- as.factor(hallucinations_predictions$Sex)
hallucinations_predictions$PostInstOFFmore <- as.factor(hallucinations_predictions$PostInstOFFmore)
hallucinations_predictions$GAITOFFmore <- as.factor(hallucinations_predictions$GAITOFFmore)
hallucinations_predictions$LCTrespMORE <- as.factor(hallucinations_predictions$LCTrespMORE)
hallucinations_predictions$Neuropsychologic <- as.factor(hallucinations_predictions$Neuropsychologic)

newhallucinationsback = cbind( factor2ind(hallucinations_predictions$Phenotype, "Tremor"))

modnewhallucinationsback = crr(hallucinations_predictions$ftimehallucinations, hallucinations_predictions$statushallucinations, newhallucinationsback)

summary(modnewhallucinationsback)


                 


# Dementia Regression Backwards Stepping ------------------
dementia_predictions = read.csv("dementia_predictions.csv", sep = ";", na.strings=c("", " ","NA"))
dementia_predictions$Sex <- as.factor(dementia_predictions$Sex)
dementia_predictions$PostInstOFFmore <- as.factor(dementia_predictions$PostInstOFFmore)
dementia_predictions$GAITOFFmore <- as.factor(dementia_predictions$GAITOFFmore)
dementia_predictions$LCTrespMORE <- as.factor(dementia_predictions$LCTrespMORE)
dementia_predictions$Neuropsychologic <- as.factor(dementia_predictions$Neuropsychologic)


newdementiaback = cbind(dementia_predictions$AgeSurgery, 
                                             dementia_predictions$DiseaseDurationSurgery, 
                                             dementia_predictions$UPDRS2off,
                                             dementia_predictions$UPDRS2on,
                                             factor2ind(dementia_predictions$PostInstOFFmore, "ZeroOne"))
modnewdementiaback = crr(dementia_predictions$ftimedementia, dementia_predictions$statusdementia, newdementiaback)

summary(modnewdementiaback)



                 

# Hospitalization Regression Backwards Stepping ---------
nursing_predictions = read.csv("nursing_predictions.csv", sep = ";", na.strings=c("", " ","NA"))
nursing_predictions$Sex <- as.factor(nursing_predictions$Sex)
nursing_predictions$PostInstOFFmore <- as.factor(nursing_predictions$PostInstOFFmore)
nursing_predictions$GAITOFFmore <- as.factor(nursing_predictions$GAITOFFmore)
nursing_predictions$LCTrespMORE <- as.factor(nursing_predictions$LCTrespMORE)
nursing_predictions$Neuropsychologic <- as.factor(nursing_predictions$Neuropsychologic)

newnursingback = cbind(nursing_predictions$AgeSurgery, 
                                           nursing_predictions$UPDRS2off)
 

modnewnursingback = crr(nursing_predictions$ftimenursing, nursing_predictions$statusnursing, newnursingback)

summary(modnewnursingback)

                 







#To compare the performance of different models (i.e. backwards stepping models for each event--------
install.packages("aod")
library(aod)

## checking whether the explanatory variables in a given model are significant or not 
wald.test(mod17$var, mod17$coef, Terms = 1:2)

# Creating different models
# ftime: time to event, 
# status: event happened or competing risk (death) occurred
# [ x,   "all_rows" : "which_columns" ]

mod1 = crr(filename$ftime, filename$Status, x[,c(4,7)]) # vars 4,7
mod2 = crr(filename$ftime, filename$Status, x[,c(4,2,3)]) # vars 4,2,3
mod3 = crr(filename$ftime, filename$Status, x[,c(4,5,2)]) # vars 4,5,2

modsel.crr(mod1, mod2, mod3)

# The aforementioned model comparisons will only work in the absence of missing values 
# or if NAs occur in a systematic way
# "n" needs to be the same for all models 
# if different patients have different NAs on different variables, modesel.crr won't work









# Correlation Matrix -------------
Correlation_Matrix_Input_Data <- read.csv("Correlation_Matrix_Input_Data.csv", sep = ";", header = T)

cor_data = rcorr(as.matrix(Correlation_Matrix_Input_Data),  type=c("spearman"))

corr_values <- as.data.frame(cor_data$r)
p_values <- as.data.frame(cor_data$P)

ggcorrplot(corr_values, method=c("circle"),  colors = c("yellow", "white", "darkslateblue"), outline.color = "white", hc.order = TRUE, type = "lower", lab = TRUE)
ggcorrplot(p_values, method=c("circle"),  colors = c("white", "white", "cadetblue2"), outline.color = "white", hc.order = TRUE, type = "lower", lab = TRUE, show.legend = F)













# Logistic Regression Falls ---------------------

# Falls
Multiple_Logistic_Regression_Falls <- read.csv("Multiple_Logistic_Regression_Falls.csv", sep = ";", header = T)
Multiple_Logistic_Regression_Falls <-Multiple_Logistic_Regression_Falls %>% drop_na()

Falls_fit_allVars <- glm(Falls ~ ageatsurgery + diseaseduration + LEDD + UPDRS2ON + UPDRS2OFF + UPDRS3OFF + UPDRS3ON + PostInstOFF + GaitOFF+ LD.response + MMSE, data = Multiple_Logistic_Regression_Falls, family = binomial, na.action = na.pass)

summary(Falls_fit_allVars)


roc(Multiple_Logistic_Regression_Falls$Falls, as.vector(fitted.values(Falls_fit_allVars)), percent=F,   boot.n=1000, ci.alpha=0.9, stratified=FALSE, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, ci=TRUE,  col="deepskyblue4", main = paste("ROC curve for Falls \n ~All Variables \n") )

Falls_fit_UPDRS2OFF_PostInstVars <- glm(Falls ~ UPDRS2OFF + PostInstOFF, data = Multiple_Logistic_Regression_Falls, family = binomial, na.action = na.pass)

summary(Falls_fit_UPDRS2OFF_PostInstVars)



roc(Multiple_Logistic_Regression_Falls$Falls, as.vector(fitted.values(Falls_fit_UPDRS2OFF_PostInstVars)), percent=F,   boot.n=1000, ci.alpha=0.9, stratified=FALSE, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, ci=TRUE,  col="deeppink4", main = paste("ROC curve for Falls \n ~UPDRS 2 OFF + Postural Instability \n") )




# Partition data: 80% vs 20% split for FALLS -------------

training.samples.falls <- Multiple_Logistic_Regression_Falls$Falls %>% 
  createDataPartition(p = 0.8, list = FALSE)

train.data.falls  <- Multiple_Logistic_Regression_Falls[training.samples.falls, ]
test.data.falls <- Multiple_Logistic_Regression_Falls[-training.samples.falls, ]

model <- glm( Falls ~., data = train.data.falls, family = binomial)

summary(model)


probabilities <- model %>% predict(test.data.falls, type = "response")

predicted.classes <- ifelse(probabilities > 0.5, 1, 0)

mean(predicted.classes == test.data.falls$Falls) # 0.5333333
 

model <- glm(Falls ~ ageatsurgery + diseaseduration + LEDD + UPDRS2ON + UPDRS2OFF + UPDRS3OFF + UPDRS3ON + PostInstOFF + GaitOFF+ LD.response + MMSE, data = train.data.falls, family = binomial, na.action = na.pass)

summary(model)$coef

probabilities <- model %>% predict(test.data.falls, type = "response")

contrasts(as.factor(test.data.falls$Falls))

predicted.classes <- ifelse(probabilities > 0.5, 1, 0)

mean(predicted.classes == test.data.falls$Falls) # 0.5333333








# Logistic Regression Freezing ------------------------------

# Freezing
Multiple_Logistic_Regression_Freezing <- read.csv("Multiple_Logistic_Regression_Freezing.csv", sep = ";", header = T)
Multiple_Logistic_Regression_Freezing <- Multiple_Logistic_Regression_Freezing %>% drop_na()

Freezing_fit_allVars <- glm(Freezing ~ ageatsurgery + diseaseduration + LEDD + UPDRS2ON + UPDRS2OFF + UPDRS3OFF + UPDRS3ON + PostInstOFF + GaitOFF+ LD.response + MMSE, data = Multiple_Logistic_Regression_Freezing, family = binomial, na.action = na.pass)

summary(Freezing_fit_allVars)

 
# roc(Multiple_Logistic_Regression_Freezing$Freezing, as.vector(fitted.values(Freezing_fit_allVars)))
# Area under the curve: 0.7435

roc(Multiple_Logistic_Regression_Freezing$Freezing, as.vector(fitted.values(Freezing_fit_allVars)), percent=F,   boot.n=1000, ci.alpha=0.9, stratified=FALSE, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, ci=TRUE,  col="deepskyblue4", main = paste("ROC curve for Freezing \n ~All Variables \n") )

Freezing_fit_GAITOFF <- glm(Freezing ~ GaitOFF, data = Multiple_Logistic_Regression_Freezing, family = binomial, na.action = na.pass)

summary(Freezing_fit_GAITOFF)


roc(Multiple_Logistic_Regression_Freezing$Freezing, as.vector(fitted.values(Freezing_fit_GAITOFF)), percent=F,   boot.n=1000, ci.alpha=0.9, stratified=FALSE, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, ci=TRUE,  col="deeppink4", main = paste("ROC curve for Freezing \n ~GAIT OFF \n") )




# Partition data: 80% vs 20% split for FREEZING -------------

training.samples.freezing <- Multiple_Logistic_Regression_Freezing$Freezing %>% 
  createDataPartition(p = 0.8, list = FALSE)

train.data.freezing  <- Multiple_Logistic_Regression_Freezing[training.samples.freezing, ]
test.data.freezing <- Multiple_Logistic_Regression_Freezing[-training.samples.freezing, ]

model <- glm( Freezing ~., data = train.data.freezing, family = binomial)

summary(model)


probabilities <- model %>% predict(test.data.freezing, type = "response")

predicted.classes <- ifelse(probabilities > 0.5, 1, 0)

mean(predicted.classes == test.data.freezing$Freezing) # 0.6666667


model <- glm(Freezing ~ ageatsurgery + diseaseduration + LEDD + UPDRS2ON + UPDRS2OFF + UPDRS3OFF + UPDRS3ON + PostInstOFF + GaitOFF+ LD.response + MMSE, data = train.data.freezing, family = binomial, na.action = na.pass)

summary(model)$coef



probabilities <- model %>% predict(test.data.freezing, type = "response")

contrasts(as.factor(test.data.freezing$Freezing))

predicted.classes <- ifelse(probabilities > 0.5, 1, 0)

mean(predicted.classes == test.data.freezing$Freezing) # 0.6666667





# Plotting Figures -----------------------------------------------

# Cumulative Incidence for each Event --------------------------------
Milestones_Summary <- read.csv("Milestones_Summary.csv", sep = ";", header = T)
event_order <- c("Falls", "Freezing", "Dementia", "Hallucinations", "Institutionalization")

Milestones_Summary <- Milestones_Summary %>% mutate(Event = factor(Event, levels = event_order))

Milestones_Summary %>% 
  ggplot(mapping = aes(x = Follow_up_months, y = Proportion)) +
  geom_step(aes(color = Event), show.legend = FALSE) +
  geom_stepconfint(aes(ymin = conf.low, ymax = conf.high, fill = Event), alpha = 0.6) +
  ylim(0,1)+
  labs(x = "\n Follow-up (Months)", y = "Cumulative Incidence \n") +
  scale_fill_viridis_d()+
  scale_colour_viridis_d()+
  theme_minimal()+
  facet_wrap(~Event, scales = "free")


# Summary With All Events Cumulative Incidence  ----------------------------
Milestones_Each <- read.csv("Milestones_Each.csv", sep = ";", header = T)
Milestones_Each <- Milestones_Each %>% mutate(Event = factor(Event, levels = event_order))

Milestones_Each %>% 
  mutate(Follow_up_months = log10(Follow_up_months)) %>%
  ggplot(mapping = aes(x = Follow_up_months, y = Proportion)) +
  geom_line(aes(color = Event), show.legend = FALSE, size=2) +
  ylim(0,1)+
  labs(x = NULL, y = NULL) +
  scale_fill_viridis_d()+
  scale_colour_viridis_d()+
  theme_minimal()


# Events Incidence Paired Groups (with vs without accompanying event) ---------------
Milestones_Grouped <- read.csv("Milestones_Grouped.csv", sep = ";", header = T)

Milestones_Grouped[,3] <- as.factor(Milestones_Grouped[,3])
Milestones_Grouped[,4] <- as.factor(Milestones_Grouped[,4])

Milestones_Grouped %>% 
  ggplot(mapping = aes(x = Follow_up_months, y = Proportion)) +
  geom_step(aes(color = Group), show.legend = FALSE, size=1.5) +
  geom_stepconfint(aes(ymin = conf.low, ymax = conf.high, fill = Group), alpha = 0.4) +
  ylim(0,1)+
  labs(x = "\n Follow-up (Months)", y = "Cumulative Incidence \n") +
  scale_fill_viridis_d()+
  scale_colour_viridis_d()+
  theme_minimal()+
  facet_wrap(~Event, scales = "free")


# Survival ~ Milestone event development  ---------------
Milestones_Mortality_per_Group <- read.csv("Milestones_Mortality_per_Group.csv", sep = ";", header = T)

Milestones_Mortality_per_Group[,3] <- as.factor(Milestones_Mortality_per_Group[,3])
Milestones_Mortality_per_Group[,4] <- as.factor(Milestones_Mortality_per_Group[,4])

event_order_2 <- c("Overall_Survival", "Survival_Freezing", "Survival_Falls", "Survival_Hallucinations", "Survival_Dementia" , "Survival_Institutionalization")

Milestones_Mortality_per_Group <- Milestones_Mortality_per_Group %>% mutate(Event = factor(Event, levels = event_order_2))

Milestones_Mortality_per_Group %>% 
  group_by(Event, Group) %>%
  fill(Proportion) %>%
  fill(conf.high) %>% fill(conf.high, .direction = c("up")) %>%
  fill(conf.low) %>% fill(conf.low, .direction = c("up")) %>%
  ungroup() %>%
  ggplot(mapping = aes(x = Follow_up_months, y = Proportion)) +
  geom_step(aes(color = Group), show.legend = FALSE, size=2) +
  geom_stepconfint(aes(ymin = conf.low, ymax = conf.high, fill = Group), alpha = 0.2) +
  ylim(0,100)+xlim(0,110)+
  labs(x = "\n Follow-up (Months)", y = "Survival (%) \n") +
  scale_fill_viridis_d(option = "D")+
  scale_colour_viridis_d(option = "D")+
  theme_minimal()+
  facet_wrap(~Event, scales = "free")




# Check for drug usage  -------------

MilestonesDrugUsage <- read.csv("MilestonesDrugUsage.csv", sep = ";", header = T)

for(x in 4:17) { print(sum(MilestonesDrugUsage[,x] != 0)) }



# Stimulation settings per group  -------------

Stimul_Settings <- read.csv("Stimul_Settings.csv", sep = ";", header = T)

Stimul_Settings <- Stimul_Settings[,3:10]

names(Stimul_Settings)

# Right Side
Stimul_Settings %>% group_by(modeSTNrightEND) %>% count()
Stimul_Settings %>% group_by(modeSTNrightEND) %>% summarise(n = mean(voltageSTNrightEND, na.rm=T))
Stimul_Settings %>% group_by(modeSTNrightEND) %>% summarise(n = sd(voltageSTNrightEND, na.rm=T))
Stimul_Settings %>% group_by(modeSTNrightEND) %>% summarise(n = mean(frequencySTNrightEND, na.rm=T))
Stimul_Settings %>% group_by(modeSTNrightEND) %>% summarise(n = sd(frequencySTNrightEND, na.rm=T))
Stimul_Settings %>% group_by(modeSTNrightEND) %>% summarise(n = mean(pulsewidthSTNrightEND, na.rm=T))
Stimul_Settings %>% group_by(modeSTNrightEND) %>% summarise(n = sd(pulsewidthSTNrightEND, na.rm=T))

# Left Side
Stimul_Settings %>% group_by(modeSTNLeftEND) %>% count()
Stimul_Settings %>% group_by(modeSTNLeftEND) %>% summarise(n = mean(voltageSTNLeftEND, na.rm=T))
Stimul_Settings %>% group_by(modeSTNLeftEND) %>% summarise(n = sd(voltageSTNLeftEND, na.rm=T))
Stimul_Settings %>% group_by(modeSTNLeftEND) %>% summarise(n = mean(frequencySTNLeftEND, na.rm=T))
Stimul_Settings %>% group_by(modeSTNLeftEND) %>% summarise(n = sd(frequencySTNLeftEND, na.rm=T))
Stimul_Settings %>% group_by(modeSTNLeftEND) %>% summarise(n = mean(pulsewidthSTNLeftEND, na.rm=T))
Stimul_Settings %>% group_by(modeSTNLeftEND) %>% summarise(n = sd(pulsewidthSTNLeftEND, na.rm=T))

