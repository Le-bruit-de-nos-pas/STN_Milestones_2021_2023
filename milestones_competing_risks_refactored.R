############################################################
# Milestones & Competing Risks Analyses
# Parkinson's disease DBS cohort
############################################################

## Packages -------------------------------------------------

required_pkgs <- c(
  "survival", "survminer", "lubridate", "cmprsk", "tidyverse",
  "utile.visuals", "Hmisc", "corrplot", "ggcorrplot",
  "pROC", "caret", "mgcv"
)

invisible(lapply(required_pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    stop(sprintf("Package '%s' is required. Please install it.", pkg))
  }
}))

## Helper functions -----------------------------------------

# Cumulative incidence (Scrucca et al. 2007)
CumIncidence <- function(ftime, fstatus, group, t, strata, rho = 0,
                         cencode = 0, subset, na.action = na.omit, level,
                         xlab = "Time", ylab = "Probability",
                         col, lty, lwd, digits = 4) {

  if (!require("cmprsk")) {
    stop(
      "Package 'cmprsk' is required and must be installed.\n",
      "See help(install.packages) and then run:\n",
      "  install.packages(\"cmprsk\")"
    )
  }

  mf <- match.call(expand.dots = FALSE)
  mf[[1]] <- as.name("list")
  mf$t <- mf$digits <- mf$col <- mf$lty <- mf$lwd <- mf$level <-
    mf$xlab <- mf$ylab <- NULL
  mf <- eval(mf, parent.frame())

  g <- max(1, length(unique(mf$group)))
  s <- length(unique(mf$fstatus))

  if (missing(t)) {
    time <- pretty(c(0, max(mf$ftime)), 6)
    ttime <- time <- time[time < max(mf$ftime)]
  } else {
    ttime <- time <- t
  }

  fit <- do.call("cuminc", mf)
  tfit <- timepoints(fit, time)

  cat("\n+", paste(rep("-", 67), collapse = ""), "+", sep = "")
  cat("\n| Cumulative incidence function estimates from competing risks data |")
  cat("\n+", paste(rep("-", 67), collapse = ""), "+\n", sep = "")

  tests <- NULL
  if (g > 1) {
    tests <- data.frame(fit$Tests[, c(1, 3, 2)], check.names = FALSE)
    colnames(tests) <- c("Statistic", "df", "p-value")
    tests$`p-value` <- format.pval(tests$`p-value`)
    cat("Test equality across groups:\n")
    print(tests, digits = digits)
  }

  cat("\nEstimates at time points:\n")
  print(tfit$est, digits = digits)

  cat("\nStandard errors:\n")
  print(sqrt(tfit$var), digits = digits)

  if (missing(level)) {
    if (missing(t)) {
      time <- sort(unique(c(ftime, time)))
      x <- timepoints(fit, time)
    } else {
      x <- tfit
    }

    col <- if (missing(col)) rep(seq_len(s - 1), each = g) else col
    lty <- if (missing(lty)) rep(seq_len(g), times = s - 1) else lty
    lwd <- if (missing(lwd)) rep(1, g * (s - 1)) else lwd

    matplot(
      time, base::t(x$est), type = "s", ylim = c(0, 1),
      xlab = xlab, ylab = ylab, xaxs = "i", yaxs = "i",
      col = col, lty = lty, lwd = lwd
    )
    legend(
      "topleft", legend = rownames(x$est), x.intersp = 2,
      bty = "n", xjust = 1, col = col, lty = lty, lwd = lwd
    )

    out <- list(test = tests, est = tfit$est, se = sqrt(tfit$var))
  } else {
    if (level < 0 || level > 1)
      stop("level must be a value in the range [0, 1]")

    oldpar <- par(ask = TRUE)
    on.exit(par(oldpar))

    if (missing(t)) {
      time <- sort(unique(c(ftime, time)))
      x <- timepoints(fit, time)
    } else {
      x <- tfit
    }

    z <- qnorm(1 - (1 - level) / 2)
    lower <- x$est ^ exp(-z * sqrt(x$var) / (x$est * log(x$est)))
    upper <- x$est ^ exp( z * sqrt(x$var) / (x$est * log(x$est)))

    col <- if (missing(col)) rep(seq_len(s - 1), each = g) else rep(col, g * (s - 1))
    lwd <- if (missing(lwd)) rep(1, g * (s - 1)) else rep(lwd, g * (s - 1))

    for (j in seq_len(nrow(x$est))) {
      matplot(
        time, cbind(x$est[j, ], lower[j, ], upper[j, ]),
        type = "s",
        xlab = xlab, ylab = ylab, xaxs = "i", yaxs = "i",
        ylim = c(0, 1), col = col[j], lwd = lwd[j], lty = c(1, 3, 3)
      )
      legend(
        "topleft", legend = rownames(x$est)[j],
        bty = "n", xjust = 1
      )
    }

    i <- match(ttime, time)
    ci <- array(NA_real_, c(2, length(i), nrow(lower)))
    ci[1, , ] <- base::t(lower[, i])
    ci[2, , ] <- base::t(upper[, i])
    dimnames(ci) <- list(c("lower", "upper"), ttime, rownames(lower))

    cat(paste0("\n", level * 100, "% pointwise confidence intervals:\n\n"))
    print(ci, digits = digits)

    out <- list(
      test = tests,
      est  = x$est,
      se   = sqrt(tfit$var),
      ci   = ci
    )
  }

  invisible(out)
}

# Dummy encoding for factors
factor2ind <- function(x, baseline) {
  xname <- deparse(substitute(x))
  n <- length(x)
  x <- as.factor(x)
  if (!missing(baseline)) x <- relevel(x, baseline)

  X <- matrix(0, n, length(levels(x)))
  X[(seq_len(n)) + n * (unclass(x) - 1)] <- 1
  X[is.na(x), ] <- NA
  dimnames(X) <- list(names(x), paste(xname, levels(x), sep = ":"))
  X[, -1, drop = FALSE]
}

# Model selection for crr models
modsel.crr <- function(object, ..., d = log(object$n)) {
  if (!inherits(object, "crr"))
    stop("object is not of class 'crr'")

  objects  <- list(object, ...)
  nmodels  <- length(objects)
  modnames <- paste0(
    "Model ", format(seq_len(nmodels)), ": ",
    lapply(objects, function(x) x$call),
    collapse = "\n"
  )

  mod0 <- object
  mod0$loglik <- mod0$loglik.null
  mod0$coef   <- NULL
  mod0$call$cov1 <- mod0$call$cov2 <- NULL

  objects  <- c(list(mod0), objects)
  nmodels  <- nmodels + 1
  modnames <- c("Model 0: Null model", modnames)

  ns  <- sapply(objects, function(x) x$n)
  dfs <- sapply(objects, function(x) length(x$coef))

  if (any(ns != ns[1]))
    stop("models were not all fitted to the same dataset")

  out <- matrix(NA_real_, nrow = nmodels, ncol = 5)

  loglik <- sapply(objects, function(x) x$loglik)
  crit   <- sapply(objects, function(x) -2 * x$loglik + d * length(x$coef))

  out[, 1] <- ns
  out[, 2] <- loglik
  out[, 3] <- dfs
  out[, 4] <- crit
  out[, 5] <- crit - min(crit)

  critname <- if (d == log(object$n)) {
    "BIC"
  } else if (d == 2) {
    "AIC"
  } else {
    "Criterion"
  }

  colnames(out) <- c("Num.obs", "logLik", "Df.fit", critname, paste(critname, "diff"))
  rownames(out) <- 0:(nmodels - 1)

  title   <- "Model selection table\n"
  topnote <- modnames

  structure(
    as.data.frame(out),
    heading = c(title, topnote),
    class   = c("anova", "data.frame")
  )
}

## 1. Global survival ---------------------------------------

GlobalTimeToDeath <- read.csv("GlobalTimeToDeath.csv", sep = ";")

# Overall survival curve
fit_surv_global <- survfit(Surv(time, status) ~ 1, data = GlobalTimeToDeath)

ggsurvplot(
  fit = fit_surv_global,
  xlab = "Follow-up (Months)",
  ylab = "Overall survival probability",
  legend = "none",
  break.x.by = 12,
  risk.table = TRUE,
  tables.col = "strata"
)

# Survival probabilities at selected times (months)
summary(fit_surv_global, times = c(12, 24, 48, 96))

# By sex
fit_surv_sex <- survfit(Surv(time, status) ~ sex, data = GlobalTimeToDeath)

ggsurvplot(
  fit = fit_surv_sex,
  xlab = "Follow-up (Months)",
  ylab = "Overall survival probability",
  legend = "none",
  break.x.by = 12,
  risk.table = TRUE,
  tables.col = "strata",
  conf.int = TRUE
)

coxph(Surv(time, status) ~ sex, data = GlobalTimeToDeath)

## 2. Cox regression: time to death -------------------------

GlobalDeath <- read.csv(
  "GlobalDeath.csv",
  sep = ";",
  na.strings = c("", " ", "NA")
)

# Univariate models (non-adjusted)
cox_vars <- c(
  "AgeSurgery", "Sex", "DiseaseDurationSurgery", "LEDDpreop",
  "LCTresp", "AgeSurgeryMORE", "DiseaseDurationSurgeryMORE",
  "LCTrespMORE", "GAITOFFmore", "PostInstOFFmore", "HYOFFmore",
  "factor(Neuropsychologic, levels = c('Normal', 'MCI'))",
  "UPDRS3off", "UPDRS3on", "UPDRS2off", "UPDRS2on",
  "factor(Phenotype, levels = c('Tremor', 'PIGD', 'Indeterminate'))"
)

for (v in cox_vars) {
  f <- as.formula(paste("Surv(time, status) ~", v))
  m <- coxph(f, data = GlobalDeath)
  print(v)
  print(summary(m))
}

# Multivariable model building
res.coxGLOBAL  <- coxph(
  Surv(time, status) ~ AgeSurgery +
    factor(Sex, levels = c("F", "M")) +
    DiseaseDurationSurgery + LEDDpreop +
    factor(GAITOFFmore, levels = c("ZeroOROne", "EqualORMore2")) +
    factor(PostInstOFFmore, levels = c("ZeroOne", "EqualMore2")) +
    factor(Neuropsychologic, levels = c("Normal", "MCI")) +
    UPDRS3off + UPDRS3on + UPDRS2off + UPDRS2on +
    factor(Phenotype, levels = c("Tremor", "PIGD", "Indeterminate")),
  data = GlobalDeath
)
summary(res.coxGLOBAL)

res.coxGLOBAL2 <- update(
  res.coxGLOBAL,
  . ~ . - factor(Neuropsychologic, levels = c("Normal", "MCI"))
)
summary(res.coxGLOBAL2)

res.coxGLOBAL3 <- update(
  res.coxGLOBAL2,
  . ~ . - factor(GAITOFFmore, levels = c("ZeroOROne", "EqualORMore2"))
)
summary(res.coxGLOBAL3)

res.coxGLOBAL4 <- update(
  res.coxGLOBAL3,
  . ~ . - factor(Sex, levels = c("F", "M"))
)
summary(res.coxGLOBAL4)

res.coxGLOBAL5 <- update(res.coxGLOBAL4, . ~ . - DiseaseDurationSurgery)
summary(res.coxGLOBAL5)

res.coxGLOBAL6 <- update(res.coxGLOBAL5, . ~ . - LEDDpreop)
summary(res.coxGLOBAL6)

res.coxGLOBAL7 <- update(res.coxGLOBAL6, . ~ . - UPDRS2on)
summary(res.coxGLOBAL7)

res.coxGLOBAL8 <- update(res.coxGLOBAL7, . ~ . - UPDRS2off)
summary(res.coxGLOBAL8)

res.coxGLOBAL9 <- update(res.coxGLOBAL8, . ~ . - UPDRS2on)
summary(res.coxGLOBAL9)

res.coxGLOBAL10 <- update(res.coxGLOBAL9, . ~ . - UPDRS3off)
summary(res.coxGLOBAL10)

res.coxGLOBAL11 <- update(res.coxGLOBAL10, . ~ . - UPDRS3on)
summary(res.coxGLOBAL11)

## 3. Competing risks: crude cumulative incidence -----------

read_milestone_data <- function(file) {
  read.csv(file, sep = ";")
}

# Falls
falls <- read_milestone_data("falls.csv")
with(falls, {
  disfalls <- factor(disfalls, levels = 0, labels = "PD")
  print(table(disfalls, statusfalls))
  fitfalls <- CumIncidence(
    ftimefalls, statusfalls, cencode = 0,
    xlab = "Months",
    t = c(0, 3, 6, 9, 12, 24, 36, 48, 72, 96, 97),
    level = 0.95
  )
})

# Freezing
freezing <- read_milestone_data("freezing.csv")
with(freezing, {
  disfreezing <- factor(disfreezing, levels = 0, labels = "PD")
  print(table(disfreezing, statusfreezing))
  fitfreezing <- CumIncidence(
    ftimefreezing, statusfreezing, cencode = 0,
    xlab = "Months",
    t = c(0, 3, 6, 9, 12, 24, 36, 48, 72, 96, 97),
    level = 0.95
  )
})

# Hallucinations
hallucinations <- read_milestone_data("hallucinations.csv")
with(hallucinations, {
  dishallucinations <- factor(dishallucinations, levels = 0, labels = "PD")
  print(table(dishallucinations, statushallucinations))
  fithallucinations <- CumIncidence(
    ftimehallucinations, statushallucinations, cencode = 0,
    xlab = "Months",
    t = c(0, 3, 6, 9, 12, 24, 36, 48, 72, 96, 97),
    level = 0.95
  )
})

# Dementia
dementia <- read_milestone_data("dementia.csv")
with(dementia, {
  disdementia <- factor(disdementia, levels = 0, labels = "PD")
  print(table(disdementia, statusdementia))
  fitdementia <- CumIncidence(
    ftimedementia, statusdementia, cencode = 0,
    xlab = "Months",
    t = c(0, 3, 6, 9, 12, 24, 36, 48, 72, 96, 97),
    level = 0.95
  )
})

# Institutionalization
nursinghome <- read_milestone_data("nursinghome.csv")
with(nursinghome, {
  disnursinghome <- factor(disnursinghome, levels = 0, labels = "PD")
  print(table(disnursinghome, statusnursinghome))
  fitnursinghome <- CumIncidence(
    ftimenursinghome, statusnursinghome, cencode = 0,
    xlab = "Months",
    t = c(0, 3, 6, 9, 12, 24, 36, 48, 72, 96, 97),
    level = 0.95
  )
})

## 4. Competing risks regression (crr) ----------------------

read_predictions <- function(file) {
  read.csv(file, sep = ";", na.strings = c("", " ", "NA"))
}

# Falls regression
falls_predictions <- read_predictions("falls_predictions.csv")
falls_predictions <- falls_predictions %>%
  mutate(
    Sex             = factor(Sex),
    PostInstOFFmore = factor(PostInstOFFmore),
    GAITOFFmore     = factor(GAITOFFmore),
    LCTrespMORE     = factor(LCTrespMORE)
  )

newfalls <- cbind(
  falls_predictions$AgeSurgery,
  falls_predictions$DiseaseDurationSurgery,
  factor2ind(falls_predictions$Sex, "F"),
  factor2ind(falls_predictions$PostInstOFFmore, "ZeroOne"),
  falls_predictions$UPDRS3off,
  falls_predictions$UPDRS2off,
  factor2ind(falls_predictions$LCTrespMORE, "Less50")
)

mod_falls <- crr(
  falls_predictions$ftimefalls,
  falls_predictions$statusfalls,
  newfalls
)
summary(mod_falls)

# Freezing regression
freezing_predictions <- read_predictions("freezing_predictions.csv")
freezing_predictions <- freezing_predictions %>%
  mutate(
    Sex             = factor(Sex),
    PostInstOFFmore = factor(PostInstOFFmore),
    GAITOFFmore     = factor(GAITOFFmore),
    LCTrespMORE     = factor(LCTrespMORE)
  )

newfreezing <- cbind(
  freezing_predictions$AgeSurgery,
  freezing_predictions$DiseaseDurationSurgery,
  factor2ind(freezing_predictions$Sex, "F"),
  factor2ind(freezing_predictions$PostInstOFFmore, "ZeroOne"),
  factor2ind(freezing_predictions$GAITOFFmore, "ZeroOROne"),
  freezing_predictions$UPDRS3off,
  freezing_predictions$UPDRS2off,
  factor2ind(freezing_predictions$LCTrespMORE, "Less50")
)

mod_freezing <- crr(
  freezing_predictions$ftimefreezing,
  freezing_predictions$statusfreezing,
  newfreezing
)
summary(mod_freezing)

# Hallucinations regression
hallucinations_predictions <- read_predictions("hallucinations_predictions.csv")
hallucinations_predictions <- hallucinations_predictions %>%
  mutate(
    Sex             = factor(Sex),
    PostInstOFFmore = factor(PostInstOFFmore),
    GAITOFFmore     = factor(GAITOFFmore),
    LCTrespMORE     = factor(LCTrespMORE),
    Neuropsychologic = factor(Neuropsychologic)
  )

newhallucinations <- cbind(
  hallucinations_predictions$AgeSurgery,
  hallucinations_predictions$DiseaseDurationSurgery,
  factor2ind(hallucinations_predictions$Sex, "F"),
  hallucinations_predictions$UPDRS3off,
  hallucinations_predictions$UPDRS2off,
  hallucinations_predictions$UPDRS1preop,
  factor2ind(hallucinations_predictions$LCTrespMORE, "Less50"),
  factor2ind(hallucinations_predictions$Neuropsychologic, "Normal"),
  hallucinations_predictions$LEDDpreop
)

mod_hallucinations <- crr(
  hallucinations_predictions$ftimehallucinations,
  hallucinations_predictions$statushallucinations,
  newhallucinations
)
summary(mod_hallucinations)

# Dementia regression
dementia_predictions <- read_predictions("dementia_predictions.csv")
dementia_predictions <- dementia_predictions %>%
  mutate(
    Sex             = factor(Sex),
    PostInstOFFmore = factor(PostInstOFFmore),
    GAITOFFmore     = factor(GAITOFFmore),
    LCTrespMORE     = factor(LCTrespMORE),
    Neuropsychologic = factor(Neuropsychologic)
  )

newdementia <- cbind(
  dementia_predictions$AgeSurgery,
  dementia_predictions$DiseaseDurationSurgery,
  factor2ind(dementia_predictions$Sex, "F"),
  dementia_predictions$UPDRS3off,
  dementia_predictions$UPDRS2off,
  dementia_predictions$UPDRS1preop,
  factor2ind(dementia_predictions$LCTrespMORE, "Less50"),
  factor2ind(dementia_predictions$Neuropsychologic, "Normal")
)

mod_dementia <- crr(
  dementia_predictions$ftimedementia,
  dementia_predictions$statusdementia,
  newdementia
)
summary(mod_dementia)

# Institutionalization regression
nursing_predictions <- read_predictions("nursing_predictions.csv")
nursing_predictions <- nursing_predictions %>%
  mutate(
    Sex             = factor(Sex),
    PostInstOFFmore = factor(PostInstOFFmore),
    GAITOFFmore     = factor(GAITOFFmore),
    LCTrespMORE     = factor(LCTrespMORE),
    Neuropsychologic = factor(Neuropsychologic)
  )

newnursing <- cbind(
  nursing_predictions$AgeSurgery,
  nursing_predictions$DiseaseDurationSurgery,
  factor2ind(nursing_predictions$Sex, "F"),
  nursing_predictions$UPDRS3off,
  nursing_predictions$UPDRS2off,
  nursing_predictions$UPDRS1preop,
  factor2ind(nursing_predictions$LCTrespMORE, "Less50"),
  factor2ind(nursing_predictions$Neuropsychologic, "Normal")
)

mod_nursing <- crr(
  nursing_predictions$ftimenursing,
  nursing_predictions$statusnursing,
  newnursing
)
summary(mod_nursing)

## 5. Backward models (final predictors) --------------------

# Falls
newfalls_back <- cbind(
  falls_predictions$AgeSurgery,
  falls_predictions$UPDRS2off,
  falls_predictions$LEDDpreop
)
mod_falls_back <- crr(
  falls_predictions$ftimefalls,
  falls_predictions$statusfalls,
  newfalls_back
)
summary(mod_falls_back)

# Freezing
newfreezing_back <- cbind(
  freezing_predictions$AgeSurgery,
  factor2ind(freezing_predictions$GAITOFFmore, "ZeroOROne"),
  freezing_predictions$UPDRS3on,
  freezing_predictions$UPDRS2off,
  freezing_predictions$UPDRS2on
)
mod_freezing_back <- crr(
  freezing_predictions$ftimefreezing,
  freezing_predictions$statusfreezing,
  newfreezing_back
)
summary(mod_freezing_back)

# Hallucinations
newhallucinations_back <- cbind(
  factor2ind(hallucinations_predictions$Phenotype, "Tremor")
)
mod_hallucinations_back <- crr(
  hallucinations_predictions$ftimehallucinations,
  hallucinations_predictions$statushallucinations,
  newhallucinations_back
)
summary(mod_hallucinations_back)

# Dementia
newdementia_back <- cbind(
  dementia_predictions$AgeSurgery,
  dementia_predictions$DiseaseDurationSurgery,
  dementia_predictions$UPDRS2off,
  dementia_predictions$UPDRS2on,
  factor2ind(dementia_predictions$PostInstOFFmore, "ZeroOne")
)
mod_dementia_back <- crr(
  dementia_predictions$ftimedementia,
  dementia_predictions$statusdementia,
  newdementia_back
)
summary(mod_dementia_back)

# Institutionalization
newnursing_back <- cbind(
  nursing_predictions$AgeSurgery,
  nursing_predictions$UPDRS2off
)
mod_nursing_back <- crr(
  nursing_predictions$ftimenursing,
  nursing_predictions$statusnursing,
  newnursing_back
)
summary(mod_nursing_back)

## 6. Logistic regression & ROC (falls, freezing) ----------

# Falls logistic regression
Multiple_Logistic_Regression_Falls <- read.csv(
  "Multiple_Logistic_Regression_Falls.csv",
  sep = ";",
  header = TRUE
) %>%
  drop_na()

Falls_fit_all <- glm(
  Falls ~ ageatsurgery + diseaseduration + LEDD + UPDRS2ON +
    UPDRS2OFF + UPDRS3OFF + UPDRS3ON + PostInstOFF +
    GaitOFF + LD.response + MMSE,
  data = Multiple_Logistic_Regression_Falls,
  family = binomial,
  na.action = na.pass
)
summary(Falls_fit_all)

roc(
  Multiple_Logistic_Regression_Falls$Falls,
  as.vector(fitted.values(Falls_fit_all)),
  percent      = FALSE,
  boot.n       = 1000,
  ci.alpha     = 0.9,
  stratified   = FALSE,
  plot         = TRUE,
  grid         = TRUE,
  show.thres   = TRUE,
  legacy.axes  = TRUE,
  reuse.auc    = TRUE,
  print.auc    = TRUE,
  ci           = TRUE,
  col          = "deepskyblue4",
  main         = "ROC curve for Falls\n~ All Variables\n"
)

Falls_fit_reduced <- glm(
  Falls ~ UPDRS2OFF + PostInstOFF,
  data = Multiple_Logistic_Regression_Falls,
  family = binomial,
  na.action = na.pass
)
summary(Falls_fit_reduced)

roc(
  Multiple_Logistic_Regression_Falls$Falls,
  as.vector(fitted.values(Falls_fit_reduced)),
  percent      = FALSE,
  boot.n       = 1000,
  ci.alpha     = 0.9,
  stratified   = FALSE,
  plot         = TRUE,
  grid         = TRUE,
  show.thres   = TRUE,
  legacy.axes  = TRUE,
  reuse.auc    = TRUE,
  print.auc    = TRUE,
  ci           = TRUE,
  col          = "deeppink4",
  main         = "ROC curve for Falls\n~ UPDRS2 OFF + Postural Instability\n"
)

# Train-test split: Falls
set.seed(123)
training.samples.falls <- Multiple_Logistic_Regression_Falls$Falls %>%
  createDataPartition(p = 0.8, list = FALSE)

train.data.falls <- Multiple_Logistic_Regression_Falls[training.samples.falls, ]
test.data.falls  <- Multiple_Logistic_Regression_Falls[-training.samples.falls, ]

model_falls_full <- glm(
  Falls ~ .,
  data = train.data.falls,
  family = binomial
)
summary(model_falls_full)

prob_falls <- predict(model_falls_full, test.data.falls, type = "response")
pred_falls <- ifelse(prob_falls > 0.5, 1, 0)
mean(pred_falls == test.data.falls$Falls)

## Freezing logistic regression
Multiple_Logistic_Regression_Freezing <- read.csv(
  "Multiple_Logistic_Regression_Freezing.csv",
  sep = ";",
  header = TRUE
) %>%
  drop_na()

Freezing_fit_all <- glm(
  Freezing ~ ageatsurgery + diseaseduration + LEDD + UPDRS2ON +
    UPDRS2OFF + UPDRS3OFF + UPDRS3ON + PostInstOFF +
    GaitOFF + LD.response + MMSE,
  data = Multiple_Logistic_Regression_Freezing,
  family = binomial,
  na.action = na.pass
)
summary(Freezing_fit_all)

roc(
  Multiple_Logistic_Regression_Freezing$Freezing,
  as.vector(fitted.values(Freezing_fit_all)),
  percent      = FALSE,
  boot.n       = 1000,
  ci.alpha     = 0.9,
  stratified   = FALSE,
  plot         = TRUE,
  grid         = TRUE,
  show.thres   = TRUE,
  legacy.axes  = TRUE,
  reuse.auc    = TRUE,
  print.auc    = TRUE,
  ci           = TRUE,
  col          = "deepskyblue4",
  main         = "ROC curve for Freezing\n~ All Variables\n"
)

Freezing_fit_GAITOFF <- glm(
  Freezing ~ GaitOFF,
  data = Multiple_Logistic_Regression_Freezing,
  family = binomial,
  na.action = na.pass
)
summary(Freezing_fit_GAITOFF)

roc(
  Multiple_Logistic_Regression_Freezing$Freezing,
  as.vector(fitted.values(Freezing_fit_GAITOFF)),
  percent      = FALSE,
  boot.n       = 1000,
  ci.alpha     = 0.9,
  stratified   = FALSE,
  plot         = TRUE,
  grid         = TRUE,
  show.thres   = TRUE,
  legacy.axes  = TRUE,
  reuse.auc    = TRUE,
  print.auc    = TRUE,
  ci           = TRUE,
  col          = "deeppink4",
  main         = "ROC curve for Freezing\n~ GAIT OFF\n"
)

# Train-test split: Freezing
set.seed(123)
training.samples.freezing <- Multiple_Logistic_Regression_Freezing$Freezing %>%
  createDataPartition(p = 0.8, list = FALSE)

train.data.freezing <- Multiple_Logistic_Regression_Freezing[training.samples.freezing, ]
test.data.freezing  <- Multiple_Logistic_Regression_Freezing[-training.samples.freezing, ]

model_freezing_full <- glm(
  Freezing ~ .,
  data = train.data.freezing,
  family = binomial
)
summary(model_freezing_full)

prob_freezing <- predict(model_freezing_full, test.data.freezing, type = "response")
pred_freezing <- ifelse(prob_freezing > 0.5, 1, 0)
mean(pred_freezing == test.data.freezing$Freezing)

## 7. Correlation matrix ------------------------------------

Correlation_Matrix_Input_Data <- read.csv(
  "Correlation_Matrix_Input_Data.csv",
  sep = ";",
  header = TRUE
)

cor_data  <- rcorr(as.matrix(Correlation_Matrix_Input_Data), type = "spearman")
corr_vals <- as.data.frame(cor_data$r)
p_vals    <- as.data.frame(cor_data$P)

ggcorrplot(
  corr_vals,
  method = "circle",
  colors = c("yellow", "white", "darkslateblue"),
  outline.color = "white",
  hc.order = TRUE,
  type = "lower",
  lab = TRUE
)

ggcorrplot(
  p_vals,
  method = "circle",
  colors = c("white", "white", "cadetblue2"),
  outline.color = "white",
  hc.order = TRUE,
  type = "lower",
  lab = TRUE,
  show.legend = FALSE
)

## 8. Plotting cumulative incidence & survival --------------

# Per-event cumulative incidence
Milestones_Summary <- read.csv("Milestones_Summary.csv", sep = ";", header = TRUE)
event_order <- c("Falls", "Freezing", "Dementia", "Hallucinations", "Institutionalization")

Milestones_Summary <- Milestones_Summary %>%
  mutate(Event = factor(Event, levels = event_order))

Milestones_Summary %>%
  ggplot(aes(x = Follow_up_months, y = Proportion)) +
  geom_step(aes(color = Event), show.legend = FALSE) +
  geom_stepconfint(aes(ymin = conf.low, ymax = conf.high, fill = Event), alpha = 0.6) +
  ylim(0, 1) +
  labs(x = "\nFollow-up (Months)", y = "Cumulative Incidence\n") +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  theme_minimal() +
  facet_wrap(~Event, scales = "free")

# All events in single panel
Milestones_Each <- read.csv("Milestones_Each.csv", sep = ";", header = TRUE) %>%
  mutate(Event = factor(Event, levels = event_order))

Milestones_Each %>%
  mutate(Follow_up_months = log10(Follow_up_months)) %>%
  ggplot(aes(x = Follow_up_months, y = Proportion)) +
  geom_line(aes(color = Event), show.legend = FALSE, size = 2) +
  ylim(0, 1) +
  labs(x = NULL, y = NULL) +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  theme_minimal()

# Paired groups (with vs without event)
Milestones_Grouped <- read.csv("Milestones_Grouped.csv", sep = ";", header = TRUE)
Milestones_Grouped[[3]] <- as.factor(Milestones_Grouped[[3]])
Milestones_Grouped[[4]] <- as.factor(Milestones_Grouped[[4]])

Milestones_Grouped %>%
  ggplot(aes(x = Follow_up_months, y = Proportion)) +
  geom_step(aes(color = Group), show.legend = FALSE, size = 1.5) +
  geom_stepconfint(aes(ymin = conf.low, ymax = conf.high, fill = Group), alpha = 0.4) +
  ylim(0, 1) +
  labs(x = "\nFollow-up (Months)", y = "Cumulative Incidence\n") +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  theme_minimal() +
  facet_wrap(~Event, scales = "free")

# Survival by milestone group
Milestones_Mortality_per_Group <- read.csv(
  "Milestones_Mortality_per_Group.csv",
  sep = ";",
  header = TRUE
)

Milestones_Mortality_per_Group[[3]] <- as.factor(Milestones_Mortality_per_Group[[3]])
Milestones_Mortality_per_Group[[4]] <- as.factor(Milestones_Mortality_per_Group[[4]])

event_order_2 <- c(
  "Overall_Survival", "Survival_Freezing", "Survival_Falls",
  "Survival_Hallucinations", "Survival_Dementia", "Survival_Institutionalization"
)

Milestones_Mortality_per_Group <- Milestones_Mortality_per_Group %>%
  mutate(Event = factor(Event, levels = event_order_2))

Milestones_Mortality_per_Group %>%
  group_by(Event, Group) %>%
  fill(Proportion) %>%
  fill(conf.high) %>%
  fill(conf.high, .direction = "up") %>%
  fill(conf.low) %>%
  fill(conf.low, .direction = "up") %>%
  ungroup() %>%
  ggplot(aes(x = Follow_up_months, y = Proportion)) +
  geom_step(aes(color = Group), show.legend = FALSE, size = 2) +
  geom_stepconfint(aes(ymin = conf.low, ymax = conf.high, fill = Group), alpha = 0.2) +
  ylim(0, 100) +
  xlim(0, 110) +
  labs(x = "\nFollow-up (Months)", y = "Survival (%)\n") +
  scale_fill_viridis_d(option = "D") +
  scale_colour_viridis_d(option = "D") +
  theme_minimal() +
  facet_wrap(~Event, scales = "free")

## 9. Drug usage summary ------------------------------------

MilestonesDrugUsage <- read.csv("MilestonesDrugUsage.csv", sep = ";", header = TRUE)

drug_nonzero_counts <- sapply(
  MilestonesDrugUsage[, 4:17],
  function(x) sum(x != 0, na.rm = TRUE)
)
print(drug_nonzero_counts)

## 10. Stimulation settings summary ------------------------

Stimul_Settings <- read.csv("Stimul_Settings.csv", sep = ";", header = TRUE)
Stimul_Settings <- Stimul_Settings[, 3:10]

# Right side
Stimul_Settings %>%
  group_by(modeSTNrightEND) %>%
  summarise(
    n           = n(),
    mean_V      = mean(voltageSTNrightEND, na.rm = TRUE),
    sd_V        = sd(voltageSTNrightEND, na.rm = TRUE),
    mean_freq   = mean(frequencySTNrightEND, na.rm = TRUE),
    sd_freq     = sd(frequencySTNrightEND, na.rm = TRUE),
    mean_pw     = mean(pulsewidthSTNrightEND, na.rm = TRUE),
    sd_pw       = sd(pulsewidthSTNrightEND, na.rm = TRUE),
    .groups     = "drop"
  )

# Left side
Stimul_Settings %>%
  group_by(modeSTNLeftEND) %>%
  summarise(
    n           = n(),
    mean_V      = mean(voltageSTNLeftEND, na.rm = TRUE),
    sd_V        = sd(voltageSTNLeftEND, na.rm = TRUE),
    mean_freq   = mean(frequencySTNLeftEND, na.rm = TRUE),
    sd_freq     = sd(frequencySTNLeftEND, na.rm = TRUE),
    mean_pw     = mean(pulsewidthSTNLeftEND, na.rm = TRUE),
    sd_pw       = sd(pulsewidthSTNLeftEND, na.rm = TRUE),
    .groups     = "drop"
  )
