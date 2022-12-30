---
name: Survival
topic: Survival Analysis
maintainer: Arthur Allignol, Aurelien Latouche
email: arthur.allignol@gmail.com
version: 2022-03-07
source: https://github.com/cran-task-views/Survival/
---


Survival analysis, also called event history analysis in social science,
or reliability analysis in engineering, deals with time until occurrence
of an event of interest. However, this failure time may not be observed
within the relevant time period, producing so-called censored
observations.

This task view aims at presenting the useful R packages for the analysis
of time to event data.

Please let the maintainers know if something is inaccurate or missing,
either via e-mail or by submitting an issue or pull request in the GitHub
repository linked above.


Standard Survival Analysis
==========================

### Estimation of the Survival Distribution

-   ***Kaplan-Meier:*** The `survfit` function from the
    `r pkg("survival", priority = "core")` package computes
    the Kaplan-Meier estimator for truncated and/or censored data.
    `r pkg("rms", priority = "core")` (replacement of the
    Design package) proposes a modified version of the `survfit`
    function. The `r pkg("prodlim")` package implements a
    fast algorithm and some features not included in
    `r pkg("survival")`. Various confidence intervals and
    confidence bands for the Kaplan-Meier estimator are implemented in
    the `r pkg("km.ci")` package. `plot.Surv` of package
    `r pkg("eha", priority = "core")` plots the Kaplan-Meier
    estimator. The `r pkg("NADA")` package includes a
    function to compute the Kaplan-Meier estimator for left-censored
    data. `svykm` in `r pkg("survey")` provides a weighted
    Kaplan-Meier estimator. The `kaplan-meier` function in
    `r pkg("spatstat")` computes the Kaplan-Meier estimator
    from histogram data. The `KM` function in package
    `r pkg("rhosp")` plots the survival function using a
    variant of the Kaplan-Meier estimator in a hospitalisation risk
    context. The `r pkg("survPresmooth")` package computes
    presmoothed estimates of the main quantities used for right-censored
    data, i.e., survival, hazard and density functions. The
    `r pkg("asbio")` package permits to compute the
    Kaplan-Meier estimator following Pollock et al. (1998). The
    `r pkg("bpcp")` package provides several functions for
    computing confidence intervals of the survival distribution (e.g.,
    beta product confidence procedure). The
    `r pkg("kmc")` package implements the Kaplan-Meier
    estimator with constraints. The `r pkg("landest")`
    package allows landmark estimation and testing of survival
    probabilities. The `r pkg("jackknifeKME")` package
    computes the original and modified jackknife estimates of
    Kaplan-Meier estimators. The `r pkg("tranSurv")` package
    permits to estimate a survival distribution in the presence of
    dependent left-truncation and right-censoring. The
    `r pkg("condSURV")` package provides methods for
    estimating the conditional survival function for ordered
    multivariate failure time data. The `r pkg("gte")`
    package implements the generalised Turnbull estimator proposed by
    Dehghan and Duchesne for estimating the conditional survival
    function with interval-censored data.
-   ***Non-Parametric maximum likelihood estimation (NPMLE):*** The  
    `r bioc("Icens")` package provides several ways to compute the NPMLE
    of the survival distribution for various censoring and truncation
    schemes. `r pkg("MLEcens")` can also be used to compute the MLE
    for interval-censored data. `r pkg("dblcens")` permits to compute
    the NPMLE of the cumulative distribution function for left- and
    right-censored data.  The `icfit` function in package 
    `r pkg("interval")` computes the NPMLE for interval-censored
    data. The `r pkg("DTDA")` package implements several algorithms
    permitting to analyse possibly doubly truncated survival data.  
    `r pkg("npsurv")` computes the NPMLE of a survival function for
    general interval-censored data. The `r pkg("csci")` package
    provides confidence intervals for the cumulative distribution
    function of the event time for current status data, including a
    new method that is valid (i.e., exact).

-   ***Parametric:*** The `r pkg("fitdistrplus")` package
    permits to fit an univariate distribution by maximum likelihood.
    Data can be interval censored. The `r pkg("vitality")`
    package provides routines for fitting models in the vitality family
    of mortality models.

### Hazard Estimation

-   The `r pkg("muhaz", priority = "core")` package permits
    to estimate the hazard function through kernel methods for
    right-censored data.
-   The `epi.insthaz` function from `r pkg("epiR")` computes
    the instantaneous hazard from the Kaplan-Meier estimator.
-   `r pkg("polspline")`, `r pkg("gss")` and
    `r pkg("logspline")` allow to estimate the hazard
    function using splines.
-   The `r pkg("bshazard")` package provides non-parametric
    smoothing of the hazard through B-splines.

### Testing

-   The `survdiff` function in `r pkg("survival")` compares
    survival curves using the Fleming-Harrington G-rho family of test.
    `r pkg("NADA")` implements this class of tests for
    left-censored data.
-   The `r pkg("maxcombo")` package compares survival curves using the
    max-combo test, which is often based on the Fleming-Harrington
    G-rho family of tests and is designed to have higher power than
    the logrank test in the scenario of non-proportional hazards such
    as those resulting from delayed treatment effects.
-   `r pkg("clinfun")` implements a permutation version of
    the logrank test and a version of the logrank that adjusts for
    covariates.
-   The `r pkg("exactRankTests")` implements the
    shift-algorithm by Streitberg and Roehmel for computing exact
    conditional p-values and quantiles, possibly for censored data.
-   `SurvTest` in the `r pkg("coin")` package implements the
    logrank test reformulated as a linear rank test.
-   The `r pkg("maxstat")` package performs tests using
    maximally selected rank statistics.
-   The `r pkg("interval")` package implements logrank and
    Wilcoxon type tests for interval-censored data.
-   `r bioc("survcomp")` compares 2 hazard ratios.
-   The `r pkg("TSHRC")` implements a two stage procedure
    for comparing hazard functions.
-   The `r pkg("FHtest")` package offers several tests based
    on the Fleming-Harrington class for comparing surival curves with
    right- and interval-censored data.
-   The short term and long term hazard ratio model for two samples
    survival data can be found in the `r pkg("YPmodel")`
    package.
-   The `r pkg("controlTest")` implements a nonparametric
    two-sample procedure for comparing the median survival time.
-   The `r pkg("survRM2")` package performs two-sample
    comparison of the restricted mean survival time
-   The `r pkg("emplik2")` package permits to compare two
    samples with censored data using empirical likelihood ratio tests.
-   The `r pkg("KONPsurv")` package provides powerful nonparametric
    K-sample tests for right-censored data. The tests are consistent
    against any differences between the hazard functions of the groups.

### Regression Modelling

-   ***Cox model:*** The `coxph` function in the
    `r pkg("survival")` package fits the Cox model. `cph` in
    the `r pkg("rms")` package and the
    `r pkg("eha")` package propose some extensions to the
    `coxph` function. The package `r pkg("coxphf")`
    implements the Firth's penalised maximum likelihood bias reduction
    method for the Cox model. An implementation of weighted estimation
    in Cox regression can be found in `r pkg("coxphw")`. The
    `r pkg("coxrobust")` package proposes a robust
    implementation of the Cox model. `timecox` in package
    `r pkg("timereg", priority = "core")` fits Cox models
    with possibly time-varying effects. 
    A Cox model model can be fitted to data from complex
    survey design using the `svycoxph` function in
    `r pkg("survey")`. The
    `r pkg("multipleNCC")` package fits Cox models using a
    weighted partial likelihood for nested case-control studies. The
    `r pkg("MIICD")` package implements Pan's (2000)
    multiple imputation approach to Cox models for interval censored
    data. The `r pkg("ICsurv")` package fits Cox models for
    interval-censored data through an EM algorithm. The
    `r pkg("dynsurv")` package fits time-varying coefficient
    models for interval censored and right censored survival data using
    a Bayesian Cox model, a spline based Cox model or a transformation
    model. The `r pkg("OrdFacReg")` package implements the
    Cox model using an active set algorithm for dummy variables of
    ordered factors. The `r pkg("survivalMPL")` package fits
    Cox models using maximum penalised likelihood and provide a non
    parametric smooth estimate of the baseline hazard function. A Cox
    model with piecewise constant hazards can be fitted using the
    `r pkg("pch")` package. The
    `r pkg("icenReg")` package implements several models for
    interval-censored data, e.g., Cox, proportional odds, and
    accelerated failure time models. A Cox type Self-Exciting Intensity
    model can be fitted to right-censored data using the
    `r pkg("coxsei")` package. The
    `r pkg("SurvLong")` contains methods for estimation of
    proportional hazards models with intermittently observed
    longitudinal covariates. The `r pkg("plac")` package
    provides routines to fit the Cox model with left-truncated data
    using augmented information from the marginal of the truncation
    times.\
    The proportionality assumption can be checked using the `cox.zph`
    function in `r pkg("survival")`. The `coxphCPE` function
    in `r pkg("clinfun")` calculates concordance probability
    estimate for the Cox model. The `coxphQuantile` in the latter
    package draws a quantile curve of the survival distribution as a
    function of covariates. The `r pkg("multcomp")` package
    computes simultaneous tests and confidence intervals for the Cox
    model and other parametric survival models. The
    `r pkg("lsmeans")` package permits to obtain
    least-squares means (and contrasts thereof) from linear models. In
    particular, it provides support for the `coxph`, `survreg` and
    `coxme` functions. The `r bioc("multtest")` package on
    Bioconductor proposes a resampling based multiple hypothesis testing
    that can be applied to the Cox model. Testing coefficients of Cox
    regression models using a Wald test with a sandwich estimator of
    variance can be done using the `r pkg("saws")` package.
    The `r pkg("rankhazard")` package permits to plot
    visualisation of the relative importance of covariates in a
    proportional hazards model. The `r pkg("smoothHR")`
    package provides hazard ratio curves that allows for nonlinear
    relationship between predictor and survival. The `r pkg("PHeval")`
    package proposes tools to check the proportional hazards assumption
    using a standardised score process. The `r pkg("ELYP")`
    package implements empirical likelihood analysis for the Cox Model
    and Yang-Prentice (2005) Model.
-   ***Parametric Proportional Hazards Model:*** `survreg` (from
    `r pkg("survival")`) fits a parametric proportional
    hazards model. The `r pkg("eha")` and
    `r pkg("mixPHM")` packages implement a proportional
    hazards model with a parametric baseline hazard. The `pphsm` in
    `r pkg("rms")` translates an AFT model to a proportional
    hazards form. The `r pkg("polspline")` package includes
    the `hare` function that fits a hazard regression model, using
    splines to model the baseline hazard. Hazards can be, but not
    necessarily, proportional. The `r pkg("flexsurv")`
    package implements the model of Royston and Parmar (2002). The model
    uses natural cubic splines for the baseline survival function, and
    proportional hazards, proportional odds or probit functions for
    regression. The `r pkg("SurvRegCensCov")` package allows
    estimation of a Weibull Regression for a right-censored endpoint,
    one interval-censored covariate, and an arbitrary number of
    non-censored covariates.
-   ***Accelerated Failure Time (AFT) Models:*** The `survreg` function
    in package `r pkg("survival")` can fit an accelerated
    failure time model. A modified version of `survreg` is implemented
    in the `r pkg("rms")` package ( `psm` function). It
    permits to use some of the `r pkg("rms")`
    functionalities. The `r pkg("eha")` package also
    proposes an implementation of the AFT model (function `aftreg`). An
    AFT model with an error distribution assumed to be a mixture of
    G-splines is implemented in the `r pkg("smoothSurv")`
    package. The `r pkg("NADA")` package proposes the front
    end of the `survreg` function for left-censored data. The
    `r pkg("simexaft")` package implements the
    Simulation-Extrapolation algorithm for the AFT model, that can be
    used when covariates are subject to measurement error. A robust
    version of the accelerated failure time model can be found in
    `r pkg("RobustAFT")`. The
    `r pkg("coarseDataTools")` package fits AFT models for
    interval censored data. An alternative weighting scheme for
    parameter estimation in the AFT model is proposed in the
    `r pkg("imputeYn")` package. The `r pkg("aftgee")` package
    implements recently developed inference procedures for the AFT
    models with both the rank-based approach and the least squares approach.
    `r pkg("imputeYn")` package.
-   ***Additive Models:*** Both `r pkg("survival")` and
    `r pkg("timereg")` fit the additive hazards model of
    Aalen in functions `aareg` and `aalen`, respectively.
    `r pkg("timereg")` also proposes an implementation of
    the Cox-Aalen model (that can also be used to perform the Lin, Wei
    and Ying (1994) goodness-of-fit for Cox regression models) and the
    partly parametric additive risk model of McKeague and Sasieni. A
    version of the Cox-Aalen model for interval censored data is
    available in the `r pkg("coxinterval")` package. The
    `r pkg("uniah")` package fits shape-restricted additive
    hazards models. The `r pkg("addhazard")` package
    contains tools to fit additive hazards model to random sampling,
    two-phase sampling and two-phase sampling with auxiliary
    information.
-   ***Buckley-James Models:*** The `bj` function in
    `r pkg("rms")` and `BJnoint` in
    `r pkg("emplik")` compute the Buckley-James model,
    though the latter does it without an intercept term. The
    `r pkg("bujar")` package fits the Buckley-James model
    with high-dimensional covariates (L2 boosting, regression trees and
    boosted MARS, elastic net).
-   ***Other models:*** Functions like `survreg` can fit other types of
    models depending on the chosen distribution, *e.g.* , a tobit model.
    The `r pkg("AER")` package provides the `tobit`
    function, which is a wrapper of `survreg` to fit the tobit model. An
    implementation of the tobit model for cross-sectional data and panel
    data can be found in the `r pkg("censReg")` package. The
    `r pkg("timereg")` package provides implementation of
    the proportional odds model and of the proportional excess hazards
    model. The `r pkg("invGauss")` package fits the inverse
    Gaussian distribution to survival data. The model is based on
    describing time to event as the barrier hitting time of a Wiener
    process, where drift towards the barrier has been randomized with a
    Gaussian distribution. The `r pkg("pseudo")` package
    computes the pseudo-observation for modelling the survival function
    based on the Kaplan-Meier estimator and the restricted mean. The
    `r pkg("fastpseudo")` package dose the same for the
    restricted mean survival time. `r pkg("flexsurv")` fits
    parametric time-to-event models, in which any parametric
    distribution can be used to model the survival probability, and
    where one of the parameters is a linear function of covariates. The
    `Icens` function in package `r pkg("Epi")` provides a
    multiplicative relative risk and an additive excess risk model for
    interval-censored data. The `r pkg("VGAM")` package can
    fit vector generalised linear and additive models for censored data.
    The `r pkg("gamlss.cens")` package implements the
    generalised additive model for location, scale and shape that can be
    fitted to censored data. The `locfit.censor` function in
    `r pkg("locfit")` produces local regression estimates.
    The `crq` function included in the `r pkg("quantreg")`
    package implements a conditional quantile regression model for
    censored data. The `r pkg("JM")` package fits shared
    parameter models for the joint modelling of a longitudinal response
    and event times. The temporal process regression model is
    implemented in the `r pkg("tpr")` package. Aster models,
    which combine aspects of generalized linear models and Cox models,
    are implemented in the `r pkg("aster")` and
    `r pkg("aster2")` packages. The
    `r pkg("concreg")` package implements conditional
    logistic regression for survival data as an alternative to the Cox
    model when hazards are non-proportional. The
    `r pkg("surv2sampleComp")` packages proposes some
    model-free contrast comparison measures such as difference/ratio of
    cumulative hazards, quantiles and restricted mean. The
    `r pkg("rstpm2")` package provides link-based survival
    models that extend the Royston-Parmar models, a family of flexible
    parametric models. The `r pkg("TransModel")` package
    implements a unified estimation procedure for the analysis of
    censored data using linear transformation models. The
    `r pkg("ICGOR")` fits the generalized odds rate hazards
    model to interval-censored data while `r pkg("GORCure")`
    generalized odds rate mixture cure model to interval-censored data.
    The `r pkg("thregI")` package permits to fit a threshold
    regression model for interval-censored data based on the
    first-hitting-time of a boundary by the sample path of a Wiener
    diffusion process. The `r pkg("miCoPTCM")` package fits
    semiparametric promotion time cure models with possibly mis-measured
    covariates. The `r pkg("smcure")` package permits to fit
    semiparametric proportional hazards and accelerated failure time
    mixture cure models. The case-base sampling approach for fitting
    flexible hazard regression models to survival data with single event
    type or multiple competing causes via logistic and multinomial
    regression can be found in package `r pkg("casebase")`.
    The `r pkg("intsurv")` package fits regular Cox cure rate model via
    an EM algorithm, regularized Cox cure rate model with elastic net
     penalty, and weighted concordance index for cure models.


Multistate Models
=================

-   ***General Multistate Models:*** The `coxph` function from package
    `r pkg("survival")` can be fitted for any transition of
    a multistate model. It can also be used for comparing two transition
    hazards, using correspondence between multistate models and
    time-dependent covariates. Besides, all the regression methods
    presented above can be used for multistate models as long as they
    allow for left-truncation.\
    The `r pkg("mvna")` package provides convenient
    functions for estimating and plotting the cumulative transition
    hazards in any multistate model, possibly subject to right-censoring
    and left-truncation. The `r pkg("etm")` package
    estimates and plots transition probabilities for any multistate
    models. It can also estimate the variance of the Aalen-Johansen
    estimator, and handles left-truncated data. The
    `r pkg("mstate", priority = "core")` package permits to
    estimate hazards and probabilities, possibly depending on
    covariates, and to obtain prediction probabilities in the context of
    competing risks and multistate models. The
    `r pkg("msm")` package contains functions for fitting
    general continuous-time Markov and hidden Markov multistate models
    to longitudinal data. Transition rates and output processes can be
    modelled in terms of covariates. The `r pkg("msmtools")`
    package provides utilities to facilitate the modelling of
    longitudinal data under a multistate framework using the
    `r pkg("msm")` package.The
    `r pkg("SemiMarkov")` package can be used to fit
    semi-Markov multistate models in continuous time. The distribution
    of the waiting times can be chosen between the exponential, the
    Weibull and exponentiated Weibull distributions. The
    `r pkg("TPmsm")` package permits to estimate transition
    probabilities of an illness-death model or three-state progressive
    model. The `r pkg("gamboostMSM")` package extends the
    `r pkg("mboost")` package to estimation in the
    mulstistate model framework, while the `r pkg("penMSM")`
    package proposes L1 penalised estimation.
    The `r pkg("SmoothHazard")` package fits proportional
    hazards models for the illness-death model with possibly
    interval-censored data for transition toward the transient state.
    Left-truncated and right-censored data are also allowed. The model
    is either parametric (Weibull) or semi-parametric with M-splines
    approximation of the baseline intensities. The
    `r pkg("TP.idm")` package implement the estimator of
    Una-Alvarez and Meira-Machado (2015) for non-Markov illness-death
    models.\
    The `r pkg("Epi")` package implements Lexis objects as a
    way to represent, manipulate and summarise data from multistate
    models. The `r pkg("LexisPlotR")` package, based on
    **ggplot2** , permits to draw Lexis diagrams. The
    `r pkg("TraMineR")` package is intended for analysing
    state or event sequences that describe life courses.
    `r pkg("asbio")` computes the expected numbers of
    individuals in specified age classes or life stages given
    survivorship probabilities from a transition matrix.
-   ***Competing risks:*** The package
    `r pkg("cmprsk", priority = "core")` estimates the
    cumulative incidence functions, but they can be compared in more
    than two samples. The package also implements the Fine and Gray
    model for regressing the subdistribution hazard of a competing risk.
    `r pkg("crrSC")` extends the
    `r pkg("cmprsk")` package to stratified and clustered
    data. The `r pkg("kmi")` package performs a Kaplan-Meier
    multiple imputation to recover missing potential censoring
    information from competing risks events, permitting to use standard
    right-censored methods to analyse cumulative incidence functions. Package
    `r pkg("pseudo")` computes pseudo observations for
    modelling competing risks based on the cumulative incidence
    functions. `r pkg("timereg")` does flexible regression
    modelling for competing risks data based on the on the
    inverse-probability-censoring-weights and direct binomial regression
    approach. `r pkg("riskRegression")` implements risk
    regression for competing risks data, along with other extensions of
    existing packages useful for survival analysis and competing risks
    data. The `r pkg("Cprob")` package estimates the
    conditional probability of a competing event, aka., the conditional
    cumulative incidence. It also implements a proportional-odds model
    using either the temporal process regression or the pseudo-value
    approaches. Packages `r pkg("survival")` (via `survfit`)
    and `r pkg("prodlim")` can also be used to estimate the
    cumulative incidence function. The
    `r pkg("NPMLEcmprsk")` package implements the
    semi-parametric mixture model for competing risks data.
	The `r pkg("CFC")` package permits to
    perform Bayesian, and non-Bayesian, cause-specific competing risks
    analysis for parametric and non-parametric survival functions. The
    `r pkg("gcerisk")` package provides some methods for
    competing risks data. Estimation, testing and regression modeling of
    subdistribution functions in the competing risks setting using
    quantile regressions can be had in `r pkg("cmprskQR")`.
    The `r pkg("intccr")` package permits to fit the Fine
    and Gray model as well other models that belong to the class of
    semiparametric generalized odds rate transformation models to
    interval-censored competing risks data.
-   ***Recurrent event data:*** `coxph` from the
    `r pkg("survival")` package can be used to analyse
    recurrent event data. The `cph` function of the
    `r pkg("rms")` package fits the Anderson-Gill model for
    recurrent events, model that can also be fitted with the
    `r pkg("frailtypack")` package. The latter also permits
    to fit joint frailty models for joint modelling of recurrent events
    and a terminal event. The `r pkg("condGEE")` package
    implements the conditional GEE for recurrent event gap times. The
    `r pkg("reda")` package provides function to fit gamma
    frailty model with either a piecewise constant or a spline as the
    baseline rate function for recurrent event data, as well as some
    miscellaneous functions for recurrent event data. Several regression
    models for recurrent event data are implemented in the
    `r pkg("reReg")` package. The
    `r pkg("spef")` package includes functions for fitting
    semiparametric regression models for panel count survival data.

Relative Survival
=================

-   The `r pkg("relsurv")` package proposes several
    functions to deal with relative survival data. For example,
    `rs.surv` computes a relative survival curve. `rs.add` fits an
    additive model and `rsmul` fits the Cox model of Andersen et al. for
    relative survival, while `rstrans` fits a Cox model in transformed
    time.
-   The `r pkg("timereg")` package permits to fit relative
    survival models like the proportional excess and additive excess
    models.
-   The `r pkg("mexhaz")` package allows fitting an hazard
    regression model using different shapes for the baseline hazard. The
    model can be used in the relative survival setting (excess mortality
    hazard) as well as in the overall survival setting (overall
    mortality hazard).
-   The `r pkg("flexrsurv")` package implements the models
    of Remontet et al. (2007) and Mahboubi et al. (2011).
-   The `r pkg("survexp.fr")` package computes relative
    survival, absolute excess risk and standardized mortality ratio
    based on French death rates.
-   The `r pkg("MRsurv")` package permits to fit
    multiplicative regression models for relative survival.
-   The `r pkg("GJRM")` package permits to fit link-based additive models for the excess hazard that allows for the inclusion of many types of covariate effects, including spatial and time-dependent effects, using any type of smoother, such as thin plate, cubic splines, tensor products and Markov random fields.

Random Effect Models
====================

-   ***Frailties:*** Frailty terms can be added in `coxph` and `survreg`
    functions in package `r pkg("survival")`. A
    mixed-effects Cox model is implemented in the
    `r pkg("coxme")` package. The `two.stage` function in
    the `r pkg("timereg")` package fits the
    Clayton-Oakes-Glidden model. The `r pkg("frailtypack")` package
    fits proportional hazards models with a shared Gamma frailty to
    right-censored and/or left-truncated data using a penalised
    likelihood on the hazard function. The package also fits additive
    and nested frailty models that can be used for, e.g., meta-analysis
    and for hierarchically clustered data (with 2 levels of clustering),
    respectively. The Cox model using
    h-likelihood estimation for the frailty terms can be fitted using
    the `r pkg("frailtyHL")` package. The `r pkg("frailtySurv")` package
    simulates and fits semiparametric shared frailty models under a wide
    range of frailty distributions. The `r pkg("PenCoxFrail")` package provides
    a regularisation approach for Cox frailty models through
    penalisation. The `r pkg("mexhaz")` enables modelling of
    the excess hazard regression model with time-dependent and/or
    non-linear effect(s) and a random effect defined at the cluster
    level. The `r pkg("frailtyEM")` package contains
    functions for fitting shared frailty models with a semi-parametric
    baseline hazard with the Expectation-Maximization algorithm.
    Supported data formats include clustered failures with left
    truncation and recurrent events in gap-time or Andersen-Gill format
-   ***Joint modelling of time-to-event and longitudinal data:*** The
    `r pkg("joineR")` package allows the analysis of
    repeated measurements and time-to-event data via joint random
    effects models. The `r pkg("joint.Cox")` package
    performs Cox regression and dynamic prediction under the joint
    frailty-copula model between tumour progression and death for
    meta-analysis. `r pkg("JointModel")` fits semiparametric
    regression model for longitudinal responses and a semiparametric
    transformation model for time-to-event data. The
    `r pkg("joineRML")` package fits the joint model
    proposed by Henderson and colleagues (2000)
    [doi:10.1093/biostatistics/1.4.465](http://dx.doi.org/10.1093/biostatistics/1.4.465)
    , but extended to the case of multiple continuous longitudinal
    measures. The `r pkg("rstanarm")` package fits joint
    models for one or more longitudinal outcomes (continuous, binary or
    count data) and a time-to-event, estimated under a Bayesian
    framework.

Multivariate Survival
=====================

Multivariate survival refers to the analysis of unit, e.g., the survival
of twins or a family. To analyse such data, we can estimate the joint
distribution of the survival times

-   ***Joint modelling:*** Both `r bioc("Icens")` and
    `r pkg("MLEcens")` can estimate bivariate survival data
    subject to interval censoring.
-   The `r pkg("mets")` package implements various
    statistical models for multivariate event history data, e.g.,
    multivariate cumulative incidence models, bivariate random effects
    probit models, Clayton-Oakes model.
-   The `r pkg("MST")` package constructs trees for
    multivariate survival data using marginal and frailty models.
-   The `r pkg("SurvCorr")` package permits to estimate
    correlation coefficients with associated confidence limits for
    bivariate, partially censored survival times.

Bayesian Models
===============

-   The `r pkg("bayesSurv")` package proposes an
    implementation of a bivariate AFT model.
-   The package `r pkg("BMA")` computes a Bayesian model
    averaging for Cox proportional hazards models.
-   `NMixMCMC` in `r pkg("mixAK")` performs an MCMC
    estimation of normal mixtures for censored data.
-   A MCMC for Gaussian linear regression with left-, right- or
    interval-censored data can be fitted using the `MCMCtobit` in
    `r pkg("MCMCpack")`.
-   The `weibullregpost` function in `r pkg("LearnBayes")`
    computes the log posterior density for a Weibull proportional-odds
    regression model.
-   The `r pkg("MCMCglmm")` fits generalised linear mixed
    models using MCMC to right-, left- and interval censored data.
-   The `r pkg("BaSTA")` package aims at drawing inference
    on age-specific mortality from capture-recapture/recovery data when
    some or all records have missing information on times of birth and
    death. Covariates can also be included in the model.
-   The `r pkg("JMbayes")` package performs joint modelling
    of longitudinal and time-to-event data under a bayesian approach.
-   The `r pkg("rstanarm")` package fits a joint model for
    one or more longitudinal outcomes (continuous, binary or count data)
    and a time-to-event under a Bayesian framework.
-   Bayesian parametric and semi-parametric estimation for
    semi-competing risks data is available via the
    `r pkg("SemiCompRisks")` package.
-   The `r pkg("psbcGroup")` package implements penalized
    semi-parametric Bayesian Cox models with elastic net, fused lasso
    and group lasso priors.
-   The `r pkg("PReMiuM")` package implements Bayesian
    clustering using a Dirichlet process mixture model to censored
    responses.
-   The `r pkg("spBayesSurv")` package provides Bayesian
    model fitting for several survival models including spatial copula,
    linear dependent Dirichlet process mixture model, anova Dirichlet
    process mixture model, proportional hazards model and marginal
    spatial proportional hazards model.
-   The `r pkg("IDPSurvival")` package implements
    non-parametric survival analysis techniques using a prior
    near-ignorant Dirichlet Process.
-   The `r pkg("ICBayes")` packages permits to fit Bayesian
    semiparametric regression survival models (proportional hazards
    model, proportional odds model, and probit model) to
    interval-censored time-to-event data
-   The `r pkg("BayesPiecewiseICAR")` package fits a
    piecewise exponential hazard to survival data using a Hierarchical
    Bayesian model.

Machine learning
================

-   ***Recursive partitioning:*** `r pkg("rpart")`
    implements CART-like trees that can be used with censored outcomes.
    The `r pkg("party")` package implements recursive
    partitioning for survival data. `r pkg("LogicReg")` can
    perform logic regression. 
    The `r pkg("LTRCtrees")` package provides recursive
    partition algorithms designed for fitting survival tree with
    left-truncated and right censored data. The package also includes an
    implementation of recursive partitioning (conditional inference
    trees) for interval-censored data.
    `r pkg("bnnSurvival")` implements a bootstrap aggregated
    version of the k-nearest neighbors survival probability prediction
    method.
-   ***Random forest:*** Package `r pkg("ipred")` implements
    bagging for survival data. The
    `r pkg("randomForestSRC")` package fits random forest to
    survival data, while a variant of the random forest is implemented
    in `r pkg("party")`. A faster implementation can be
    found in package `r pkg("ranger")`. An alternative
    algorithm for random forests is implemented in
    `r pkg("icRSF")`.
-   ***Regularised and shrinkage methods:*** The
    `r pkg("glmnet")` package provides procedures for
    fitting the entire lasso or elastic-net regularization path for Cox
    models. The `r pkg("glmpath")` package implements a L1
    regularised Cox proportional hazards model. An L1 and L2 penalised
    Cox models are available in `r pkg("penalized")`. The
    `r pkg("pamr")` package computes a nearest shrunken
    centroid for survival gene expression data. The
    `r pkg("lpc")` package implements the lassoed principal
    components method. The `r pkg("ahaz")` package
    implements the LASSO and elastic net estimator for the additive risk
    model. The `r pkg("fastcox")` package implements the
    Lasso and elastic-net penalized Cox's regression using the cockail
    algorithm. The `r pkg("SGL")` package permits to fit Cox
    models with a combination of lasso and group lasso regularisation.
    The `r pkg("hdnom")` package implements 9 types of
    penalised Cox regression methods and provides methods for model
    validation, calibration, comparison, and nomogram visualisation.
	The `r pkg("Cyclops")`
    package implements cyclic coordinate descent for the Cox
    proportional hazards model.
-   ***Boosting:*** Gradient boosting for the Cox model is implemented
    in the `r pkg("gbm")` package. The
    `r pkg("mboost")` package includes a generic gradient
    boosting algorithm for the construction of prognostic and diagnostic
    models for right-censored data. `r pkg("xgboost")` includes methods
    for Cox regression (right censored survival data) and AFT models
    (right-, left-, interval-, and uncensored).
-   ***Other:*** The `r pkg("superpc")` package implements
    the supervised principal components for survival data. The
    `r pkg("compound.Cox")` package fits Cox proportional
    hazards model using the compound covariate method.
    `r pkg("plsRcox")` provides partial least squares
    regression and various techniques for fitting Cox models in high
    dimensionnal settings.

Predictions and Prediction Performance
======================================

-   The `r pkg("pec")` package provides utilities to plot
    prediction error curves for several survival models. The
    `r pkg("riskRegression")` package now includes most of
    the functionality of the `r pkg("pec")` package.
-   `r pkg("peperr")` implements prediction error techniques
    which can be computed in a parallelised way. Useful for
    high-dimensional data.
-   The `r pkg("timeROC")` package permits to estimate
    time-dependent ROC curves and time-dependent AUC with censored data,
    possibly with competing risks.
-   `r pkg("survivalROC")` computes time-dependent ROC
    curves and time-dependent AUC from censored data using Kaplan-Meier
    or Akritas's nearest neighbour estimation method (Cumulative
    sensitivity and dynamic specificity).
-   `r pkg("tdROC")` can be used to compute time-dependent
    ROC curve from censored survival data using nonparametric weight
    adjustments.
-   `r pkg("risksetROC")` implements time-dependent ROC
    curves, AUC and integrated AUC of Heagerty and Zheng (Biometrics,
    2005).
-   Various time-dependent true/false positive rates and
    Cumulative/Dynamic AUC are implemented in the
    `r pkg("survAUC")` package.
-   The `r bioc("survcomp")` package provides several
    functions to assess and compare the performance of survival models.
-   C-statistics for risk prediction models with censored survival data
    can be computed via the `r pkg("survC1")` package.
-   The `r pkg("survIDINRI")` package implements the
    integrated discrimination improvement index and the category-less
    net reclassification index for comparing competing risks prediction
    models.
-   The `r pkg("compareC")` package permits to compare C
    indices with right-censored survival outcomes
-   The `r pkg("APtools")` package provide tools to estimate
    the average positive predictive values and the AUC for risk scores
    or marker.

Power Analysis
==============

-   The `r pkg("NPHMC")` permits to calculate sample size
    based on proportional hazards mixture cure models.
-   The `r pkg("powerSurvEpi")` package provides power and
    sample size calculation for survival analysis (with a focus towards
    epidemiological studies).
-   Power analysis and sample size calculation for SNP association
    studies with time-to-event outcomes can be done using the
    `r pkg("survSNP")` package.

Simulation
==========

-   The `r pkg("genSurv")` package permits to generate data
    wih one binary time-dependent covariate and data stemming from a
    progressive illness-death model.
-   The `r pkg("PermAlgo")` package permits the user to
    simulate complex survival data, in which event and censoring times
    could be conditional on an user-specified list of (possibly
    time-dependent) covariates.
-   The `r pkg("prodlim")` package proposes some functions
    for simulating complex event history data.
-   The `r pkg("gems")` package also permits to simulate and
    analyse multistate models. The package allows for a general
    specification of the transition hazard functions, for non-Markov
    models and for dependencies on the history.
-   The `r pkg("simMSM")` package provides functions for
    simulating complex multistate models data with possibly nonlinear
    baseline hazards and nonlinear covariate effects.
-   The `r pkg("simPH")` package implements tools for
    simulating and plotting quantities of interest estimated from
    proportional hazards models.
-   The `r pkg("survsim")` package permits to simulate
    simple and complex survival data such as recurrent event data and
    competing risks.
-   The `r pkg("simsurv")` package enables the user to
    simulate survival times from standard parametric survival
    distributions (exponential, Weibull, Gompertz), 2-component mixture
    distributions, or a user-defined hazard or log hazard function. Time
    dependent effects (i.e. non-proportional hazards) can be included by
    interacting covariates with linear time or some transformation of
    time.
-   The `r pkg("MicSim")` package provides routines for
    performing continuous-time microsimulation for population
    projection. The basis for the microsimulation are a multistate
    model, Markov or non-Markov, for which the transition intensities
    are specified, as well as an initial cohort.
-   The `r pkg("SimHaz")` package permits to simulate data
    with a dichotomous time-dependent exposure.
-   The `r pkg("SimSCRPiecewise")` package can be used to
    simulate univariate and semi-competing risks data given covariates
    and piecewise exponential baseline hazards.
-   The `r pkg("SimSurvNMarker")` package provides functions
    to simulate from joint survival and potentially multivariate marker
    models. User-defined basis expansions in time can be passed which
    effect the log hazard, the markers, and the association between the
    two.

Graphics
========

This section tries to list some specialised plot functions that might be
useful in the context of event history analysis.

-   The `r pkg("rms")` package proposes functions for
    plotting survival curves with the at risk table aligned to the x
    axis. `r pkg("prodlim")` extends this to the competing
    risks model.
-   The `plot.Hist` function in `r pkg("prodlim")` permits
    to draw the states and transitions that characterize a multistate
    model.
-   The `r pkg("Epi")` package provides many plot functions
    for representing multistate data, in particular Lexis diagrams.
-   The `r pkg("FamEvent")` generates time-to-event outcomes
    for families that habour genetic mutation under different sampling
    designs and estimates the penetrance functions for family data with
    ascertainment correction.
-   `r pkg("vsd")` provides graphical shim for visual survival data analysis. 

Miscellaneous
=============

-   The `r pkg("survminer")` package contains the function
    `ggsurvplot` for drawing survival curves with the "number at risk"
    table. Other functions are also available for visual examinations of
    Cox model assumptions.
-   The `r pkg("InformativeCensoring")` package multiple
    imputation methods for dealing with informative censoring.
-   The `r pkg("discSurv")` provides data transformations,
    estimation utilities, predictive evaluation measures and simulation
    functions for discrete time survival analysis.
-   `r pkg("dynpred")` is the companion package to "Dynamic
    Prediction in Clinical Survival Analysis".
-   Package `r pkg("boot")` proposes the `censboot` function
    that implements several types of bootstrap techniques for
    right-censored data.
-   The `r pkg("currentSurvival")` package estimates the
    current cumulative incidence and the current leukaemia free survival
    function.
-   The `r pkg("survJamda")` package provides functions for
    performing meta-analyses of gene expression data and to predict
    patients' survival and risk assessment.
-   The `r pkg("KMsurv")` package includes the data sets
    from Klein and Moeschberger (1997). The package
    `r pkg("SMPracticals")` that accompanies Davidson (2003)
    and `r pkg("DAAG")` that accompanies Maindonald, J.H.
    and Braun, W.J. (2003, 2007) also contain survival data sets.
-   The `r pkg("SvyNom")` package permits to construct,
    validate and calibrate nomograms stemming from complex
    right-censored survey data.
-   The `r pkg("logconcens")` package compute the MLE of a
    density (log-concave) possibly for interval censored data.
-   The `r pkg("coarseDataTools")` package implements an EM
    algorithm to estimate the relative case fatality ratio between two
    groups.
-   The `r pkg("GSSE")` package proposes a fully efficient
    sieve maximum likelihood method to estimate genotype-specific
    distribution of time-to-event outcomes under a nonparametric model
-   power and sample size calculation based on the difference in
    restricted mean survival times can be performed using the
    `r pkg("SSRMST")` package.
-   The `r pkg("survMisc")` provides miscellaneous routines
    to help in the analysis of right-censored survival data.
-   Accompanying data sets to the book *Applied Survival Analysis Using
    R* can be found in package `r pkg("asaur")`.



### Links
-   [GAMLSS](http://www.gamlss.com/)
-   [Tutorial in competing risks and multistate models](https://onlinelibrary.wiley.com/doi/10.1002/sim.2712)
-   [Proportional-Hazards Regression for Survival Data. Appendix to An R and S-PLUS Companion to Applied Regression.](https://socialsciences.mcmaster.ca/jfox/Books/Companion/appendices/Appendix-Cox-Regression.pdf)
-   [Journal of Statistical Software. Special Volume: Competing Risks and Multi-State Models](https://www.jstatsoft.org/v38)
