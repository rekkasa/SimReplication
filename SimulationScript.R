# ==============================================================================
# Required packages
# ==============================================================================
# install.packages("remotes")
#
# remotes::install_github("rekkasa/SimulateHte")
# remotes::install_github("rekkasa/SmoothHte")
# remotes::install_github("rekkasa/SimulationEvaluationHte")
# remotes::install_github("ohdsi/ParallelLogger)

library(dplyr)

# Contains all settings for the main simulations
analysisIds <- readr::read_csv(
  "https://raw.githubusercontent.com/mi-erasmusmc/HteSimulationRCT/main/data/processed/analysisIds.csv",
  col_types =  readr::cols(
    .default = readr::col_double(),
    base = readr::col_character(),
    type = readr::col_character(),
    harm = readr::col_character()
  )
)

selectedScenario <- 1    # This can be any of the possible 648 scenarios
idSettings <- analysisIds %>%
  dplyr::filter(scenario == selectedScenario)


# Settings for constant treatment-related harm
if (idSettings$base != "absent") {
  harm <- dplyr::case_when(
    idSettings$harm == "moderate-positive" ~ idSettings$averageTrueBenefit / 4,
    idSettings$harm == "strong-positive" ~ idSettings$averageTrueBenefit / 2,
    idSettings$harm == "negative" ~ -idSettings$averageTrueBenefit / 4,
    TRUE ~ 0
  )
} else {
  harm <- dplyr::case_when(
    idSettings$harm == "moderate-positive" ~ .01,
    idSettings$harm == "strong-positive" ~ .02,
    idSettings$harm == "negative" ~ -.01,
    TRUE ~ 0
  )
}

# Functions for generating linear, quadratic, and non-monotonic deviations
# from constant relative treatment effect
createF1 <- function(c) function(x) x - c
createF2 <- function(c) function(x) (x - c)^2


# createModelSettings defines the true treatment effect function. The settings
# below define the function:
#        lp1 = g0 + g1 * (lp0 - c) + g2 * (lp0 - c)^2
# where lp1 is the true linear predictor in the treatment arm and
# lp0 is the true linear predictor in the control arm (see paper)
treatmentEffectSettings <- SimulateHte::createTreatmentEffectSettings(
  type = "lp",
  harm = harm,
  modelSettings = SimulateHte::createModelSettings(
    constant = idSettings$g0,
    modelMatrix = matrix(c(1, 1)),
    transformationSettings = list(
      createF1(idSettings$c),
      createF2(idSettings$c)
    ),
    coefficients = c(
      idSettings$g1,
      idSettings$g2
    )
  )
)

# Distribution settings for baseline covariates
databaseSettings <- SimulateHte::createDatabaseSettings(
  numberOfObservations = as.numeric(as.character(idSettings$sampleSize)),
  numberOfCovariates = 8,
  covariateDistributionSettings = list(
    SimulateHte::createNormalDistributionSettings(),
    SimulateHte::createNormalDistributionSettings(),
    SimulateHte::createNormalDistributionSettings(),
    SimulateHte::createNormalDistributionSettings(),
    SimulateHte::createBinomialDistributionSettings(prob = .2),
    SimulateHte::createBinomialDistributionSettings(prob = .2),
    SimulateHte::createBinomialDistributionSettings(prob = .2),
    SimulateHte::createBinomialDistributionSettings(prob = .2)
  )
)


# Baseline risk settings (lp0). Settings here define:
#         lp0 = b0 + b1 * x1 + ... + b8 * x8
# p0 = exp(lp0) / (1 + exp(lp0)) -> baseline risk under control
baselineRiskSettings <- SimulateHte::createBaselineRiskSettings(
  type = "binary",
  modelSettings = SimulateHte::createModelSettings(
    constant = idSettings %>% pull(b0),
    modelMatrix = diag(8),
    transformationSettings = list(
      identity,
      identity,
      identity,
      identity,
      identity,
      identity,
      identity,
      identity
    ),
    coefficients = idSettings %>% select(paste0("b", 1:8)) %>% unlist()
  )
)

# Treatment probability (here it's set to 0.5)
propensitySettings <- SimulateHte::createPropensitySettings(
  type = "binary",
  modelSettings = SimulateHte::createModelSettings(
    constant = 0,
    modelMatrix = diag(0),
    transformationSettings = NULL
  )
)

# Combine all settings
simulationSettings <- list(
  databaseSettings = databaseSettings,
  propensitySettings = propensitySettings,
  baselineRiskSettings = baselineRiskSettings,
  treatmentEffectSettings = treatmentEffectSettings
)

analysisSettings <- SimulationEvaluationHte::createAnalysisSettings(
  threads = 1, # define number of threads for simulations (1 run per thread)
  seed = 19910930,
  replications = 500,
  validationSize = 5e5,
  analysisId = paste(
    "scenario",
    idSettings$scenario,
    sep = "_"
  ),
  description = "description",
  saveDirectory = "data/raw"
)


# Definition of methods used in the paper
smoothSettings <- list(
  constant = SmoothHte::createHteSettings(
    label = "constant_treatment_effect",
    settings = SmoothHte::createModelBasedSettings(
      type = "treatment",
      model = "logistic"
    )
  ),
  stratified = SmoothHte::createHteSettings(
    settings = SmoothHte::createStratifiedSettings(),
    label = "stratified"
  ),
  constantLp = SmoothHte::createHteSettings(
    label = "linear_predictor",
    settings = SmoothHte::createModelBasedSettings()
  ),
  rcs3 = SmoothHte::createHteSettings(
    settings = SmoothHte::createRcsSettings(),
    label = "rcs_3_knots"
  ),
  rcs4 = SmoothHte::createHteSettings(
    settings = SmoothHte::createRcsSettings(nKnots = 4),
    label = "rcs_4_knots"
  ),
  rcs5 = SmoothHte::createHteSettings(
    settings = SmoothHte::createRcsSettings(nKnots = 5),
    label = "rcs_5_knots"
  ),
  adaptive = SmoothHte::createHteSettings(
    settings = SmoothHte::createAdaptiveSettings(
      list(
        rcs3 = SmoothHte::createRcsSettings(),
        rcs4 = SmoothHte::createRcsSettings(nKnots = 4),
        rcs5 = SmoothHte::createRcsSettings(nKnots = 5),
        linearSettings = SmoothHte::createModelBasedSettings(),
        constantSettings = SmoothHte::createModelBasedSettings(
          type = "treatment",
          model = "logistic"
        )
      )
    ),
    label = "adaptive"
  )
)

predictionSettings <- SimulationEvaluationHte::createPredictionSettings(
  args = list(
    formula = "outcome ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + treatment",
    family = "binomial"
  ),
  fun = "glm"
)

res <- SimulationEvaluationHte::runAnalysis(
  analysisSettings = analysisSettings,
  simulationSettings = simulationSettings,
  predictionSettings = predictionSettings,
  smoothSettings = smoothSettings
)
