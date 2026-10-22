#' Fit two-stage spatial SDM
#'
#' Stage 1: abiotic model (e.g. pH ~ landuse + terrain + spatial field)  
#' Stage 2: species distribution model using abiotic predictions.
#'
#' Deze functie is een sjabloon gebaseerd op jouw script; pas de formule(s)
#' en response aan jouw case aan.
#'
#' @param prep list zoals geretourneerd door \code{prepare_model_data()},
#'   met minimaal elementen `plots`, `grid`, `xy.plots`, `xy.grid`.
#' @param spde_obj list zoals geretourneerd door \code{build_spde_and_projectors()},
#'   met elementen `spde`, `A.plots`, `A.grid`.
#' @param abiotic_response character, naam van de kolom in \code{prep$plots}
#'   die als response in het eerste (abiotics) model dient, bv. `"pH.z"`.
#' @param species_response character, naam van de kolom in \code{prep$plots}
#'   met de soortrespons (aanwezigheid/abundantie), bv. `"species_presence"`.
#' @param abiotic_formula formula voor stage 1, zonder response.
#'   Default: \code{~ landuse + terrain + Dagnummer.z}.
#' @param species_formula formula voor stage 2, zonder abiotic term; de
#'   abiotic voorspelling wordt automatisch toegevoegd.
#'   Default: \code{~ landuse + terrain + Dagnummer.z}.
#' @param family_stage1 INLA family voor abiotic model (default: "gaussian").
#' @param family_stage2 INLA family voor species model (default: "binomial").
#' @param control.compute list doorgegeven aan \code{INLA::inla()}, bv.
#'   \code{list(dic = TRUE, waic = TRUE, cpo = TRUE)}.
#'
#' @return list met:
#'   \item{fit_stage1}{INLA object voor abiotic model}
#'   \item{fit_stage2}{INLA object voor species model}
#'   \item{abiotic_grid_pred}{vector of matrix met voorspelde abiotics op grid}
#' @export
fit_two_stage_sdm <- function(prep,
                              spde_obj,
                              abiotic_response = "pH.z",
                              species_response = "species_presence",
                              abiotic_formula = ~ landuse + terrain + Dagnummer.z,
                              species_formula = ~ landuse + terrain + Dagnummer.z,
                              family_stage1 = "gaussian",
                              family_stage2 = "binomial",
                              control.compute = list(dic = TRUE,
                                                     waic = TRUE,
                                                     cpo = TRUE)) {

  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("Package 'INLA' is required but not installed.")
  }

  sf.plots <- prep$plots
  sf.grid  <- prep$grid

  # --- STAGE 1: ABIOTICS MODEL ---------------------------------------------

  # Maak indices voor het ruimtelijke veld
  n.spde <- spde_obj$spde$n.spde
  field_index <- INLA::inla.spde.make.index("field", n.spde)

  # Projectormatrices
  A.plots <- spde_obj$A.plots
  A.grid  <- spde_obj$A.grid

  # Stack voor plots (observed abiotics)
  formula_stage1 <- stats::as.formula(
    paste(
      abiotic_response,
      "~",
      deparse(abiotic_formula),
      "+ f(field, model = spde_obj$spde)"
    )
  )

  stack_stage1 <- INLA::inla.stack(
    data = list(
      y = sf.plots[[abiotic_response]]
    ),
    A = list(A.plots, 1),
    effects = list(
      field_index,
      sf.plots
    ),
    tag = "stage1"
  )

  fit_stage1 <- INLA::inla(
    formula_stage1,
    data = INLA::inla.stack.data(stack_stage1),
    family = family_stage1,
    control.compute = control.compute
  )

  # Voorspellingen op grid (abiotics)
  # We maken een nieuwe stack met dezelfde latent field + grid covariates
  stack_stage1_grid <- INLA::inla.stack(
    data = list(
      y = NA
    ),
    A = list(A.grid, 1),
    effects = list(
      field_index,
      sf.grid
    ),
    tag = "stage1_grid"
  )

  # Combineer stacks (observed + prediction)
  stack_stage1_combined <- INLA::inla.stack(stack_stage1, stack_stage1_grid)

  fit_stage1_pred <- INLA::inla(
    formula_stage1,
    data = INLA::inla.stack.data(stack_stage1_combined),
    family = family_stage1,
    control.predictor = list(compute = TRUE),
    control.compute = control.compute
  )

  # Extract voorspellingen voor grid
  idx_grid <- INLA::inla.stack.index(stack_stage1_combined, "stage1_grid")$data
  abiotic_grid_pred <- fit_stage1_pred$summary.linear.predictor[idx_grid, "mean"]

  # Voeg voorspelde abiotic waarde toe aan grid en (optioneel) plots
  sf.grid$abiotic_pred <- abiotic_grid_pred

  # Voor stage 2 kun je kiezen:
  # - gebruik observed abiotics in plots
  # - of voorspelde waarden op plots (hier houden we observed voor plots)
  sf.plots$abiotic_pred <- sf.plots[[abiotic_response]]

  # --- STAGE 2: SPECIES DISTRIBUTION MODEL ---------------------------------

  formula_stage2 <- stats::as.formula(
    paste(
      species_response,
      "~ abiotic_pred +",
      deparse(species_formula),
      "+ f(field, model = spde_obj$spde)"
    )
  )

  # Stack voor species (alleen op plots, want daar heb je soortdata)
  stack_stage2 <- INLA::inla.stack(
    data = list(
      y = sf.plots[[species_response]]
    ),
    A = list(A.plots, 1),
    effects = list(
      field_index,
      sf.plots
    ),
    tag = "stage2"
  )

  fit_stage2 <- INLA::inla(
    formula_stage2,
    data = INLA::inla.stack.data(stack_stage2),
    family = family_stage2,
    control.compute = control.compute
  )

  list(
    fit_stage1       = fit_stage1,
    fit_stage2       = fit_stage2,
    abiotic_grid_pred = abiotic_grid_pred
  )
}