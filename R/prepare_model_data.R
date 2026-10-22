#' Prepare input data for the spatial models
#'
#' Takes raw sf objects for plots and grid, and prepares:
#' - landuse recoding
#' - flags (field visits, abiotics)
#' - z-scored abiotics
#' - terrain coding
#'
#' @param sf.plots sf object with plot data (including columns X, Y, bodemgebruik, Dagnummer, abiotics, etc.)
#' @param sf.grid sf object with grid data (including columns X, Y, bodemgebruik, etc.)
#' @param pH.mean numeric, mean pH used for z-scoring
#' @param pH.sd numeric, sd pH used for z-scoring
#'
#' @return list with elements `plots` and `grid` (prepared sf objects)
#' @export
prepare_model_data <- function(sf.plots,
                               sf.grid,
                               pH.mean = 7,
                               pH.sd = 1) {

  # Landuse recoding
  recode_landuse <- function(x) {
    dplyr::recode(x,
      "Bos" = "Forest",
      "Bebouwd" = "Built",
      "Droog natuurlijk terrein" = "Natural terrain (dry)",
      "Hoofdweg" = "Highway",
      "Landbouw en overig agrarisch" = "Agricultural",
      "Nat natuurlijk terrein" = "Natural terrain (wet)",
      "Recreatie" = "Recreational",
      "Spoorterrein" = "Railway",
      "Water" = "Water"
    )
  }

  sf.grid$landuse <- recode_landuse(sf.grid$bodemgebruik)
  sf.plots$landuse <- recode_landuse(sf.plots$bodemgebruik)

  sf.grid$landuse <- forcats::fct_relevel(
    sf.grid$landuse,
    "Forest", "Natural terrain (dry)", "Natural terrain (wet)",
    "Agricultural", "Recreational", "Highway", "Railway", "Water", "Built"
  )
  sf.plots$landuse <- forcats::fct_relevel(
    sf.plots$landuse,
    levels(sf.grid$landuse)
  )

  # Flags
  sf.plots$is_fieldvisit <- !is.na(sf.plots$Dagnummer)
  sf.plots$has_abiotics  <- !is.na(sf.plots$releve_nr)

  # Z-scoring helpers
  to_zscore <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  zclip <- function(x) ifelse(x > 3, 3, ifelse(x < -3, -3, x))

  n <- sum(sf.plots$is_fieldvisit)

  sf.plots <- sf.plots %>%
    dplyr::mutate(
      Pq_grootte = tidyr::replace_na(Pq_grootte, 25),
      Dagnummer = tidyr::replace_na(Dagnummer, 200),
      Dagnummer.z = (Dagnummer - 200) / 30.44,
      pH.z = zclip((pH - pH.mean) / pH.sd),
      ORG_STOF.z = zclip(to_zscore(log(ORG_STOF))),
      CN.z = zclip(to_zscore(log(`C/N`))),
      N_totaal.z = zclip(to_zscore(log(N_totaal))),
      P_totaal.z = zclip(to_zscore(log(P_totaal))),
      K_CaCl2.z = zclip(to_zscore(log(ifelse(K_CaCl2 <= 0, 0.1, K_CaCl2))))
    )

  sf.grid <- sf.grid %>%
    dplyr::mutate(
      Dagnummer = 200,
      Pq_grootte = 25,
      Dagnummer.z = (Dagnummer - 200) / 30.44,
      pH.z = zclip((pH - pH.mean) / pH.sd)
    )

  # Terrain coding
  code_terrain <- function(x) {
    switch(x,
      "Bos" = -4,
      "Droog natuurlijk terrein" = -3,
      "Nat natuurlijk terrein" = -2,
      "Water" = -1,
      "Recreatie" = 0,
      "Spoorterrein" = 1,
      "Hoofdweg" = 2,
      "Bebouwd" = 3,
      "Landbouw en overig agrarisch" = 4
    )
  }

  sf.plots$terrain <- sapply(as.character(sf.plots$bodemgebruik), code_terrain)
  sf.grid$terrain  <- sapply(as.character(sf.grid$bodemgebruik), code_terrain)

  # Coordinates matrix

  xy.plots <- sf.plots %>% as.data.frame() %>% dplyr::select(X, Y)
  xy.grid  <- sf.grid %>% as.data.frame() %>% dplyr::select(X, Y)

  list(
    plots = sf.plots,
    grid  = sf.grid,
    xy.plots = xy.plots,
    xy.grid  = xy.grid
  )
}