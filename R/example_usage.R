library(sf)
library(dplyr)

set.seed(123)

# 1.1 Simuleer plot-coördinaten
n_plots <- 200
plots_df <- data.frame(
  X = runif(n_plots, 0, 10000),
  Y = runif(n_plots, 0, 10000)
)

# 1.2 Simuleer grid-coördinaten (bijv. raster-achtige punten)
n_grid <- 1000
grid_df <- data.frame(
  X = runif(n_grid, 0, 10000),
  Y = runif(n_grid, 0, 10000)
)

# 1.3 Landuse-categorieën (gebruik jouw Nederland-lijst als basis)
landuse_levels_nl <- c(
  "Bos",
  "Droog natuurlijk terrein",
  "Nat natuurlijk terrein",
  "Water",
  "Recreatie",
  "Spoorterrein",
  "Hoofdweg",
  "Bebouwd",
  "Landbouw en overig agrarisch"
)

plots_df$bodemgebruik <- sample(landuse_levels_nl, n_plots, replace = TRUE)
grid_df$bodemgebruik  <- sample(landuse_levels_nl, n_grid, replace = TRUE)

# 1.4 Maak sf-objecten (plots en grid)
sf.plots <- st_as_sf(plots_df, coords = c("X", "Y"), crs = 28992) # bv. RD-coördinaten
sf.grid  <- st_as_sf(grid_df,  coords = c("X", "Y"), crs = 28992)

# Voor gemak: X/Y weer als kolommen bewaren
sf.plots$X <- plots_df$X
sf.plots$Y <- plots_df$Y
sf.grid$X  <- grid_df$X
sf.grid$Y  <- grid_df$Y

# 1.5 Simuleer pH als ruimtelijk patroon + landuse-effect + ruis
# Maak een simpele "centrum"-afstand
center_x <- 5000
center_y <- 5000

dist_plots <- sqrt((sf.plots$X - center_x)^2 + (sf.plots$Y - center_y)^2)
dist_grid  <- sqrt((sf.grid$X  - center_x)^2 + (sf.grid$Y  - center_y)^2)

# Baseline pH afhankelijk van afstand tot centrum (lager als verder weg)
sf.plots$pH <- 7 - dist_plots / 10000 + rnorm(n_plots, sd = 0.2)
sf.grid$pH  <- 7 - dist_grid  / 10000 + rnorm(n_grid,  sd = 0.2)

# 1.6 Simuleer een soort (presence/absence) die van pH en landuse afhangt
# Zet eerst een numerieke terrain-code (zoals in je script)
code_terrain <- function(x) {
  switch(x,
    "Bos"                        = -4,
    "Droog natuurlijk terrein"   = -3,
    "Nat natuurlijk terrein"     = -2,
    "Water"                      = -1,
    "Recreatie"                  = 0,
    "Spoorterrein"               = 1,
    "Hoofdweg"                   = 2,
    "Bebouwd"                    = 3,
    "Landbouw en overig agrarisch" = 4
  )
}

sf.plots$terrain <- sapply(as.character(sf.plots$bodemgebruik), code_terrain)
sf.grid$terrain  <- sapply(as.character(sf.grid$bodemgebruik),  code_terrain)

# Logit-lijn: soort komt meer voor bij hogere pH, in bos/natuurlijk terrein
logit <- -1 + 1.5 * (sf.plots$pH - 7) - 0.2 * sf.plots$terrain
prob  <- 1 / (1 + exp(-logit))
sf.plots$species_presence <- rbinom(n_plots, size = 1, prob = prob)

# Voeg wat andere kolommen toe die je data-prep functie verwacht
sf.plots$Pq_grootte  <- sample(c(25, 50, NA), n_plots, replace = TRUE)
sf.plots$Dagnummer   <- sample(c(150:250, NA), n_plots, replace = TRUE)
sf.plots$releve_nr   <- sample(c(1:5, NA), n_plots, replace = TRUE)
sf.plots$ORG_STOF    <- exp(rnorm(n_plots, 2, 0.4))
sf.plots$`C/N`       <- exp(rnorm(n_plots, 1, 0.3))
sf.plots$N_totaal    <- exp(rnorm(n_plots, -1, 0.5))
sf.plots$P_totaal    <- exp(rnorm(n_plots, -0.5, 0.5))
sf.plots$K_CaCl2     <- exp(rnorm(n_plots, -2, 0.5))

# Voor grid hoeven we alleen pH, bodemgebruik, X, Y en evt. Dagnummer/Pq_grootte
sf.grid$Dagnummer  <- 200
sf.grid$Pq_grootte <- 25

#2. Data voorbereiden met prepare_model_data()

prep <- prepare_model_data(
  sf.plots = sf.plots,
  sf.grid  = sf.grid,
  pH.mean  = 7,
  pH.sd    = 1
)

str(prep$plots)
str(prep$grid)
head(prep$xy.plots)
head(prep$xy.grid)

# 3.1 Bouw de mesh
mesh <- build_mesh(
  xy.plots = prep$xy.plots,
  xy.grid  = prep$xy.grid,
  mindist  = 500  # iets ruimer, want het domein is 0-10000
)

plot(mesh)   # snelle check: mesh ziet er logisch uit?

# 3.2 Bouw SPDE-model en projector-matrices
spde_obj <- build_spde_and_projectors(
  mesh     = mesh,
  xy.plots = prep$xy.plots,
  xy.grid  = prep$xy.grid,
  alpha    = 2
)

str(spde_obj$spde)
dim(spde_obj$A.plots)
dim(spde_obj$A.grid)

res_two_stage <- fit_two_stage_sdm(
  prep      = prep,
  spde_obj  = spde_obj,
  abiotic_response  = "pH.z",
  species_response  = "species_presence",
  abiotic_formula   = ~ landuse + terrain + Dagnummer.z,
  species_formula   = ~ landuse + terrain + Dagnummer.z,
  family_stage1     = "gaussian",
  family_stage2     = "binomial",
  control.compute   = list(dic = TRUE, waic = TRUE)
)

# Resultaten inspecteren
names(res_two_stage)

# Fixed effects in abiotic model (stage 1)
res_two_stage$fit_stage1$summary.fixed

# Fixed effects in species model (stage 2)
res_two_stage$fit_stage2$summary.fixed

# Waic/dic
res_two_stage$fit_stage1$waic
res_two_stage$fit_stage2$waic

abiotic_grid_pred <- res_two_stage$abiotic_grid_pred
prep$grid$abiotic_pred <- abiotic_grid_pred

# Snel plotten met ggplot2 als check
library(ggplot2)

ggplot(prep$grid) +
  geom_point(aes(X, Y, color = abiotic_pred), size = 1) +
  scale_color_viridis_c() +
  coord_equal() +
  labs(title = "Voorspelde pH (z-score) op grid uit stage 1")