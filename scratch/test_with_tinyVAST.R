


#devtools::install_github("james-thorson/lanczosRTMB@dev", force=TRUE)
#setwd( R'(C:\Users\James.Thorson\Desktop\Work files\AFSC\2025-05 -- Rasumusson ODFG)' )
setwd( R'(C:\Users\jtuth\Desktop\Work\AFSC\2025-05 -- Rasumusson ODFG)' )

library(sdmTMB)
library(sdmTMBextra)
library(sf)
library(RTMB)
library(Matrix)
library(lanczosRTMB)

load("ODFW_RockfishData.RData")
domain = st_read( file.path("SA_0to80m", "SA_0to80m.shp") )
domain = st_geometry(domain)

# Add UTM
get_crs(dat1,c("Lon_M","Lat_M"))
dat1 <- add_utm_columns(dat1,c("Lon_M","Lat_M"))


library(tinyVAST)

# Choose grid
sf_grid = readRDS( file.path(R'(C:\Users\jtuth\Desktop\Work\AFSC\2026-01 -- Rasmusson ODFG)', "sf_grid_250m.RDS") )
sf_grid = st_transform( sf_grid, crs = st_crs('+proj=utm +zone=10 +units=km') )

# Rescale covariates for numerical stability
df_data = as.data.frame(dat1)
sf_data = st_as_sf( df_data[,c("X","Y")], coords = c("X","Y"), crs = st_crs(sf_grid) )
within = st_within( sf_data, sf_grid )
df_data = df_data[ which(!is.na(as.integer(within))), ]
df_data$Depth = (df_data$Depth - 40) / 10
df_data$HR = (df_data$HR - 10) / 10

nngp_domain = make_nngp_domain(
  sf_areal = sf_grid,
  nn = 4
)

tv_nngp <- tinyVAST(
  formula = density_m2 ~ 0 + Substrate2 + O2 + s(Depth), #+ s(Depth,k=5) + Substrate2 + O2 + 0 ,  # s(HR,k=5)  had log_lambda -> Inf
  data = df_data,
  spatial_domain = nngp_domain,
  family = tweedie("log"),
  space_columns = c("X","Y"),
  time_column = "Year",
  space_term = "",
  control = tinyVASTcontrol(
    #iter.max = 3,
    #eval.max = 3,
    gmrf_parameterization = "projection",
    getsd = FALSE
  )
)
tv_nngp

########################
# Try in lanczos_MakeADFun
########################

source( file.path( R'(C:\Users\jtuth\Documents\GitHub\lanczosRTMB\R)', "lanczos.R") )

lanobj = lanczos_MakeADFun(
  obj = tv_nngp$obj,
  parameters = tv_nngp$tmb_inputs$tmb_par,
  random = tv_nngp$tmb_inputs$tmb_random,
  k = 30,
  make_gr = TRUE
)

#
opt = nlminb(
  lanobj$par,
  lanobj$fn,
  lanobj$gr,
  control = list( trace = 1, iter.max = 1e5, eval.max = 1e5 )
)

