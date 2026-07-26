
#pak::pak("james-thorson-noaa/ecostate")

library(ecostate)

# Load, inspect, and attach data
data(gulf_of_alaska)
attach(gulf_of_alaska)

# Biomass-dynamics steps
n_step = 50

# Age-structured dynamics steps
STEPS_PER_YEAR = 24

# Constant expected recruitment (matching assessment models)
# So h = 0.999 and back-calculate SpawnX
SpawnX = c( "Walleye pollock" = 4 / (5-1/0.999), "Sablefish" = 4 / (5-1/0.999) )

# Default vulnerability as fixed or starting values
X = array( 2,
           dim = rep(length(taxa),2),
           dimnames = list(Prey=taxa,Predator=taxa) )

# Unassimilated food
U = array( 0.2, dim=length(taxa), dimnames=list(taxa) )

# Estimate catchability coefficient for all surveys
fit_Q = c( "Walleye pollock adult", "Sablefish adult",
           "Euphausiids", "Large copepods" )

# Fit biomass-dynamics process errors for zooplankton
fit_eps = c( "Euphausiids", "Large copepods" )   #

# Fit recruitment deviations for age-structured populations
fit_phi = c("Walleye pollock", "Sablefish")

# Fit equilibrium biomass for age-structured populations
# using fishery as depletion experiment to identify scale
fit_B = c("Walleye pollock juv", "Sablefish adult")

# Fit PB with a prior
fit_PB = c( "Walleye pollock adult", "Sablefish adult" )

# Fit vulnerability for adult age-structured populations with a prior
fit_X = c("Walleye pollock adult", "Sablefish adult")

log_prior = function(p){
  lpp = 0
  # Normal(0,0.5) prior on diag(Xprime_ij), where X = exp(Xprime) + 1
  lpp = lpp + sum(dnorm( diag(p$Xprime_ij), mean=0, sd=0.5, log=TRUE ))
  # Normal prior on log(q) for adult pollock, matching stock assessment
  lpp = lpp + dnorm( p$logq_i[match('Walleye pollock adult',taxa)],
                     mean=log(0.85), sd=0.1, log=TRUE )
  # Tight normal prior on log(PB) for sablefish, matching assessment M value
  lpp = lpp + dnorm( p$logPB_i['Sablefish adult'],
                     mean=log(0.1), sd=0.1, log=TRUE )
  # Tight normal prior on log(PB) for pollock, matching assessment M value
  lpp = lpp + dnorm( p$logPB_i['Walleye pollock adult'],
                    mean=log(0.3), sd=0.1, log=TRUE )
  return(lpp)
}

out = ecostate(
  taxa = taxa,
  years = years,
  catch = catch_data,
  biomass = biomass_data,
  agecomp = agecomp_data,
  PB = P_over_B,
  QB = Q_over_B,
  DC = Diet_proportions,
  B = B,
  EE = EE,
  X = X,
  type = type,
  U = U,
  fit_B = fit_B,
  fit_Q = fit_Q,
  fit_PB = fit_PB,
  fit_eps = fit_eps,
  log_prior = log_prior,
  control = ecostate_control(
    n_steps = n_step,
    profile = NULL, # Penalized likelihood so use empty set
    random = NULL, # Penalized likelihood so use empty set
    nlminb_loops = 0,
    getsd = FALSE ),
  settings = stanza_settings(
    taxa = taxa,
    stanza_groups = stanza_groups,
    K = K,
    Wmat = Wmat,
    Amat = Amat,
    d = d,
    Amax = Amax,
    SpawnX = SpawnX,
    fit_phi = fit_phi,
    STEPS_PER_YEAR = STEPS_PER_YEAR,
    comp_weight = "multinom",
    Leading = Leading,
    Wmatslope = Wmatslope)
)

map = out$tmb_inputs$map
tmb_par = out$obj$env$parList()

# Estimate vulnerability for adult fish with prior
map$Xprime_ij = factor( ifelse(taxa %in% fit_X, seq_along(taxa), NA)[col(array(dim=rep(length(taxa),2)))] )

# Penalized likelihood:  Fixed SD for process errors = 1
map$logtau_i = factor(rep(NA,length(map$logtau_i)))
tmb_par$logtau_i = ifelse( is.na(tmb_par$logtau_i), NA, log(1) )

# Penalized likelihood:  Fixed SD for recruitment deviations = 1
map$logpsi_g2 = factor(rep(NA,length(map$logpsi_g2)))
tmb_par$logpsi_g2 = ifelse( is.na(tmb_par$logpsi_g2), NA, log(1) )

out = ecostate(
  taxa = taxa,
  years = years,
  catch = catch_data,
  biomass = biomass_data,
  agecomp = agecomp_data,
  PB = P_over_B,
  QB = Q_over_B,
  DC = Diet_proportions,
  B = B,
  EE = EE,
  X = X,
  type = type,
  U = U,
  fit_B = fit_B,
  fit_Q = fit_Q,
  fit_PB = fit_PB,
  fit_eps = fit_eps,
  log_prior = log_prior,
  control = ecostate_control(
    n_steps = n_step,
    profile = NULL, # Penalized likelihood so use empty set
    random = NULL, # Penalized likelihood so use empty set
    derived_quantities = c(), # Turn off for speed
    map = map,       # Pass map back in
    tmb_par = tmb_par,  # Pass parameters back in
    getsd = FALSE ),
  settings = stanza_settings(
    taxa = taxa,
    stanza_groups = stanza_groups,
    K = K,
    Wmat = Wmat,
    Amat = Amat,
    d = d,
    Amax = Amax,
    SpawnX = SpawnX,
    fit_phi = fit_phi,
    STEPS_PER_YEAR = STEPS_PER_YEAR,
    comp_weight = "multinom",
    Leading = Leading,
    Wmatslope = Wmatslope)
)

###########################
# TRY WITH lanczosRTMB
###########################

library( lanczosRTMB )

tape = GetTape(out$obj)
grad = tape$jacfun()
Hq = make_Hq(
  tape = tape,
  method =  "FD-on-reverse",
  x0 = out$obj$par
)

start_time = Sys.time()
opt = newton_CG(
  par = out$obj$par,
  fn = tape,
  gr = grad,
  ustep = 0.001,
  #maxit_CG = 20,
  #e_ratio = 0.1,
  #power = 0.2,
  Hq = Hq
)
opt$run_time = Sys.time() - start_time

Hq = make_Hq(
  tape = tape,
  method =  "FD-on-reverse",
  x0 = out$opt$par
)
start_time = Sys.time()
logdet = lanczos_logdet(
  Hq = Hq,
  k = 40,
  m = 3
)
run_time = Sys.time() - start_time

lanobj = lanczos_MakeADFun(
  func = out$tmb_inputs$func,
  parameters = out$tmb_inputs$p,
  map = out$tmb_inputs$map,
  random = ecostate_control()$random,
  profile = ecostate_control()$profile
)

func = out$tmb_inputs$func
parameters = out$tmb_inputs$p
map = out$tmb_inputs$map
random = ecostate_control()$random
profile = ecostate_control()$profile
k = 40
m = 3
method = "newton_CG"
seed = 123
make_gr = TRUE
pu_update = c("implicit", "FD", "exact")[1]

silent = TRUE
#CG(
#  b = out$obj$par,
#  #x = 0*out$obj$par + 0.01,
#  Hq = Hq
#)
#
#library(Matrix)
#b = out$obj$par + 0.01
##Hq
#x = 0 * b
#Minv = Diagonal(n = length(b))
#max.it = length(b)
#e = 1e-10
#stop_if_nonPD = TRUE
#silent = TRUE

