# ============================================================
# Trivariate phylogenetic mixed model: fliC, T6SS, bacteriocins
# ============================================================

library(MCMCglmm)
library(ape)
library(coda)

# --- 1. Load data and phylogeny ---
dat <- read.csv("Merged_counts.csv")        # columns: species_name, fliC, T6SS_count, bact_count
phy <- read.nexus("4125_pruned_tree1.nexus")  # phylo object

# --- 1a. Randomly select 500 species if >500 ---
set.seed(123)
if(nrow(dat) > 500){
  subdata <- dat[sample(nrow(dat), 500), ]
} else {
  subdata <- dat
}

# --- 1b. Prune tree ---
subspecies <- subdata$species_name
phy <- drop.tip(phy, setdiff(phy$tip.label, subspecies))
subdata$species_name <- factor(subdata$species_name, levels = phy$tip.label)

cat("Number of species for MCMCglmm:", nrow(subdata), "\n")

# --- 2. Prepare responses ---

# fliC: threshold (binary)
subdata$fliC_num <- as.numeric(as.character(subdata$fliC))  # 0/1

# T6SS: multinomial counts (raw integers)
subdata$T6SS_success <- subdata$T6SS_count
subdata$T6SS_fail    <- 11 - subdata$T6SS_count   # 11, Max no of T6SS core proteins

# Bacteriocins: Poisson counts (raw integers)
subdata$bact_count <- subdata$bact_count

# --- 3. Inverse phylogenetic matrix ---
Ainv <- inverseA(phy, nodes = "TIPS", scale = TRUE)$Ainv

# --- 4. Define priors ---
prior_tri <- list(
  G = list(
    G1 = list(
      V = diag(3),
      nu = 4,
      alpha.mu = rep(0,3),
      alpha.V = diag(3)*1000
    )
  ),
  R = list(
    V = diag(3),
    nu = 4,
    fix = 2   # fixes residual variance for threshold trait (fliC)
  )
)

# --- 5. Fit trivariate model ---
m_tri <- MCMCglmm(
  cbind(fliC_num, cbind(T6SS_success, T6SS_fail), bact_count) ~ trait - 1,
  random = ~ us(trait):species_name,
  rcov = ~ us(trait):units,
  family = c("threshold", "multinomial2", "poisson"),
  ginverse = list(species_name = Ainv),
  data = subdata,
  prior = prior_tri,
  nitt   = 500000,
  burnin = 100000,
  thin   = 500,
  pr = TRUE,
  verbose = TRUE
)

# --- 6. Check convergence ---
effectiveSize(m_tri$VCV)
effectiveSize(m_tri$Sol)
autocorr.diag(m_tri$VCV)
autocorr.diag(m_tri$Sol)
plot(m_tri$VCV)
plot(m_tri$Sol)

# --- 7. Extract phylogenetic correlations ---
VCV <- m_tri$VCV

v_fli_phy  <- VCV[,"traitfliC_num:traitfliC_num.species_name"]
v_T6_phy   <- VCV[,"traitT6SS_success:traitT6SS_success.species_name"]
v_bact_phy <- VCV[,"traitbact_count:traitbact_count.species_name"]

cov_fli_T6  <- VCV[,"traitfliC_num:traitT6SS_success.species_name"]
cov_fli_bact <- VCV[,"traitfliC_num:traitbact_count.species_name"]
cov_T6_bact <- VCV[,"traitT6SS_success:traitbact_count.species_name"]

rho_fli_T6_phy <- cov_fli_T6 / sqrt(v_fli_phy * v_T6_phy)
rho_fli_bact_phy <- cov_fli_bact / sqrt(v_fli_phy * v_bact_phy)
rho_T6_bact_phy <- cov_T6_bact / sqrt(v_T6_phy * v_bact_phy)

cat("\n--- Phylogenetic correlations (median + 95% CI) ---\n")
cat("fliC - T6SS:\n"); print(quantile(rho_fli_T6_phy, c(0.025,0.5,0.975)))
cat("fliC - bacteriocins:\n"); print(quantile(rho_fli_bact_phy, c(0.025,0.5,0.975)))
cat("T6SS - bacteriocins:\n"); print(quantile(rho_T6_bact_phy, c(0.025,0.5,0.975)))

# --- 8. Phylogenetic signal (lambda) ---
resid_var_fli <- VCV[,"traitfliC_num:traitfliC_num.units"]
resid_var_T6 <- VCV[,"traitT6SS_success:traitT6SS_success.units"]
resid_var_bact <- VCV[,"traitbact_count:traitbact_count.units"]

lambda_fli <- v_fli_phy / (v_fli_phy + resid_var_fli)
lambda_T6 <- v_T6_phy / (v_T6_phy + resid_var_T6)
lambda_bact <- v_bact_phy / (v_bact_phy + resid_var_bact)

cat("\n--- Phylogenetic signal (lambda) ---\n")
cat("fliC lambda:\n"); print(quantile(lambda_fli, c(0.025,0.5,0.975)))
cat("T6SS lambda:\n"); print(quantile(lambda_T6, c(0.025,0.5,0.975)))
cat("Bacteriocins lambda:\n"); print(quantile(lambda_bact, c(0.025,0.5,0.975)))

# --- 9. Save results ---
save(m_tri, file = "trivariate_mcmcglmm_output_clean.RData")
