## =====================================================================
##  DiabetesCtmcToOdeV31_Analysis.R
##  Companion analysis for manuscript "DiabetesCtmcToOdeV31.tex"
##
##  Equilibrium-constrained Maximum Caliber (MaxCal) calibration of the
##  diabetes CTMC, exact-semigroup (ODE) and Gillespie (SSA) simulation,
##  phase-type survival summaries, and sensitivity analyses.
##
##  Base R only (no packages needed).  Run from the project root:
##       Rscript R/DiabetesCtmcToOdeV31_Analysis.R
##  Inputs : Data/DiabetesCtmcToOdeV31_nhanes_average_target.csv, Data/DiabetesCtmcToOdeV31_literature_edge_inputs.csv,
##           Data/DiabetesCtmcToOdeV31_nhanes_cycle_ratios_digitised.csv (cycle-specific ratios for Figure 3)
##  Outputs: Data/DiabetesCtmcToOdeV31_derived_*.csv, Figures/DiabetesCtmcToOdeV31_*.png,
##           Data/DiabetesCtmcToOdeV31_derived_numbers_for_tex.txt (all manuscript numbers, 4 significant digits)
##
##  V31 changes (figures only; the computations are those of V22):
##   - Figure 3 (NHANES ratios) is drawn here, with x-axis label "Time [years]".
##   - Figure 7 (elasticities) has a colour legend and Greek pi labels.
##   - All figure labels use "Exact completion" (no internal version names)
##     and plotmath symbols for Q^(0), Q*, x_D(t); cell values use 4 significant digits.
## =====================================================================

options(digits = 10)
set.seed(42)

dir.create("Data", showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)

## ---------------------------------------------------------------------
## 0. Inputs
## ---------------------------------------------------------------------
target <- read.csv("Data/DiabetesCtmcToOdeV31_nhanes_average_target.csv", stringsAsFactors = FALSE)
pi_raw <- setNames(target$average_proportion, target$state)
pi_cal <- pi_raw / sum(pi_raw)            # renormalise (reported values sum to 1.000001)

lit <- read.csv("Data/DiabetesCtmcToOdeV31_literature_edge_inputs.csv", stringsAsFactors = FALSE)
## conversion: mean time T -> q = 1/T ; probability p over horizon T -> q = -log(1-p)/T
lit$rate <- ifelse(lit$type == "mean_time", 1 / lit$horizon_years,
                   -log(1 - lit$probability) / lit$horizon_years)
lit$rate_naive <- ifelse(lit$type == "mean_time", NA, lit$probability / lit$horizon_years)
q0 <- setNames(lit$rate, lit$edge)

S_cal <- c("H", "P", "T", "Mi", "Ma", "B")          # calibration (regenerative) chain
edges <- data.frame(
  edge = c("HP", "PT", "TP", "TMi", "TMa", "MiT", "MiB", "MiD", "MaB", "MaD", "BD"),
  from = c("H",  "P",  "T",  "T",   "T",   "Mi",  "Mi",  "Mi",  "Ma",  "Ma",  "B"),
  to   = c("P",  "T",  "P",  "Mi",  "Ma",  "T",   "B",   "H",   "B",   "H",   "H"),
  stringsAsFactors = FALSE)
## In the calibration chain, death edges are re-routed to the healthy reservoir H.
K  <- nrow(edges)
iF <- match(edges$from, S_cal); iT <- match(edges$to, S_cal)

## Stationarity operator: (pi^T Q)_k = sum_e pi_from(e) q_e [1{to=k} - 1{from=k}]
Amat <- function(pi) {
  A <- matrix(0, length(S_cal), K)
  for (e in 1:K) { A[iT[e], e] <- A[iT[e], e] + pi[iF[e]]; A[iF[e], e] <- A[iF[e], e] - pi[iF[e]] }
  A
}
A <- Amat(pi_cal)

buildQcal <- function(q) {
  Q <- matrix(0, 6, 6, dimnames = list(S_cal, S_cal))
  for (e in 1:K) Q[iF[e], iT[e]] <- Q[iF[e], iT[e]] + q[e]
  diag(Q) <- -rowSums(Q); Q
}

## V16 "minimal exact positive calibration" (Table 3 of V16)
q_v16 <- c(HP = 0.010324, PT = 0.202006, TP = 0.032100, TMi = 0.042600, TMa = 0.097225,
           MiT = 0.023500, MiB = 0.010948, MiD = 0.154379, MaB = 0.064900, MaD = 0.2, BD = 0.2)

## ---------------------------------------------------------------------
## 1. Structural checks
## ---------------------------------------------------------------------
out <- list()
out$rankA <- qr(A)$rank
out$resid_v16 <- max(abs(A %*% q_v16))
q_lit_full <- c(HP = NA, q0)
## Best possible residual of the literature rates with q_HP free (least squares in q_HP)
r_lit <- A[, -1] %*% q0[edges$edge[-1]]
qHP_ls <- -sum(A[, 1] * r_lit) / sum(A[, 1]^2)
out$resid_lit_best <- max(abs(A[, 1] * qHP_ls + r_lit))

## Fixed-exit-rate (all lambda at literature) feasibility: balance at Ma forces
## q_TMa = pi_Ma * lambda_Ma / pi_T, which must not exceed lambda_T.
lam0 <- tapply(q0[edges$edge[-1]], edges$from[-1], sum)
out$lambda_lit <- lam0
out$qTMa_required <- unname(pi_cal["Ma"] * lam0["Ma"] / pi_cal["T"])
out$lambdaT_lit <- unname(lam0["T"])
out$fixed_lambda_infeasible <- out$qTMa_required > out$lambdaT_lit

## ---------------------------------------------------------------------
## 2. MaxCal calibration (relative-entropy-rate I-projection)
##    min  sum_{e in ref} pi_from [ q log(q/q0) - q + q0 ]
##    s.t. pi^T Q_cal = 0,   q >= 0,   fixed edges held at given values.
##    KKT (interior): q_e = q0_e * exp(u_to - u_from) on referenced free edges;
##                    u_to = u_from on unreferenced (reservoir) edges.
##    Solved by Newton's method on (u, unreferenced rates).
## ---------------------------------------------------------------------
maxcal <- function(pi, q0, fixed = character(0), fixed_val = NULL, noref = "HP",
                   tol = 1e-13, maxit = 200) {
  A <- Amat(pi)
  ref_free <- setdiff(edges$edge, c(fixed, noref))
  ## tie the potentials of unreferenced edges (KKT); gauge u_H = 0
  ## unknown vector: u (6 entries) and q_noref
  rates <- function(par) {
    u <- par[1:6]; q <- setNames(numeric(K), edges$edge)
    for (e in 1:K) {
      nm <- edges$edge[e]
      if (nm %in% ref_free) q[e] <- q0[nm] * exp(u[iT[e]] - u[iF[e]])
      if (nm %in% fixed)    q[e] <- fixed_val[nm]
    }
    q[noref] <- par[7:(6 + length(noref))]
    q
  }
  F <- function(par) {
    q <- rates(par); u <- par[1:6]
    bal <- as.vector(A %*% q)[-1]                 # drop the redundant node (H)
    tie <- sapply(noref, function(nm) { e <- match(nm, edges$edge); u[iT[e]] - u[iF[e]] })
    c(u[1], bal, tie)                              # gauge, 5 balances, ties
  }
  par <- c(rep(0, 6), rep(0.01, length(noref)))
  for (it in 1:maxit) {
    f <- F(par)
    if (max(abs(f)) < tol) break
    J <- sapply(seq_along(par), function(j) { h <- 1e-7 * max(1, abs(par[j]))
      p2 <- par; p2[j] <- p2[j] + h; (F(p2) - f) / h })
    step <- solve(J, -f)
    ## damped Newton with positivity safeguard on unreferenced rates
    a <- 1
    repeat {
      cand <- par + a * step
      if (all(cand[7:length(cand)] > 0) && max(abs(F(cand))) < max(abs(f)) * (1 - 1e-4 * a)) break
      a <- a / 2; if (a < 1e-12) break
    }
    par <- par + a * step
  }
  q <- rates(par)
  obj <- sum(sapply(ref_free, function(nm) { e <- match(nm, edges$edge)
    pi[iF[e]] * (q[nm] * log(q[nm] / q0[nm]) - q[nm] + q0[nm]) }))
  list(q = q, u = setNames(par[1:6], S_cal), obj = obj, iter = it,
       resid = max(abs(A %*% q)), converged = max(abs(F(par))) < 1e-10)
}

partial_obj <- function(q, pi, q0, ref = edges$edge[-1]) {
  sum(sapply(ref, function(nm) { e <- match(nm, edges$edge)
    pi[iF[e]] * (q[nm] * log(q[nm] / q0[nm]) - q[nm] + q0[nm]) }))
}

mc <- maxcal(pi_cal, q0)
q_star <- mc$q
out$maxcal_obj <- mc$obj; out$maxcal_resid <- mc$resid; out$maxcal_iter <- mc$iter
out$v16_obj <- partial_obj(q_v16, pi_cal, q0)
## KKT check: cycle products q*q/q0q0 over the doubled pairs must equal 1
out$cycle_PT_TP <- unname(q_star["PT"] * q_star["TP"] / (q0["PT"] * q0["TP"]))
out$cycle_TMi_MiT <- unname(q_star["TMi"] * q_star["MiT"] / (q0["TMi"] * q0["MiT"]))

## Brute-force check of optimality with a generic constrained optimiser on the null space
chk <- local({
  ref <- edges$edge[-1]
  N <- qr.Q(qr(t(A)), complete = TRUE)[, (qr(A)$rank + 1):K]   # null-space basis
  f <- function(th) { q <- setNames(as.vector(q_star + N %*% th), edges$edge)
    if (any(q <= 0)) return(1e10); partial_obj(q, pi_cal, q0, ref) }
  o <- optim(rep(0, ncol(N)), f, method = "Nelder-Mead", control = list(reltol = 1e-14, maxit = 2e4))
  c(value = o$value, dist = sqrt(sum((N %*% o$par)^2)))
})
out$nullspace_check <- chk

## ---------------------------------------------------------------------
## 3. Calibration table
## ---------------------------------------------------------------------
lam_of <- function(q) tapply(q[edges$edge[-1]], edges$from[-1], sum)[c("P", "T", "Mi", "Ma", "B")]
calib <- data.frame(edge = edges$edge, from = edges$from, to = edges$to,
                    literature = c(NA, q0), v16_exact = q_v16, maxcal = q_star)
calib$to[calib$edge %in% c("MiD", "MaD", "BD")] <- "D (->H in calibration)"
calib$abs_change <- calib$maxcal - calib$literature
calib$pct_change <- 100 * calib$abs_change / calib$literature
calib$v16_pct_change <- 100 * (calib$v16_exact - calib$literature) / calib$literature
write.csv(calib, "Data/DiabetesCtmcToOdeV31_derived_calibrated_rates.csv", row.names = FALSE)
hold <- data.frame(state = c("P", "T", "Mi", "Ma", "B"),
                   tau_literature = 1 / lam_of(q0 * 1 + 0 * q0)[1:5],
                   tau_v16 = 1 / lam_of(q_v16), tau_maxcal = 1 / lam_of(q_star))
hold$tau_literature <- 1 / tapply(q0, edges$from[-1], sum)[hold$state]
write.csv(hold, "Data/DiabetesCtmcToOdeV31_derived_holding_times.csv", row.names = FALSE)

## ---------------------------------------------------------------------
## 4. Simulation generators, exact semigroup, survival (phase-type)
## ---------------------------------------------------------------------
S_sim <- c("P", "T", "Mi", "Ma", "B", "D")
buildQsim <- function(q) {
  Q <- matrix(0, 6, 6, dimnames = list(S_sim, S_sim))
  Q["P", "T"] <- q["PT"]; Q["T", "P"] <- q["TP"]; Q["T", "Mi"] <- q["TMi"]; Q["T", "Ma"] <- q["TMa"]
  Q["Mi", "T"] <- q["MiT"]; Q["Mi", "B"] <- q["MiB"]; Q["Mi", "D"] <- q["MiD"]
  Q["Ma", "B"] <- q["MaB"]; Q["Ma", "D"] <- q["MaD"]; Q["B", "D"] <- q["BD"]
  diag(Q) <- -rowSums(Q); Q
}
expmTS <- function(M) {                          # scaling-and-squaring Taylor matrix exponential
  nrm <- max(rowSums(abs(M))); s <- max(0, ceiling(log2(nrm)) + 4)
  X <- M / 2^s; E <- diag(nrow(M)); term <- diag(nrow(M))
  for (k in 1:30) { term <- term %*% X / k; E <- E + term }
  for (i in seq_len(s)) E <- E %*% E
  E
}
gens <- list(Literature = buildQsim(q0), V16exact = buildQsim(q_v16), MaxCal = buildQsim(q_star))
x0s <- list(AllPrediabetes = c(1, 0, 0, 0, 0, 0),
            DD2baseline = c(0, 0.65, 0.12, 0.17, 0.06, 0))
tgrid <- seq(0, 50, by = 0.1)
exact_traj <- function(Q, x0) t(sapply(tgrid, function(t) as.vector(x0 %*% expmTS(Q * t))))

surv <- data.frame()
for (g in names(gens)) for (ic in names(x0s)) {
  Q <- gens[[g]]; Tm <- Q[1:5, 1:5]; a <- x0s[[ic]][1:5]
  Ninv <- solve(-Tm)
  occ <- as.vector(a %*% Ninv)
  Sfun <- function(t) sum(a %*% expmTS(Tm * t))
  surv <- rbind(surv, data.frame(generator = g, initial = ic,
     mean_time_to_death = sum(occ),
     occ_P = occ[1], occ_T = occ[2], occ_Mi = occ[3], occ_Ma = occ[4], occ_B = occ[5],
     mort_5 = 1 - Sfun(5), mort_10 = 1 - Sfun(10), mort_20 = 1 - Sfun(20), mort_30 = 1 - Sfun(30),
     median_survival = uniroot(function(t) Sfun(t) - 0.5, c(0.01, 500))$root))
}
write.csv(surv, "Data/DiabetesCtmcToOdeV31_derived_survival_summaries.csv", row.names = FALSE)

## ---------------------------------------------------------------------
## 5. Gillespie SSA cohort (N = 2000 independent individuals), seed 42
## ---------------------------------------------------------------------
ssa_cohort <- function(Q, x0, N = 2000) {
  counts <- matrix(0, length(tgrid), 6)
  init <- sample(1:6, N, replace = TRUE, prob = x0)
  for (n in 1:N) {
    s <- init[n]; t <- 0; path_t <- 0; path_s <- s
    while (s != 6 && t < max(tgrid)) {
      lam <- -Q[s, s]; t <- t + rexp(1, lam)
      s <- sample(1:6, 1, prob = pmax(Q[s, ], 0) * (1:6 != s))
      path_t <- c(path_t, t); path_s <- c(path_s, s)
    }
    idx <- findInterval(tgrid, path_t)
    st <- path_s[idx]
    counts[cbind(seq_along(tgrid), st)] <- counts[cbind(seq_along(tgrid), st)] + 1
  }
  counts / N
}
set.seed(42)
ssa <- list()
diag_ssa <- data.frame()
for (ic in names(x0s)) {
  ex <- exact_traj(gens$MaxCal, x0s[[ic]])
  sm <- ssa_cohort(gens$MaxCal, x0s[[ic]], 2000)
  ssa[[ic]] <- list(exact = ex, ssa = sm)
  se <- sqrt(pmax(ex * (1 - ex), 0) / 2000)
  diag_ssa <- rbind(diag_ssa, data.frame(initial = ic,
     sup_abs_dev = max(abs(sm - ex)),
     ## normal approximation only where N x (1 - x) >= 10 (binomial counts not too small)
     frac_points_outside_95band = { ok <- 2000 * ex * (1 - ex) >= 10
       mean((abs(sm - ex) > 1.96 * se)[ok]) },
     max_abs_z = { ok <- 2000 * ex * (1 - ex) >= 10; max((abs(sm - ex) / se)[ok]) },
     ## forward-Euler (dt = 0.1) discretisation error, for reference to V16
     euler_sup_err = { x <- x0s[[ic]]; E <- matrix(0, length(tgrid), 6); E[1, ] <- x
       for (k in 2:length(tgrid)) { x <- x + 0.1 * as.vector(x %*% gens$MaxCal); E[k, ] <- x }
       max(abs(E - ex)) }))
  write.csv(data.frame(t = tgrid, exact = ex, ssa = sm), sprintf("Data/DiabetesCtmcToOdeV31_derived_trajectories_MaxCal_%s.csv", ic), row.names = FALSE)
}
write.csv(diag_ssa, "Data/DiabetesCtmcToOdeV31_derived_ssa_diagnostics.csv", row.names = FALSE)

## ---------------------------------------------------------------------
## 6. Sensitivity analyses
## ---------------------------------------------------------------------
sens <- data.frame()
add_sens <- function(label, pi, q0x, fixed = character(0), fixed_val = NULL) {
  r <- tryCatch(maxcal(pi, q0x, fixed, fixed_val), error = function(e) NULL)
  if (is.null(r) || !r$converged || any(r$q <= 0)) {
    sens <<- rbind(sens, data.frame(scenario = label, feasible = FALSE, t(setNames(rep(NA, K), edges$edge)),
                                    objective = NA, mean_death_prediab = NA, mean_death_DD2 = NA)); return(invisible())
  }
  Q <- buildQsim(r$q); Tm <- Q[1:5, 1:5]
  md <- function(a) sum(a %*% solve(-Tm))
  sens <<- rbind(sens, data.frame(scenario = label, feasible = TRUE, t(r$q), objective = r$obj,
                                  mean_death_prediab = md(x0s$AllPrediabetes[1:5]),
                                  mean_death_DD2 = md(x0s$DD2baseline[1:5])))
}
add_sens("Base case", pi_cal, q0)
for (s in c("P", "T", "Mi", "Ma", "B")) for (f in c(0.9, 1.1)) {
  p <- pi_cal; p[s] <- p[s] * f; p["H"] <- p["H"] - (sum(p) - 1)
  add_sens(sprintf("pi_%s x %.1f", s, f), p, q0)
}
q0b <- q0; q0b[c("MaD", "BD")] <- 1 / 10
add_sens("Mortality prior Ma,B: mean 10 y", pi_cal, q0b)
q0c <- q0; q0c["PT"] <- 1 / 5
add_sens("Prior P->T: mean 5 y", pi_cal, q0c)
add_sens("Mortality edges fixed at literature", pi_cal, q0,
         fixed = c("MiD", "MaD", "BD"), fixed_val = q0[c("MiD", "MaD", "BD")])
add_sens("Remission edges fixed at literature", pi_cal, q0,
         fixed = c("TP", "MiT"), fixed_val = q0[c("TP", "MiT")])
add_sens("Only the two 5-y placeholders fixed", pi_cal, q0,
         fixed = c("MaD", "BD"), fixed_val = q0[c("MaD", "BD")])
q0d <- q0; q0d["MiT"] <- -log(1 - 0.046) / 7       # Karter horizon corrected to 7 years
add_sens("Mi->T reference on a 7-year horizon", pi_cal, q0d)
pj <- pi_cal * c(1, 1.1, 0.9, 0.9, 1.1, 1.1); pj["H"] <- pi_cal["H"] - (sum(pj) - sum(pi_cal))
pj <- pj / sum(pj)
add_sens("Joint perturbation (+/-10% in all five)", pj, q0)
write.csv(sens, "Data/DiabetesCtmcToOdeV31_derived_sensitivity.csv", row.names = FALSE)

## Elasticities: directional derivatives along pi -> pi + eps*pi_k (e_k - e_H),
## i.e. a relative change of one disease component compensated in the healthy share.
## Computed by CENTRED differences of re-solved calibrations (step 1e-4).
elas <- sapply(c("P", "T", "Mi", "Ma", "B"), function(s) {
  h <- 1e-4
  pp <- pi_cal; pp[s] <- pp[s] * (1 + h); pp["H"] <- pp["H"] - (sum(pp) - 1)
  pm <- pi_cal; pm[s] <- pm[s] * (1 - h); pm["H"] <- pm["H"] - (sum(pm) - 1)
  (log(maxcal(pp, q0)$q) - log(maxcal(pm, q0)$q)) / (2 * h) })
write.csv(round(elas, 6), "Data/DiabetesCtmcToOdeV31_derived_elasticities.csv")
out$elas_max_abs <- max(abs(elas))
out$elas_rowsum_norm <- max(rowSums(abs(elas)))            # induced max-norm of the elasticity map

## Independent check: analytic differentiation of the bordered KKT system
## [ W  -A^T ; A  0 ] (dy, du) = ( -d_theta grad_y J + (dA)^T u ; -(dA) y )
elas_kkt <- sapply(c("P", "T", "Mi", "Ma", "B"), function(s) {
  dpi <- rep(0, 6); names(dpi) <- S_cal; dpi[s] <- pi_cal[s]; dpi["H"] <- -pi_cal[s]
  dA <- Amat(dpi)                                          # A is linear in pi
  Ar <- A[-1, ]; dAr <- dA[-1, ]                           # drop the redundant row
  W <- diag(sapply(seq_len(K), function(e)
        if (edges$edge[e] == "HP") 0 else pi_cal[iF[e]] / q_star[e]))
  u <- mc$u
  dgradJ <- sapply(seq_len(K), function(e)
        if (edges$edge[e] == "HP") 0 else dpi[iF[e]] * log(q_star[e] / q0[edges$edge[e]]))
  rhsu <- as.vector(t(dA) %*% u)                           # (dA)^T u, full node set
  M <- rbind(cbind(W, -t(Ar)), cbind(Ar, matrix(0, nrow(Ar), nrow(Ar))))
  rhs <- c(-dgradJ + rhsu, -as.vector(dAr %*% q_star))
  sol <- solve(M, rhs)
  (sol[1:K]) / q_star })
out$elas_max_discrepancy <- max(abs(elas - elas_kkt))
write.csv(round(elas_kkt, 6), "Data/DiabetesCtmcToOdeV31_derived_elasticities_kkt.csv")

## ---------------------------------------------------------------------
## 7. Figures (PNG)
## ---------------------------------------------------------------------
png_dev <- function(file, width, height, res = 200) png(file, width = width, height = height, res = res, type = "cairo")
cols <- c(P = "#9ecae1", T = "#08306b", Mi = "#e6ab02", Ma = "#d95f02", B = "#e41a1c", D = "black")
labs <- c(P = "Prediabetes", T = "T2D (uncomplicated)", Mi = "Microvascular only",
          Ma = "Macrovascular only", B = "Micro + macro", D = "Death")
plot_ode_ssa <- function(ic, file, main, lx = 11, ly = 0.84) {
  ex <- ssa[[ic]]$exact; sm <- ssa[[ic]]$ssa
  png_dev(file, width = 2000, height = 1000)
  par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
  plot(NA, xlim = c(0, 50), ylim = c(0, 1), xlab = "Time [years]", ylab = "Population proportion", main = main)
  grid()
  for (j in 1:6) {
    se <- sqrt(pmax(ex[, j] * (1 - ex[, j]), 0) / 2000)
    polygon(c(tgrid, rev(tgrid)), c(ex[, j] - 1.96 * se, rev(ex[, j] + 1.96 * se)),
            col = adjustcolor(cols[j], 0.25), border = NA)
    lines(tgrid, ex[, j], col = cols[j], lwd = 2)
    lines(tgrid, sm[, j], col = cols[j], lwd = 1.2, lty = 2)
  }
  legend("topright", legend = labs, col = cols, lwd = 2, bty = "n", cex = 0.72, ncol = 2)
  legend(x = lx, y = ly, legend = c(expression("exact " * x(0) * e^{t * Q^"*"}), "SSA cohort mean (N = 2000)", "95% band"),
         lty = c(1, 2, NA), pch = c(NA, NA, 15), col = c("grey20", "grey20", adjustcolor("grey20", 0.3)),
         bty = "n", cex = 0.75)
  ## aggregated burdens
  agg_ex <- cbind(ex[, 3] + ex[, 5], ex[, 4] + ex[, 5]); agg_sm <- cbind(sm[, 3] + sm[, 5], sm[, 4] + sm[, 5])
  plot(NA, xlim = c(0, 50), ylim = c(0, max(agg_ex, agg_sm) * 1.15), xlab = "Time [years]",
       ylab = "Population proportion", main = "Aggregated complication burden")
  grid()
  cc <- c(cols["Mi"], cols["Ma"])
  for (j in 1:2) {
    se <- sqrt(pmax(agg_ex[, j] * (1 - agg_ex[, j]), 0) / 2000)
    polygon(c(tgrid, rev(tgrid)), c(agg_ex[, j] - 1.96 * se, rev(agg_ex[, j] + 1.96 * se)),
            col = adjustcolor(cc[j], 0.25), border = NA)
    lines(tgrid, agg_ex[, j], col = cc[j], lwd = 2); lines(tgrid, agg_sm[, j], col = cc[j], lty = 2, lwd = 1.2)
  }
  legend("topright", legend = c("Microvascular (incl. combined)", "Macrovascular (incl. combined)"),
         col = cc, lwd = 2, bty = "n", cex = 0.8)
  dev.off()
}
plot_ode_ssa("AllPrediabetes", "Figures/DiabetesCtmcToOdeV31_ode_ssa_all_prediabetes.png", "MaxCal generator: all-prediabetes cohort")
plot_ode_ssa("DD2baseline", "Figures/DiabetesCtmcToOdeV31_ode_ssa_dd2_baseline.png", "MaxCal generator: DD2-type baseline", lx = 22, ly = 0.47)

## literature vs calibrated comparison
png_dev("Figures/DiabetesCtmcToOdeV31_literature_vs_calibrated.png", width = 2200, height = 1500)
par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1))
gcol <- c(Literature = "grey40", V16exact = "#1b9e77", MaxCal = "#7570b3")
glty <- c(Literature = 1, V16exact = 3, MaxCal = 2)
glab <- c(expression("Literature " * Q^{(0)}), "Exact completion", expression("MaxCal " * Q^"*"))
for (ic in names(x0s)) {
  trajs <- lapply(gens, exact_traj, x0 = x0s[[ic]])
  ttl <- if (ic == "AllPrediabetes") "all-prediabetes" else "DD2-type baseline"
  plot(NA, xlim = c(0, 50), ylim = c(0, 1), xlab = "Time [years]", ylab = expression("Survival " * 1 - x[D](t)),
       main = paste("Survival,", ttl)); grid()
  for (g in names(gens)) lines(tgrid, 1 - trajs[[g]][, 6], col = gcol[g], lty = glty[g], lwd = 2)
  legend("topright", glab, col = gcol, lty = glty, lwd = 2, bty = "n", cex = 0.8, y.intersp = 1.35)
  plot(NA, xlim = c(0, 50), ylim = c(0, max(sapply(trajs, function(x) max(x[, 2])))), xlab = "Time [years]",
       ylab = "Proportion", main = paste("Uncomplicated T2D,", ttl)); grid()
  for (g in names(gens)) lines(tgrid, trajs[[g]][, 2], col = gcol[g], lty = glty[g], lwd = 2)
  plot(NA, xlim = c(0, 50), ylim = c(0, max(sapply(trajs, function(x) max(x[, 4] + x[, 5], x[, 3] + x[, 5])))),
       xlab = "Time [years]", ylab = "Proportion", main = paste("Complication burden,", ttl), cex.main = 0.95); grid()
  legend("topright", c("micro (incl. combined)", "macro (incl. combined)"), lwd = c(1, 2.5), bty = "n", cex = 0.8)
  for (g in names(gens)) { lines(tgrid, trajs[[g]][, 3] + trajs[[g]][, 5], col = gcol[g], lty = glty[g], lwd = 1)
    lines(tgrid, trajs[[g]][, 4] + trajs[[g]][, 5], col = gcol[g], lty = glty[g], lwd = 2.5) }
}
dev.off()

## sensitivity (elasticity) heat map, with colour legend and Greek pi labels
fmt4 <- function(x) formatC(signif(x, 4), digits = 4, format = "fg", flag = "#")
png_dev("Figures/DiabetesCtmcToOdeV31_elasticities.png", width = 1800, height = 1150)
layout(matrix(1:2, 1, 2), widths = c(5.2, 1))
par(mar = c(4.5, 5, 3, 0.5))
E <- elas; lim <- max(abs(E)); ncol_pal <- 41; pal <- hcl.colors(ncol_pal, "Blue-Red 3")
image(1:ncol(E), 1:nrow(E), t(E), col = pal, zlim = c(-lim, lim),
      axes = FALSE, xlab = "Perturbed target component", ylab = "",
      main = expression("Elasticity  " * partialdiff * log * " " * q[e]^"*" * " / " * partialdiff * log * " " * pi[k]))
xl <- lapply(colnames(E), function(k) bquote(pi[.(k)]))
axis(1, at = 1:ncol(E), labels = as.expression(xl), cex.axis = 1.1, padj = 0.3)
rl <- c(HP = "H*P", PT = "P*T", TP = "T*P", TMi = "T*Mi", TMa = "T*Ma", MiT = "Mi*T", MiB = "Mi*B",
        MiD = "Mi*D", MaB = "Ma*B", MaD = "Ma*D", BD = "B*D")
yl <- lapply(rownames(E), function(r) parse(text = paste0("q[", rl[r], "]"))[[1]])
axis(2, at = 1:nrow(E), labels = as.expression(yl), las = 1)
box()
for (i in 1:nrow(E)) for (j in 1:ncol(E))
  text(j, i, fmt4(E[i, j]), cex = 0.62, col = if (abs(E[i, j]) > 0.55 * lim) "white" else "black")
## colour legend (colour bar)
par(mar = c(4.5, 0.8, 3, 4.2))
zb <- seq(-lim, lim, length.out = ncol_pal + 1)
image(1, (zb[-1] + zb[-length(zb)]) / 2, matrix((zb[-1] + zb[-length(zb)]) / 2, nrow = 1), col = pal,
      zlim = c(-lim, lim), axes = FALSE, xlab = "", ylab = "")
axis(4, at = pretty(c(-lim, lim)), las = 1, cex.axis = 0.9)
mtext("Elasticity", side = 4, line = 3, cex = 0.9)
box()
dev.off()

## NHANES cycle-specific composition among non-healthy participants (Figure 3)
## Values digitised from the V22 figure (the participant-level classification is not in the project folder).
nh <- read.csv("Data/DiabetesCtmcToOdeV31_nhanes_cycle_ratios_digitised.csv")
ncl <- c(P = "#1f77b4", T = "#ff7f0e", Ma = "#2ca02c", Mi = "#d62728", B = "#9467bd")
nlb <- c(P = "Prediabetes", T = "Diabetes without complications", Ma = "Diabetes macro only",
         Mi = "Diabetes micro only", B = "Diabetes both macro and micro")
png_dev("Figures/DiabetesCtmcToOdeV31_nhanes_state_ratios.png", width = 1800, height = 1000)
par(mar = c(4.2, 4.5, 1, 1))
plot(NA, xlim = range(nh$cycle_label), ylim = c(0, 0.47), xaxt = "n",
     xlab = "Year", ylab = "Proportion among non-healthy participants")
axis(1, at = nh$cycle_label); abline(v = nh$cycle_label, h = seq(0, 0.4, 0.1), col = "grey90")
for (k in names(ncl)) lines(nh$cycle_label, nh[[k]], type = "o", pch = 19, col = ncl[k], lwd = 2)
legend("right", legend = nlb, col = ncl, lwd = 2, pch = 19, bg = "white", cex = 0.9)
box()
dev.off()

## ---------------------------------------------------------------------
## 8. Dump numbers used in the manuscript
## ---------------------------------------------------------------------
sink("Data/DiabetesCtmcToOdeV31_derived_numbers_for_tex.txt")
cat("pi_cal (renormalised):\n"); print(round(pi_cal, 6))
cat("\nLiterature rates (exact conversion):\n"); print(round(q0, 6))
cat("\nNaive p/T rates:\n"); print(round(setNames(lit$rate_naive, lit$edge), 6))
str(out)
cat("\nCalibration table:\n"); print(calib, digits = 6)
cat("\nHolding times:\n"); print(hold, digits = 5)
cat("\nMaxCal potentials u:\n"); print(round(mc$u, 6)); cat("exp(u):\n"); print(round(exp(mc$u), 6))
cat("\nSurvival summaries:\n"); print(surv, digits = 4)
cat("\nSSA diagnostics:\n"); print(diag_ssa, digits = 4)
cat("\nSensitivity:\n"); print(sens, digits = 4)
cat("\nElasticities:\n"); print(round(elas, 3))
sink()
cat("Done.\n")

## ---------------------------------------------------------------------
## 9. SSA validation over replicate cohorts (average exceedance rate)
## ---------------------------------------------------------------------
set.seed(4242)
rep_exc <- sapply(1:20, function(r) {
  ex <- ssa$AllPrediabetes$exact; sm <- ssa_cohort(gens$MaxCal, x0s$AllPrediabetes, 2000)
  se <- sqrt(pmax(ex * (1 - ex), 0) / 2000); ok <- 2000 * ex * (1 - ex) >= 10
  mean((abs(sm - ex) > 1.96 * se)[ok]) })
cat(sprintf("Replicate-cohort mean pointwise exceedance of 95%% band: %.4f (range %.4f-%.4f)\n",
            mean(rep_exc), min(rep_exc), max(rep_exc)),
    file = "Data/DiabetesCtmcToOdeV31_derived_numbers_for_tex.txt", append = TRUE)

## ---------------------------------------------------------------------
## 10. Checks requested in the technical review of V17
## ---------------------------------------------------------------------
chk <- list()

## (a) Optimality certificate on full-precision output
chk$min_rate            <- min(q_star)
chk$balance_resid       <- max(abs(A %*% q_star))
chk$tilt_resid          <- max(sapply(setdiff(edges$edge, "HP"), function(nm) {
  e <- match(nm, edges$edge)
  abs(log(q_star[nm] / q0[nm]) - (mc$u[iT[e]] - mc$u[iF[e]])) }))
chk$free_tie_resid      <- abs(mc$u["P"] - mc$u["H"])
chk$constraint_rank     <- qr(A)$rank
## Consistency check (not an exact certificate): the Pythagorean identity
## J(q) - J(q*) = sum pi_i g(q, q*) on a perturbation q* + N v with A N = 0.
## Such a perturbation inherits the small balance residual of q*, so the agreement
## below verifies the algebra and the implementation, not exact feasibility.
Nsp <- qr.Q(qr(t(A)), complete = TRUE)[, (chk$constraint_rank + 1):K]
qpert <- as.vector(q_star + 0.02 * Nsp %*% rep(1, ncol(Nsp))); names(qpert) <- edges$edge
chk$pythagoras_gap <- if (all(qpert > 0)) {
  lhs <- partial_obj(qpert, pi_cal, q0) - partial_obj(q_star, pi_cal, q0)
  rhs <- partial_obj(qpert, pi_cal, q_star)
  abs(lhs - rhs) } else NA
## V16 printed (rounded) vector is not exactly feasible
chk$v16_balance_resid <- max(abs(A %*% q_v16))

## (b) Literature throughputs pi_i * lambda_i^(0) of the complication states
lam0all <- tapply(q0, edges$from[-1], sum)
chk$throughputs <- c(Mi = pi_cal["Mi"] * lam0all["Mi"], Ma = pi_cal["Ma"] * lam0all["Ma"],
                     B = pi_cal["B"] * lam0all["B"])
chk$lit_resid_lsq <- out$resid_lit_best

## (c) Renewal identity  e_P (-U)^{-1} = pi_L / (pi_H q_HP)
U <- buildQsim(q_star)[1:5, 1:5]
chk$renewal_occupation <- as.vector(solve(-t(U), c(1, 0, 0, 0, 0)))   # e_P (-U)^{-1}
chk$renewal_target     <- as.vector(pi_cal[c("P", "T", "Mi", "Ma", "B")] / (pi_cal["H"] * q_star["HP"]))
chk$renewal_max_err    <- max(abs(chk$renewal_occupation - chk$renewal_target))
chk$r_flux             <- unname(pi_cal["H"] * q_star["HP"])
chk$ExpTauD_renewal    <- unname((1 - pi_cal["H"]) / chk$r_flux)

## (d) Scale invariance of the calibration: references x c  =>  rates x c
mc_c <- maxcal(pi_cal, 2 * q0)
chk$scale_invariance_err <- max(abs(mc_c$q / q_star - 2))
## Invariance to the healthy share at fixed relative disease composition
s_old <- 1 - pi_cal["H"]; pi_s <- pi_cal
pi_s[c("P","T","Mi","Ma","B")] <- pi_cal[c("P","T","Mi","Ma","B")] * 1.2
pi_s["H"] <- 1 - sum(pi_s[c("P","T","Mi","Ma","B")])
mc_s <- maxcal(pi_s, q0)
chk$healthy_share_invariance_err <- max(abs(mc_s$q[-1] / q_star[-1] - 1))
chk$qHP_scaling <- unname(mc_s$q["HP"] / q_star["HP"])
chk$qHP_scaling_predicted <- unname(((1 - pi_s["H"]) / pi_s["H"]) / ((1 - pi_cal["H"]) / pi_cal["H"]))

## (e) Largest one-at-a-time change of a disease rate over the +/-10% target scenarios
sens_oat <- sens[grepl("^pi_", sens$scenario), edges$edge]
chk$max_oat_pct_change <- 100 * max(abs(sweep(as.matrix(sens_oat), 2, q_star, "/") - 1))
joint <- sens[sens$scenario == "Joint perturbation (+/-10% in all five)", edges$edge]
chk$joint_pct_change <- 100 * (as.numeric(joint) / q_star - 1); names(chk$joint_pct_change) <- edges$edge

## (f) Coverage of the normal band at a small proportion (N = 2000)
p_small <- 5e-4; band <- 1.96 * sqrt(p_small * (1 - p_small) / 2000)
kk <- 0:2000; inb <- abs(kk / 2000 - p_small) <= band
chk$normal_band_coverage_at_p0.0005 <- sum(dbinom(kk[inb], 2000, p_small))

sink("Data/DiabetesCtmcToOdeV31_derived_review_checks.txt"); str(chk); print(chk); sink()
cat("Review checks written.\n")

## ---------------------------------------------------------------------
## 11. Covariance scaling check (V19 review, item M1)
##     V(t) = diag(x) - x'x solves Vdot = Q'V + VQ + Sigma(x); the covariance
##     of the empirical proportions is V(t)/N.
## ---------------------------------------------------------------------
Qs <- gens$MaxCal
Sigma_of <- function(x) { M <- matrix(0, 6, 6)
  for (i in 1:6) for (j in 1:6) if (i != j && Qs[i, j] > 0) {
    v <- rep(0, 6); v[j] <- 1; v[i] <- v[i] - 1
    M <- M + x[i] * Qs[i, j] * (v %o% v) }
  M }
resid_cov <- sapply(c(0, 5, 20), function(tt) {
  x <- as.vector(x0s$AllPrediabetes %*% expmTS(Qs * tt))
  S <- diag(x) - outer(x, x)
  xdot <- as.vector(x %*% Qs)                       # = x Q
  Sdot <- diag(xdot) - outer(xdot, x) - outer(x, xdot)
  max(abs(Sdot - (t(Qs) %*% S + S %*% Qs + Sigma_of(x)))) })
cat(sprintf("Covariance-ODE residuals at t = 0, 5, 20 y: %s\n",
            paste(format(resid_cov, digits = 3), collapse = ", ")),
    file = "Data/DiabetesCtmcToOdeV31_derived_numbers_for_tex.txt", append = TRUE)

## ---------------------------------------------------------------------
## 12. (V31) All manuscript numbers to four significant digits
## ---------------------------------------------------------------------
s4 <- function(x) signif(x, 4)
Qlit_cal <- buildQcal(c(HP = unname(q_star["HP"]), q0))
pi_litchain <- local({ M <- rbind(t(Qlit_cal), rep(1, 6)); qr.solve(M, c(rep(0, 6), 1)) })
names(pi_litchain) <- S_cal
lam_star <- lam_of(q_star); lam_lit <- tapply(q0, edges$from[-1], sum)[c("P","T","Mi","Ma","B")]
v31 <- list(
  pi_cal = s4(pi_cal), pi_raw_sum = sum(pi_raw), healthy_pct = s4(100 * pi_cal["H"]),
  prediab_pct = s4(100 * pi_cal["P"]),
  q0 = s4(q0), q_star = s4(q_star), q_exact_completion = s4(q_v16),
  pct_change_maxcal = s4(calib$pct_change), pct_change_exact = s4(calib$v16_pct_change),
  hold = data.frame(state = hold$state, lit = s4(hold$tau_literature), maxcal = s4(hold$tau_maxcal), exact = s4(hold$tau_v16)),
  exp_u = s4(exp(mc$u)), newton_iter = mc$iter,
  inflow_TMa_factor = s4(exp(mc$u["Ma"] - mc$u["T"])),
  pi_literature_chain = s4(pi_litchain),
  ratio_Mi = s4(pi_litchain["Mi"] / pi_cal["Mi"]), ratio_Ma = s4(pi_cal["Ma"] / pi_litchain["Ma"]),
  lit_resid_best = s4(out$resid_lit_best), throughputs = s4(chk$throughputs),
  lambda_lit = s4(lam_lit), lambda_Ma_lit = s4(lam_lit["Ma"]),
  qTMa_required = s4(out$qTMa_required), lambdaT_lit = s4(out$lambdaT_lit),
  obj_maxcal = s4(mc$obj), obj_exact = s4(out$v16_obj), obj_ratio = s4(out$v16_obj / mc$obj),
  resid_exact_rounded = s4(out$resid_v16),
  min_rate = s4(chk$min_rate), balance_resid = s4(chk$balance_resid), tilt_resid = s4(chk$tilt_resid),
  nullspace_check = s4(out$nullspace_check), pythagoras_gap = s4(chk$pythagoras_gap),
  surv = transform(surv, mean_time_to_death = s4(mean_time_to_death), median_survival = s4(median_survival),
                   occ_P = s4(occ_P), occ_T = s4(occ_T), occ_Mi = s4(occ_Mi), occ_Ma = s4(occ_Ma), occ_B = s4(occ_B),
                   mort_5 = s4(100 * mort_5), mort_10 = s4(100 * mort_10), mort_20 = s4(100 * mort_20), mort_30 = s4(100 * mort_30)),
  r_flux = s4(chk$r_flux), ExpTauD_renewal = s4(chk$ExpTauD_renewal), renewal_max_err = s4(chk$renewal_max_err),
  renewal_occ = s4(chk$renewal_target),
  ssa = transform(diag_ssa, sup_abs_dev = s4(sup_abs_dev), frac = s4(100 * frac_points_outside_95band), max_abs_z = s4(max_abs_z)),
  rep_exc = s4(100 * c(mean = mean(rep_exc), min = min(rep_exc), max = max(rep_exc))),
  elas_max_abs = s4(out$elas_max_abs), elas_rowsum_norm = s4(out$elas_rowsum_norm),
  elas_kkt_discrepancy = s4(out$elas_max_discrepancy),
  max_oat_pct = s4(chk$max_oat_pct_change), joint_pct = s4(chk$joint_pct_change),
  sens = cbind(scenario = sens$scenario, as.data.frame(lapply(sens[, -(1:2)], s4))),
  normal_cov = s4(chk$normal_band_coverage_at_p0.0005),
  rare_understatement_pct = s4(100 * (1 - 0.38 / -log(1 - 0.38))),
  MiT_7y = s4(-log(1 - 0.046) / 7))
sink("Data/DiabetesCtmcToOdeV31_derived_numbers_4sf.txt"); options(width = 200); print(v31); sink()
cat("V31 4-s.f. numbers written.\n")
