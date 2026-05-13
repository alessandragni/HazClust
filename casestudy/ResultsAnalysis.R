# ============================================================
# Section 4.3
#   1. Table 5  – Cluster composition (covariates + mortality)
#   2. Figure A – KM curves stratified by cluster
#   3. Figure B – Frailty-adjusted linear predictor by cluster
#   4. Figure C – Patient similarity graph coloured by cluster
# ============================================================
library(roxygen2)
roxygenise()

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(cowplot)
library(igraph)
library(ggraph)
library(tidygraph)
library(kableExtra)   # for LaTeX table; use knitr::kable() if not available

setwd("~/Documents/DATA/POLITECNICO/PHD/CODE_REPO/HazClust")
# ── 0. Load data ────────────────────────────────────────────
sim_data <- readRDS("casestudy/data/sim_data.Rdata")
xx       <- readRDS("casestudy/result/output_gammapar0.01_frailtyingau_baseweibull_k150_c3.Rdata")

# Rebuild RESP if not already present
sim_data$RESP <- as.numeric(sim_data$BRONCHITE) +
  (as.numeric(sim_data$INSUF_RESPIRATORIA) - 1) +
  (as.numeric(sim_data$COPD)               - 1) +
  (as.numeric(sim_data$POLMONITE)           - 1)

sim_data$clusters <- as.factor(xx$clusters)

# Cluster labels (adjust order if needed after inspecting mortality rates)
levels(sim_data$clusters) <- c("Cluster 1", "Cluster 2", "Cluster 3")

# Dark2 palette – consistent across all figures
clust_colors <- c("Cluster 1" = "#1B9E77",
                  "Cluster 2" = "#D95F02",
                  "Cluster 3" = "#7570B3")


# ── 1. CLUSTER COMPOSITION TABLE ────────────────────────────
# Compute summary statistics per cluster

comp_table <- sim_data %>%
  group_by(clusters) %>%
  summarise(
    N             = n(),
    Age_mean      = mean(ETA_AL_RICOVERO, na.rm = TRUE),
    Age_sd        = sd(ETA_AL_RICOVERO,   na.rm = TRUE),
    Male_pct      = mean(SESSO == "M",    na.rm = TRUE) * 100,
    ModMCS_mean   = mean(MCS,             na.rm = TRUE),
    ModMCS_sd     = sd(MCS,               na.rm = TRUE),
    Resp_mean     = mean(RESP,            na.rm = TRUE),
    Resp_sd       = sd(RESP,              na.rm = TRUE),
    Mortality_pct = mean(cens,            na.rm = TRUE) * 100,
    Time_mean     = mean(time,            na.rm = TRUE),
    Time_sd       = sd(time,              na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Age     = sprintf("%.1f (%.1f)", Age_mean,    Age_sd),
    ModMCS  = sprintf("%.1f (%.1f)", ModMCS_mean, ModMCS_sd),
    Resp    = sprintf("%.2f (%.2f)", Resp_mean,   Resp_sd),
    Male    = sprintf("%.1f%%",      Male_pct),
    Mort    = sprintf("%.1f%%",      Mortality_pct),
    Time    = sprintf("%.1f (%.1f)", Time_mean,   Time_sd)
  ) %>%
  select(clusters, N, Age, Male, ModMCS, Resp, Mort, Time)

names(comp_table) <- c("Cluster", "N",
                       "Age, mean (SD)",
                       "Male, %",
                       "ModMCS, mean (SD)",
                       "Resp, mean (SD)",
                       "90-day mortality",
                       "Follow-up time, mean (SD)")

# Print to console
print(comp_table)

# LaTeX version for the paper
comp_table %>%
  kbl(booktabs = TRUE, align = "lrcccccc",
      caption = "Cluster composition: covariate profiles and event rates for the three clusters identified in the Enhance-Heart cohort.") %>%
  kable_styling(latex_options = c("hold_position")) %>%
  save_kable("table5_cluster_composition.tex")

# Also save as CSV for convenience
# write.csv(comp_table, "table5_cluster_composition.csv", row.names = FALSE)
# cat("Table saved.\n")


# ── 2. KAPLAN-MEIER CURVES BY CLUSTER ───────────────────────
km_fit <- survfit(Surv(time, cens) ~ clusters, data = sim_data)

fig_km <- ggsurvplot(
  km_fit,
  data          = sim_data,
  palette       = unname(clust_colors),
  linetype      = c("solid", "dashed", "dotdash"),   # B&W-safe
  size          = 0.9,
  risk.table    = TRUE,
  risk.table.height = 0.28,
  risk.table.fontsize = 3.2,
  xlab          = "Time (days)",
  ylab          = "Survival probability",
  legend.title  = "",
  legend.labs   = names(clust_colors),
  conf.int      = TRUE,
  conf.int.alpha = 0.10,
  ggtheme       = theme_classic(base_size = 11) +
    theme(legend.position = "bottom")
)

# Add log-rank p-value annotation
fig_km$plot <- fig_km$plot +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Log-rank p ", 
                          format.pval(surv_pvalue(km_fit, sim_data)$pval, 
                                      digits = 2, eps = 0.001)),
           hjust = 1.05, vjust = 1.5, size = 3.2, fontface = "italic")

# Save
# pdf("figure_km_by_cluster.pdf", width = 7, height = 6)
print(fig_km)
# dev.off()
cat("KM figure saved.\n")

# Note for the paper: if KM curves show partial overlap, add the sentence:
# "KM estimates are marginal and do not condition on the frailty term;
#  the frailty-adjusted risk separation is illustrated in Figure X."


# ── 3. FRAILTY-ADJUSTED LINEAR PREDICTOR BY CLUSTER ─────────
# This is the quantity the clustering actually partitions: x'β + log(û_g)
# Requires: estimated beta, estimated frailties per hospital, 
#           and hospital membership per patient.

# --- Extract beta estimates from model output ---
# Adjust field names to match your xx object structure
beta_hat <- xx$beta          # named vector: Age, Gender, ModMCS, Resp (check names)

# Covariate matrix (must match model formula order)
X <- model.matrix(~ ETA_AL_RICOVERO + SESSO + MCS + RESP - 1,
                  data = sim_data)
# If SESSO is coded differently (e.g. 0/1), adjust accordingly:
# X <- cbind(sim_data$ETA_AL_RICOVERO, 
#            as.numeric(sim_data$SESSO == "M"),
#            sim_data$MCS, 
#            sim_data$RESP)



computeetahat <- function(formula, cluster, strata, data, dist,
                          frailty, p,
                          transform = TRUE,
                          correct){
  
  #----- create obsdata ------------------------------------------#
  
  obs = obsdata_creator(formula = formula, data = data, cluster = cluster, 
                        strata = strata, frailty = frailty, dist = dist)
  nRpar <- obs$nRpar
  nBpar <- obs$nBpar
  nFpar <- obs$nFpar
  
  
  # ---- Assign the number of frailty parameters 'obs$nFpar' --------------- #
  # ---- and compute Sigma for the Positive Stable frailty ----------------- #
  
  if (frailty %in% c("gamma", "ingau")) {
    theta <- ifelse(transform, exp(p[1]), p[1])
  } else if (frailty == "lognormal") {
    sigma2 <- ifelse(transform, exp(p[1]), p[1])
  } else if (frailty == "possta") {
    nu <- ifelse(transform, exp(-exp(p[1])), p[1])
    D <- max(obs$dqi)
    Omega <- Omega(D, correct = correct, nu = nu)
  }
  
  
  # ---- Baseline hazard --------------------------------------------------- #
  if (frailty == 'none') obs$nFpar <- 0
  
  # baseline parameters
  if (dist %in% c("weibull", "inweibull", "frechet")) {
    if (transform) {
      pars <- cbind(rho    = exp(p[obs$nFpar + 1:obs$nstr]),
                    lambda = exp(p[obs$nFpar + obs$nstr + 1:obs$nstr]))
    } else {
      pars <- cbind(rho    = p[obs$nFpar + 1:obs$nstr],
                    lambda = p[obs$nFpar + obs$nstr + 1:obs$nstr])
    }
    beta <- p[-(1:(obs$nFpar + 2 * obs$nstr))]
  } else if (dist == "exponential") {
    if (transform) {
      pars <- cbind(lambda = exp(p[obs$nFpar + 1:obs$nstr]))
    } else {
      pars <- cbind(lambda = p[obs$nFpar + 1:obs$nstr])
    }
    beta <- p[-(1:(obs$nFpar + obs$nstr))]
  } else if (dist == "gompertz") {
    if (transform) {
      pars <- cbind(gamma  = exp(p[obs$nFpar + 1:obs$nstr]),
                    lambda = exp(p[obs$nFpar + obs$nstr + 1:obs$nstr]))
    } else {
      pars <- cbind(gamma  = p[obs$nFpar + 1:obs$nstr],
                    lambda = p[obs$nFpar + obs$nstr + 1:obs$nstr])
    }  
    beta <- p[-(1:(obs$nFpar + 2 * obs$nstr))]
  } else if (dist == "lognormal") {
    if (transform) {
      pars <- cbind(mu    = p[obs$nFpar + 1:obs$nstr],
                    sigma = exp(p[obs$nFpar + obs$nstr + 1:obs$nstr]))
    } else {
      pars <- cbind(mu    = p[obs$nFpar + 1:obs$nstr],
                    sigma = p[obs$nFpar + obs$nstr + 1:obs$nstr])
    }
    beta <- p[-(1:(obs$nFpar + 2 * obs$nstr))]
  } else if (dist == "loglogistic") {
    if (transform) {
      pars <- cbind(alpha = p[obs$nFpar + 1:obs$nstr],
                    kappa = exp(p[obs$nFpar + obs$nstr + 1:obs$nstr]))
    } else  {
      pars <- cbind(alpha = p[obs$nFpar + 1:obs$nstr],
                    kappa = p[obs$nFpar + obs$nstr + 1:obs$nstr])
    }
    beta <- p[-(1:(obs$nFpar + 2 * obs$nstr))]
  } else if (dist == "logskewnormal") {
    if (transform) {
      pars <- cbind(mu    = p[obs$nFpar + 1:obs$nstr],
                    sigma = exp(p[obs$nFpar + obs$nstr + 1:obs$nstr]),
                    alpha = exp(p[obs$nFpar + 2 * obs$nstr + 1:obs$nstr]))
    } else {
      pars <- cbind(mu    = p[obs$nFpar + 1:obs$nstr],
                    sigma = p[obs$nFpar + obs$nstr + 1:obs$nstr],
                    alpha = p[obs$nFpar + 2 * obs$nstr + 1:obs$nstr])
    }
    beta <- p[-(1:(obs$nFpar + 3 * obs$nstr))]
  }
  rownames(pars) <- levels(as.factor(obs$strata))
  
  # baseline: from string to the associated function
  dist <- eval(parse(text = dist))
  
  
  # ---- Cumulative Hazard by cluster and by strata ------------------------- #
  
  cumhaz <- NULL
  cumhaz <- matrix(unlist(
    sapply(levels(as.factor(obs$strata)),
           function(x) {t(
             cbind(dist(pars[x, ], obs$time[obs$strata == x], what = "H"
             ) * exp(as.matrix(obs$x)[
               obs$strata == x, -1, drop = FALSE] %*% as.matrix(beta)),
             obs$cluster[obs$strata == x]))
           })), ncol = 2, byrow = TRUE)
  # cumhaz_i <- cumhaz[,1]
  cumhaz <- aggregate(cumhaz[, 1], by = list(cumhaz[, 2]), 
                      FUN = sum)[, 2, drop = FALSE]
  ### NO FRAILTY
  if (frailty == "none") cumhaz <- sum(cumhaz)
  
  
  
  # ---- log[ (-1)^di L^(di)(cumhaz) ]-------------------------------------- #
  logSurv <- NULL
  if (frailty == "gamma") {
    logSurv <- mapply(fr.gamma, 
                      k = obs$di, s = as.numeric(cumhaz[[1]]), 
                      theta = rep(theta, obs$ncl), 
                      what = "logLT") 
  } else if (frailty == "ingau") {
    logSurv <- mapply(fr.ingau, 
                      k = obs$di, s = as.numeric(cumhaz[[1]]), 
                      theta = rep(theta, obs$ncl), 
                      what = "logLT") 
  } else if (frailty == "possta") {
    logSurv <- sapply(1:obs$ncl, 
                      function(x) fr.possta(k = obs$di[x], 
                                            s = as.numeric(cumhaz[[1]])[x], 
                                            nu = nu, Omega = Omega, 
                                            what = "logLT",
                                            correct = correct))
  } else if (frailty == "lognormal") {
    logSurv <- mapply(fr.lognormal, 
                      k = obs$di, s = as.numeric(cumhaz[[1]]), 
                      sigma2 = rep(sigma2, obs$ncl), 
                      what = "logLT")
  } else if (frailty == "none") {
    logSurv <- mapply(fr.none, s = cumhaz, what = "logLT")
  }
  
  
  # ---- log(\hat u_g) = log[ (-1)^(di+1) L^(di+1)(cumhaz) ] - log[ (-1)^(di) L^(di)(cumhaz) ] --- #
  DIFFlogSurv <- NULL
  if (frailty == "gamma") {
    DIFFlogSurv <- mapply(fr.gamma, 
                          k = obs$di+1, s = as.numeric(cumhaz[[1]]), 
                          theta = rep(theta, obs$ncl), 
                          what = "logLT") - logSurv
  } else if (frailty == "ingau") {
    DIFFlogSurv <- mapply(fr.ingau, 
                          k = obs$di+1, s = as.numeric(cumhaz[[1]]), 
                          theta = rep(theta, obs$ncl), 
                          what = "logLT") - logSurv
  } else if (frailty == "possta") {
    DIFFlogSurv <- sapply(1:obs$ncl, 
                          function(x) fr.possta(k = obs$di[x]+1, 
                                                s = as.numeric(cumhaz[[1]])[x], 
                                                nu = nu, Omega = Omega, 
                                                what = "logLT",
                                                correct = correct)) - logSurv
  } else if (frailty == "lognormal") {
    DIFFlogSurv <- mapply(fr.lognormal, 
                          k = obs$di+1, s = as.numeric(cumhaz[[1]]), 
                          sigma2 = rep(sigma2, obs$ncl), 
                          what = "logLT") - logSurv
  } else if (frailty == "none") {
    DIFFlogSurv <- mapply(fr.none, s = cumhaz, what = "logLT") - logSurv
  }
  
  
  # ---- X'_gi * beta and Xonly_i by cluster --------------------------------------------- #
  Xbeta <- Xbeta_i <- NULL
  if (frailty != "none")  {
    Xbeta <- matrix(unlist(
      sapply(levels(as.factor(obs$strata)),
             function(x) {
               t(cbind( as.matrix(obs$x)[
                 obs$strata == x, -1, drop = FALSE] %*% 
                   as.matrix(beta),
                 obs$cluster[obs$strata == x]))
             })), ncol = 2, byrow = TRUE)
    Xbeta_i <- Xbeta[,1] 
  } else {
    Xbeta_i <- unlist(lapply(1:nrow(pars), function(i) {
      stratum <- rownames(pars)[i]
      x_mat <- as.matrix(obs$x[obs$strata == stratum, -1, drop = FALSE])
      drop(x_mat %*% as.matrix(beta))
    }))
  }
  
  # ---- Frailty term replicated at patient level ---------------------------- #
  log_u_g <- as.vector(
    DIFFlogSurv[
      match(obs[["cluster"]], unique(obs[["cluster"]]))
    ]
  )
  
  # ---- Total frailty-adjusted predictor ------------------------------------ #
  eta_hat <- Xbeta_i + log_u_g
  
  # ---- Return all components clearly --------------------------------------- #
  return(list(
    eta_hat   = eta_hat,     # X beta + log(u_g)
    xbeta     = Xbeta_i,     # X beta only
    log_u_g   = log_u_g      # frailty contribution only
  ))
}

# ── Reconstruct the inputs used when you fitted the model ──────────────────
# These must match exactly what you passed to your fitting function.
# Adjust to match your original model call.

formula  <- Surv(time, cens) ~ ETA_AL_RICOVERO + SESSO + MCS + RESP
cluster  <- "COD_OSPEDALE"
strata   <- NULL          # or whatever strata you used, if any
dist     <- "weibull"
frailty  <- "ingau"
transform <- TRUE         # default in the function

# ── Compute frailty-adjusted linear predictor per patient ──────────────────
eta_obj <- computeetahat(
  formula   = formula,
  cluster   = cluster,
  strata    = strata,
  data      = sim_data,
  dist      = dist,
  frailty   = frailty,
  p         = xx$estimate,
  transform = transform
)

sim_data$lin_pred <- eta_obj$eta_hat
sim_data$xbeta    <- eta_obj$xbeta
sim_data$log_u_g  <- eta_obj$log_u_g




# ── Hospital-level posterior frailty + cluster composition ────────────────

hosp_frailty <- sim_data %>%
  group_by(COD_OSPEDALE) %>%
  summarise(
    log_u_hat = mean(log_u_g),
    u_hat     = exp(mean(log_u_g)),
    n_patients = n(),
    .groups = "drop"
  )

# Compute cluster proportions within hospital
hosp_cluster_prop <- sim_data %>%
  group_by(COD_OSPEDALE, clusters) %>%
  summarise(n_cluster = n(), .groups = "drop") %>%
  left_join(
    sim_data %>%
      group_by(COD_OSPEDALE) %>%
      summarise(n_total = n(), .groups = "drop"),
    by = "COD_OSPEDALE"
  ) %>%
  mutate(prop_cluster = n_cluster / n_total)

# Keep dominant cluster only
dominant_cluster <- hosp_cluster_prop %>%
  group_by(COD_OSPEDALE) %>%
  slice_max(prop_cluster, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(COD_OSPEDALE,
         dominant_cluster = clusters,
         dominant_prop = prop_cluster)

# Merge
hosp_frailty <- hosp_frailty %>%
  left_join(dominant_cluster, by = "COD_OSPEDALE") %>%
  arrange(u_hat) %>%
  mutate(
    hosp_rank = factor(row_number(), levels = row_number()),
    dominant_pct = 100 * dominant_prop
  )




fig_frailty <- ggplot(
  hosp_frailty,
  aes(
    x = hosp_rank,
    y = u_hat,
    colour = dominant_cluster,
    size = dominant_pct
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.5
  ) +
  
  geom_point(alpha = 0.9) +
  
  geom_text(
    aes(label = sprintf("%.0f%%", dominant_pct)),
    vjust = -1,
    size = 2.8,
    show.legend = FALSE
  ) +
  
  scale_colour_manual(
    values = clust_colors,
    name = "Dominant cluster"
  ) +
  
  scale_size_continuous(
    range = c(2.5, 8),
    name = "% patients in dominant cluster"
  ) +
  
  labs(
    x = "Hospital (ranked by posterior frailty)",
    y = expression(hat(u)[g]),
    title = "Posterior frailty estimates by hospital"
  ) +
  
  theme_classic(base_size = 11) +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 11, hjust = 0.5)
  )

fig_frailty




# ── Hospital-level decomposition ------------------------------------------

hosp_summary <- sim_data %>%
  group_by(COD_OSPEDALE) %>%
  summarise(
    mean_xbeta  = mean(xbeta),
    mean_log_u  = mean(log_u_g),
    u_hat       = exp(mean(log_u_g)),
    mean_eta    = mean(lin_pred),
    n_patients  = n(),
    .groups = "drop"
  )


# Dominant cluster info
hosp_summary <- hosp_summary %>%
  left_join(dominant_cluster, by = "COD_OSPEDALE") %>%
  mutate(
    dominant_pct = 100 * dominant_prop
  ) %>%
  arrange(u_hat)

# Shared ordering
hospital_order <- hosp_summary$COD_OSPEDALE


fig_A <- ggplot(
  hosp_summary,
  aes(
    x = factor(COD_OSPEDALE, levels = hospital_order),
    y = u_hat,
    colour = dominant_cluster,
    size = dominant_pct
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_point(alpha = 0.9) +
  scale_colour_manual(values = clust_colors) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    x = NULL,
    y = expression(hat(u)[g]),
    title = "Hospital frailty"
  )



fig_B <- ggplot(sim_data,
                aes(
                  x = reorder(COD_OSPEDALE, xbeta, median),
                  y = xbeta
                )) +
  geom_boxplot(outlier.size = 0.3) +
  theme_classic()


fig_C <- hosp_cluster_prop %>%
  mutate(
    COD_OSPEDALE = factor(
      COD_OSPEDALE,
      levels = hospital_order
    )
  ) %>%
  ggplot(
    aes(
      x = COD_OSPEDALE,
      y = prop_cluster,
      fill = clusters
    )
  ) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = clust_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    x = "Hospitals",
    y = "Cluster composition",
    fill = "Cluster",
    title = "Patient cluster composition"
  )

combined_fig <- cowplot::plot_grid(
  fig_A,
  fig_B,
  fig_C,
  ncol = 1,
  align = "v",
  rel_heights = c(1.2, 0.8, 1)
)

combined_fig



# “Average patient-level risk based on observed covariates was relatively homogeneous across hospitals, whereas posterior frailty estimates showed substantial heterogeneity. This suggests that the identified clusters primarily capture latent institutional or unmeasured risk structure rather than differences in observed patient severity.”




# ── Quick sanity check ─────────────────────────────────────────────────────
sim_data %>%
  group_by(clusters) %>%
  summarise(
    n      = n(),
    mean_lp = mean(lin_pred),
    sd_lp   = sd(lin_pred)
  )


# Violin + boxplot
fig_lp <- ggplot(sim_data, aes(x = clusters, y = lin_pred, fill = clusters)) +
  geom_violin(alpha = 0.4, trim = FALSE, linewidth = 0.4) +
  geom_boxplot(width = 0.18, outlier.size = 0.8, outlier.alpha = 0.4,
               fill = "white", linewidth = 0.5) +
  scale_fill_manual(values = clust_colors) +
  labs(x = NULL,
       y = expression(italic(x)[gi]^T * hat(beta) + log ~ hat(u)[g]),
       title = "Frailty-adjusted log-hazard by cluster") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(size = 11, hjust = 0.5))

fig_lp

# ggsave("figure_linpred_by_cluster.pdf", fig_lp, width = 5, height = 4.5)
# cat("Linear predictor figure saved.\n")



fig_xbeta <- ggplot(
  sim_data,
  aes(x = clusters, y = xbeta, fill = clusters)
) +
  
  geom_violin(
    alpha = 0.4,
    trim = FALSE,
    linewidth = 0.4
  ) +
  
  geom_boxplot(
    width = 0.18,
    outlier.size = 0.8,
    outlier.alpha = 0.4,
    fill = "white",
    linewidth = 0.5
  ) +
  
  scale_fill_manual(values = clust_colors) +
  
  labs(
    x = NULL,
    y = expression(italic(x)[gi]^T * hat(beta)),
    title = expression("Observed risk component " * (X * hat(beta)))
  ) +
  
  theme_classic(base_size = 11) +
  
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11, hjust = 0.5)
  )

fig_xbeta


fig_logu <- ggplot(
  sim_data,
  aes(x = clusters, y = log_u_g, fill = clusters)
) +
  
  geom_violin(
    alpha = 0.4,
    trim = FALSE,
    linewidth = 0.4
  ) +
  
  geom_boxplot(
    width = 0.18,
    outlier.size = 0.8,
    outlier.alpha = 0.4,
    fill = "white",
    linewidth = 0.5
  ) +
  
  scale_fill_manual(values = clust_colors) +
  
  labs(
    x = NULL,
    y = expression(log(hat(u)[g])),
    title = "Latent frailty contribution by cluster"
  ) +
  
  theme_classic(base_size = 11) +
  
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11, hjust = 0.5)
  )

fig_logu




plot_grid(
  fig_xbeta,
  fig_logu,
  fig_lp,
  ncol = 3
)




# Centering and scaling the components
sim_data$xbeta_c   <- as.numeric(scale(sim_data$xbeta))
sim_data$logu_c    <- as.numeric(scale(sim_data$log_u_g))
sim_data$linpred_c <- as.numeric(scale(sim_data$lin_pred))


library(ggplot2)
library(cowplot)

# 1. Scaled Observed Risk
fig_xbeta_c <- ggplot(sim_data, aes(x = clusters, y = xbeta_c, fill = clusters)) +
  geom_violin(alpha = 0.4, trim = FALSE, linewidth = 0.4) +
  geom_boxplot(width = 0.18, outlier.size = 0.8, outlier.alpha = 0.4, fill = "white") +
  scale_fill_manual(values = clust_colors) +
  labs(x = NULL, y = "Std. units", title = expression("Scaled " * X * hat(beta))) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

# 2. Scaled Latent Frailty
fig_logu_c <- ggplot(sim_data, aes(x = clusters, y = logu_c, fill = clusters)) +
  geom_violin(alpha = 0.4, trim = FALSE, linewidth = 0.4) +
  geom_boxplot(width = 0.18, outlier.size = 0.8, outlier.alpha = 0.4, fill = "white") +
  scale_fill_manual(values = clust_colors) +
  labs(x = NULL, y = "Std. units", title = expression("Scaled " * log(hat(u)[g]))) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

# 3. Scaled Linear Predictor (Total)
fig_lp_c <- ggplot(sim_data, aes(x = clusters, y = linpred_c, fill = clusters)) +
  geom_violin(alpha = 0.4, trim = FALSE, linewidth = 0.4) +
  geom_boxplot(width = 0.18, outlier.size = 0.8, outlier.alpha = 0.4, fill = "white") +
  scale_fill_manual(values = clust_colors) +
  labs(x = NULL, y = "Std. units", title = "Scaled Total Hazard") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

# Combine
plot_grid(fig_xbeta_c, fig_logu_c, fig_lp_c, ncol = 3)






# ── 4. PATIENT SIMILARITY GRAPH ─────────────────────────────
S  <- xx$S
g_sim <- graph_from_adjacency_matrix(S, mode = "undirected",
                                     weighted = TRUE, diag = FALSE)

# Attach cluster membership as node attribute
V(g_sim)$cluster <- as.character(sim_data$clusters)

set.seed(42)
fig_graph <- ggraph(g_sim, layout = "fr") +    # Fruchterman-Reingold
  geom_edge_link(alpha = 0.015, colour = "grey60") +
  geom_node_point(aes(colour = cluster), size = 0.7, alpha = 0.85) +
  scale_colour_manual(values = clust_colors, name = "") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_void(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("figure_similarity_graph.pdf", fig_graph, width = 6, height = 5.5)
cat("Graph figure saved.\n")


# ── 5. QUICK SANITY CHECK: mortality rate per cluster ────────
cat("\n--- Observed 90-day mortality by cluster ---\n")
sim_data %>%
  group_by(clusters) %>%
  summarise(N = n(),
            Deaths = sum(cens),
            Mortality_pct = mean(cens) * 100) %>%
  print()





# ── 6. POSTERIOR FRAILTY ESTIMATES BY HOSPITAL ──────────────────────────────


hosp_cluster_comp <- sim_data %>%
  group_by(COD_OSPEDALE, clusters) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(COD_OSPEDALE) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

fig_hosp_comp <- ggplot(hosp_cluster_comp,
                        aes(x = reorder(COD_OSPEDALE, prop, FUN = max),
                            y = prop,
                            fill = clusters)) +
  geom_col() +
  scale_fill_manual(values = clust_colors) +
  labs(x = "Hospital",
       y = "Proportion of patients",
       fill = "Cluster",
       title = "Patient cluster composition within hospitals") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

fig_hosp_comp

# Safe extraction: regress out Xbeta within each hospital
# log(u_hat_g) = mean of (lin_pred - Xbeta_i) within hospital g
# But get Xbeta correctly by matching parameter names exactly

# Check exact names first
print(names(xx$estimate))
# e.g. "theta" "rho" "lambda" "ETA_AL_RICOVERO" "SESSOF" "MCS" "RESP"
#                                                  ^^ might be SESSOF not SESSO_M

# Extract beta with correct names - adjust "SESSOF" to whatever appears
beta_names <- names(xx$estimate)[-(1:3)]  # drop theta, rho, lambda
beta_est   <- xx$estimate[beta_names]
print(beta_est)  # confirm these look like regression coefficients

# Build X matching exactly the order in beta_est
# Inspect what the reference level of SESSO is:
print(levels(as.factor(sim_data$SESSO)))  # e.g. "F" "M" -> reference is "F", dummy is "M"

# Build covariate matrix to match model matrix
X <- model.matrix(~ ETA_AL_RICOVERO + SESSO + MCS + RESP, data = sim_data)
# Drop intercept column, keep only the columns matching beta_names
X_beta <- X[, beta_names, drop = FALSE]
Xbeta  <- as.numeric(X_beta %*% beta_est)

# Now extract log frailty
sim_data$log_u_hat <- sim_data$lin_pred - Xbeta

# Verify: SD within each hospital should be ~0
check <- sim_data %>%
  group_by(COD_OSPEDALE) %>%
  summarise(sd_logu = sd(log_u_hat), .groups = "drop")
print(check)  # all sd_logu should be < 1e-10



hosp_frailty <- sim_data %>%
  group_by(COD_OSPEDALE) %>%
  summarise(
    u_hat            = exp(mean(log_u_hat)),
    n_patients       = n(),
    dominant_cluster = names(which.max(table(clusters))),
    .groups          = "drop"
  ) %>%
  arrange(u_hat) %>%
  mutate(hosp_rank = factor(row_number(), levels = row_number()))

# Sanity check - should have 32 rows, no NAs
print(nrow(hosp_frailty))
print(summary(hosp_frailty$u_hat))

fig_frailty <- ggplot(hosp_frailty,
                      aes(x = hosp_rank, y = u_hat, colour = dominant_cluster)) +
  geom_hline(yintercept = 1, linetype = "dashed", 
             colour = "grey50", linewidth = 0.5) +
  geom_point(size = 3) +
  scale_colour_manual(values = clust_colors, name = "Dominant cluster") +
  labs(x     = "Hospital (ranked by frailty estimate)",
       y     = expression(hat(u)[g]),
       title = "Posterior frailty estimates by hospital") +
  theme_classic(base_size = 11) +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(size = 11, hjust = 0.5))

fig_frailty



# Variance decomposition of the linear predictor
var_total  <- var(sim_data$lin_pred)
var_Xbeta  <- var(Xbeta)
var_logu   <- var(sim_data$log_u_hat)
cov_term   <- 2 * cov(Xbeta, sim_data$log_u_hat)

cat("Var(x'beta):       ", round(var_Xbeta / var_total * 100, 1), "%\n")
cat("Var(log u_hat):    ", round(var_logu  / var_total * 100, 1), "%\n")
cat("2*Cov term:        ", round(cov_term  / var_total * 100, 1), "%\n")
cat("Total:             ", round((var_Xbeta + var_logu + cov_term) / var_total * 100, 1), "%\n")
