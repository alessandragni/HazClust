# ============================================================
# Section 4.3
# ============================================================
library(roxygen2)
roxygenise()

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(cowplot)
library(patchwork)
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

clust_colors <- c("Cluster 1" = "#F4D03F",   # bright yellow      (very light)
                  "Cluster 2" = "#2E86AB",   # medium blue        (medium)
                  "Cluster 3" = "#1B1B1E")   # near black         (very dark)

library(dplyr)

# ── 1. COMPUTE P-VALUES ───────────────────────────────────────
# We run the appropriate statistical test for each variable across clusters

p_age  <- anova(lm(ETA_AL_RICOVERO ~ as.factor(clusters), data = sim_data))$"Pr(>F)"[1]
p_gen  <- chisq.test(table(sim_data$SESSO, sim_data$clusters))$p.value
p_mcs  <- kruskal.test(MCS ~ as.factor(clusters), data = sim_data)$p.value
p_resp <- kruskal.test(RESP ~ as.factor(clusters), data = sim_data)$p.value
p_stat <- chisq.test(table(sim_data$cens, sim_data$clusters))$p.value
p_time <- anova(lm(time ~ as.factor(clusters), data = sim_data))$"Pr(>F)"[1]

# Helper function to format p-values nicely (<0.001 instead of 0.000)
format_p <- function(p) {
  if (p < 0.001) return("<0.001")
  else return(sprintf("%.3f", p))
}

# ── 2. YOUR ORIGINAL AGGREGATION CODE ────────────────────────
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
    Time    = sprintf("%.2f (%.2f)", Time_mean,   Time_sd)
  ) #%>%
  #select(clusters, N, Age, Male, ModMCS, Resp, Mort, Time)

# ── 3. MERGE P-VALUES INTO A FINAL FORMATTED TABLE ───────────

# Transpose the cluster summary data to match your image layout
final_table <- data.frame(
  Variable = c("Patients", "Age", "Gender", "ModMCS", "Resp", "Status", "Time"),
  Cluster_1 = c(comp_table$N[1], comp_table$Age[1], comp_table$Male[1], comp_table$ModMCS[1], comp_table$Resp[1], comp_table$Mort[1], comp_table$Time[1]),
  Cluster_2 = c(comp_table$N[2], comp_table$Age[2], comp_table$Male[2], comp_table$ModMCS[2], comp_table$Resp[2], comp_table$Mort[2], comp_table$Time[2]),
  Cluster_3 = c(comp_table$N[3], comp_table$Age[3], comp_table$Male[3], comp_table$ModMCS[3], comp_table$Resp[3], comp_table$Mort[3], comp_table$Time[3]),
  P_value   = c("--", format_p(p_age), format_p(p_gen), format_p(p_mcs), format_p(p_resp), format_p(p_stat), format_p(p_time))
)

print(final_table)




# ── 2. KAPLAN-MEIER CURVES BY CLUSTER ───────────────────────
km_fit <- survfit(Surv(time, cens) ~ clusters, data = sim_data)

fig_km <- ggsurvplot(
  km_fit,
  data          = sim_data,
  palette       = unname(clust_colors),
  linetype      = c("solid", "dashed", "dotdash"),   # B&W-safe
  linewidth     = 0.9,
  risk.table    = TRUE,
  risk.table.height = 0.28,
  risk.table.fontsize = 3.2,
  xlab          = "Time",
  ylab          = "Survival probability",
  legend.title  = "",
  legend.labs   = names(clust_colors),
  conf.int      = TRUE,
  conf.int.alpha = 0.10,
  ggtheme       = theme_classic(base_size = 15) +
    theme(legend.position = "bottom")
)


# Add log-rank p-value annotation
fig_km$plot <- fig_km$plot +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Log-rank p ", 
                          format.pval(surv_pvalue(km_fit, sim_data)$pval, 
                                      digits = 2, eps = 0.001)),
           hjust = 1.05, vjust = 1.5, size = 6, fontface = "italic")

# Save
# pdf("figure_km_by_cluster.pdf", width = 7, height = 6)
print(fig_km)
# dev.off()
cat("KM figure saved.\n")


ggsave(
  filename = "casestudy/result/KaplanMeier.pdf",
  plot     = fig_km$plot,
  width    = 8,
  height   = 6,
  device   = "pdf"
)



# ── 3. FRAILTY-ADJUSTED LINEAR PREDICTOR BY CLUSTER ─────────

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





# ── Hospital-level decomposition ------------------------------------------


hosp_summary <- sim_data %>%
  group_by(COD_OSPEDALE) %>%
  summarise(
    mean_xbeta = mean(xbeta),
    mean_log_u = mean(log_u_g),
    u_hat      = exp(mean(log_u_g)),
    mean_eta   = mean(lin_pred),
    n_patients = n(),
    .groups = "drop"
  )

hosp_cluster_comp <- sim_data %>%
  group_by(COD_OSPEDALE, clusters) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(COD_OSPEDALE) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# dominant cluster per hospital
dominant_cluster <- hosp_cluster_comp %>%
  group_by(COD_OSPEDALE) %>%
  slice_max(prop, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    COD_OSPEDALE,
    dominant_cluster = clusters,
    dominant_prop = prop
  )

# merge into hospital summary
hosp_summary <- hosp_summary %>%
  left_join(dominant_cluster, by = "COD_OSPEDALE") %>%
  mutate(
    dominant_pct = 100 * dominant_prop
  ) %>%
  arrange(u_hat)

# shared ordering
hospital_order <- hosp_summary$COD_OSPEDALE




fig_A <- ggplot(
  hosp_summary,
  aes(x = factor(COD_OSPEDALE, levels = hospital_order),
      y = u_hat,
      colour = dominant_cluster,
      size = dominant_pct)
) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.9) +
  scale_colour_manual(values = clust_colors, name = NULL) +
  scale_size_continuous(
    range  = c(1.5, 6),
    name   = "%",
    breaks = c(50, 75, 100)
  ) +
  guides(
    colour = guide_legend(order = 1, override.aes = list(size = 3)),
    size   = guide_legend(order = 2)
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.box = "horizontal"
  ) +
  labs(
    x = NULL,
    y = expression(hat(u)[g]) #,
    #title = "Posterior frailty and cluster composition by hospital"
  )

fig_C <- hosp_cluster_comp %>%
  mutate(COD_OSPEDALE = factor(COD_OSPEDALE, levels = hospital_order)) %>%
  ggplot(aes(x = COD_OSPEDALE, y = prop, fill = clusters)) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = clust_colors, name = NULL) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_classic(base_size = 15) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  labs(
    x = "Hospitals (ranked by increasing posterior frailty)",
    y = "Proportion of patients",
    title = NULL
  )

combined_fig <- (fig_A / fig_C) +
  plot_layout(heights = c(1.2, 1))

combined_fig

ggsave(
  filename = "casestudy/result/frailty.pdf",
  plot     = combined_fig,
  width    = 7,
  height   = 7,
  device   = "pdf"
)





# ── Cluster-specific survival curves ────────────────────────────────────────
# Parameters are stored on log scale (transform=TRUE), so exponentiate first
lambda_hat <- exp(xx$estimate[which(names(xx$estimate) == "lambda")])
rho_hat    <- exp(xx$estimate[which(names(xx$estimate) == "rho")])

# Time grid up to 99th percentile of observed times
times <- seq(0, quantile(sim_data$time, 0.99), length.out = 300)

surv_cluster <- bind_rows(lapply(levels(sim_data$clusters), function(cl) {
  
  idx <- sim_data$clusters == cl
  eta <- sim_data$lin_pred[idx]   # x'β + log(û_g) — already computed
  
  # S(t | x_i, u_{g(i)}) averaged over patients in this cluster
  # outer() gives an [n_cl × 300] matrix, colMeans marginalises over patients
  S_mat <- outer(exp(eta), times,
                 function(e, t) exp(-lambda_hat * t^rho_hat * e))
  
  data.frame(
    time    = times,
    surv    = colMeans(S_mat),
    cluster = cl
  )
}))

# ── Plot ─────────────────────────────────────────────────────────────────────
ggplot(surv_cluster, aes(x = time, y = surv, colour = cluster, linetype = cluster)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = clust_colors, name = NULL) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash"), name = NULL) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  labs(
    x = "Time",
    y = "Survival probability",
    title = "Frailty-adjusted survival by cluster"
  ) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom")



ggplot(surv_cluster, aes(x = time, y = surv, 
                         colour = cluster, linetype = cluster)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = clust_colors, name = NULL) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash"), name = NULL) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  labs(
    x = "Time",
    y = "Survival probability",
    title = "Frailty-adjusted survival by cluster"
  ) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom")



# Build individual survival curves for every patient
indiv_surv <- bind_rows(lapply(seq_len(nrow(sim_data)), function(i) {
  eta_i <- sim_data$lin_pred[i]
  data.frame(
    time    = times,
    surv    = exp(-lambda_hat * times^rho_hat * exp(eta_i)),
    cluster = sim_data$clusters[i],
    patient = i
  )
}))




lastfig = ggplot() +
  # Individual curves
  geom_line(data = indiv_surv, 
            aes(x = time, y = surv, group = patient, colour = cluster),
            alpha = 0.05, linewidth = 0.2) +
  # Cluster averages on top
  geom_line(data = surv_cluster,
            aes(x = time, y = surv, colour = cluster, linetype = cluster),
            linewidth = 1.2) +
  scale_colour_manual(values = clust_colors, name = NULL) +
  scale_fill_manual(values = clust_colors, name = NULL) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash"), name = NULL) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  labs(
    x = "Time",
    y = "Survival probability" #,
    #title = "Individual and average frailty-adjusted survival by cluster"
  ) +
  theme_classic(base_size = 15) +
  theme(legend.position = "top")


lastfig

ggsave(
  filename = "casestudy/result/Survivalprob.pdf",
  plot     = lastfig,
  width    = 7,
  height   = 7,
  device   = "pdf"
)







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
  geom_violin(alpha = 0.8, trim = FALSE, linewidth = 0.4) +
  geom_boxplot(width = 0.18, outlier.size = 0.8, outlier.alpha = 0.4,
               fill = "white", linewidth = 0.5) +
  scale_fill_manual(values = clust_colors) +
  labs(x = NULL,
       y = expression(italic(x)[gi]^T * hat(beta) + log ~ hat(u)[g]),
       title = expression(italic(x)[gi]^T * hat(beta) + log ~ hat(u)[g])) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none",
        plot.title = element_text(size = 17, hjust = 0.5))

fig_lp

ggsave(
  filename = "casestudy/result/fig_lp.pdf",
  plot     = fig_lp,
  width    = 7,
  height   = 7,
  device   = "pdf"
)

# ggsave("figure_linpred_by_cluster.pdf", fig_lp, width = 5, height = 4.5)
# cat("Linear predictor figure saved.\n")



fig_xbeta <- ggplot(
  sim_data,
  aes(x = clusters, y = xbeta, fill = clusters)
) +
  
  geom_violin(
    alpha = 0.8,
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
    title = expression(italic(x)[gi]^T * hat(beta))
  ) +
  
  theme_classic(base_size = 15) +
  
  theme(
    legend.position = "none",
    plot.title = element_text(size = 17, hjust = 0.5)
  )

fig_xbeta

ggsave(
  filename = "casestudy/result/fig_xbeta.pdf",
  plot     = fig_xbeta,
  width    = 7,
  height   = 7,
  device   = "pdf"
)


fig_logu <- ggplot(
  sim_data,
  aes(x = clusters, y = log_u_g, fill = clusters)
) +
  
  geom_violin(
    alpha = 0.8,
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
  
  theme_classic(base_size = 15) +
  
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11, hjust = 0.5)
  )

fig_logu




XX = plot_grid(
  #fig_logu,
  fig_lp,
  fig_xbeta,
  ncol = 2 #3
)

XX


ggsave(
  filename = "casestudy/result/LogHazard.pdf",
  plot     = XX,
  width    = 11,
  height   = 7,
  device   = "pdf"
)



