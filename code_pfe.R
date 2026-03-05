
# NETTOYAGE ET PACKAGES 
rm(list = ls())
setwd("C:/Users/Republic Of Computer/OneDrive - GENES/Documents/pfeENSAI/FDA")

pkgs <- c("dplyr", "readr", "lubridate", "refund", "fda", "ggplot2", "tidyr", "zoo")
invisible(lapply(pkgs, function(p) {
  if (!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))

## 1. CHARGEMENT ET PRÉPARATION DES DONNÉES
don   <- read_csv2("don.csv",   show_col_types = FALSE)
meteo <- read_csv("meteo.csv",  show_col_types = FALSE)

don$Date_plantation <- as.Date(don$Date_plantation)
meteo$DATE          <- as.Date(meteo$DATE)

var_keep <- c("BINTJE","ARKA", "DESIREE", "ESCORT", "HD R1", "HD R2", "HD R3", "HD R4")

data_reg <- don %>%
  filter(annee >= 1995, status == "Temoin", fix_asymp == FALSE, variete %in% var_keep) %>%
  mutate(variete = factor(variete, levels = var_keep)) %>%
  filter(!is.na(Date_plantation))

t_domain <- 0:120


## 2. RECONSTRUCTION ET LISSAGE DES MATRICES 

f_disease <- function(a, b, c, t) a / (1 + exp((b - t) / c))
Y_mat <- t(sapply(1:nrow(data_reg), function(i) {
  f_disease(data_reg$a[i], data_reg$b[i], data_reg$c[i], t_domain)
}))

# Fonction d'extraction avec lissage B-splines cubiques
get_weather_smoothed <- function(var_name) {
  raw_mat <- t(sapply(data_reg$Date_plantation, function(p_date) {
    target  <- p_date + days(t_domain)
    vals    <- meteo %>%
      dplyr::filter(DATE %in% target) %>%
      dplyr::select(DATE, all_of(var_name))
    df_full <- data.frame(DATE = target) %>% left_join(vals, by = "DATE")
    zoo::na.approx(df_full[[2]], na.rm = FALSE, rule = 2)
  }))
  
  # Base B-splines cubiques
  basis_obj <- create.bspline.basis(
    rangeval = range(t_domain),
    norder   = 4,
    nbasis   = 30
  )
  
  # GCV : recherche du lambda optimal 
  lambda_grid <- 10^seq(-6, 2, by = 0.5)
  
  gcv_mean <- sapply(lambda_grid, function(lam) {
    fdpar    <- fdPar(basis_obj, Lfdobj = 2, lambda = lam)
    gcv_vals <- apply(raw_mat, 1, function(y_i) {
      tryCatch(
        smooth.basis(argvals = t_domain, y = y_i, fdParobj = fdpar)$gcv,
        error = function(e) NA_real_
      )
    })
    mean(gcv_vals, na.rm = TRUE)
  })
  
  lambda_opt <- lambda_grid[which.min(gcv_mean)]
  message(sprintf("  %s : lambda optimal GCV = %.2e", var_name, lambda_opt))
  
  
  # Lissage final avec lambda optimal
  fdpar_opt <- fdPar(basis_obj, Lfdobj = 2, lambda = lambda_opt)
  fd_obj    <- smooth.basis(
    argvals  = t_domain,
    y        = t(raw_mat),
    fdParobj = fdpar_opt
  )$fd
  
  return(t(eval.fd(t_domain, fd_obj)))
}

X_list <- list(
  DH = get_weather_smoothed("DH"), RR = get_weather_smoothed("RR"),
  UM = get_weather_smoothed("UM"), TM = get_weather_smoothed("TM"), V = get_weather_smoothed("V")
)

## 3. RÉGRESSION FONCTIONNELLE 

df_fda <- list(
  Y = Y_mat, variete = data_reg$variete,
  DH = X_list$DH, RR = X_list$RR, UM = X_list$UM, TM = X_list$TM, V = X_list$V
)


# On utilise limits = "s<=t" pour l'effet historique de Boschi et al.
fit_full <- pffr(Y ~ variete + 
                   ff(DH, limits="s<=t") + 
                   ff(RR, limits="s<=t") + 
                   ff(UM, limits="s<=t") + 
                   ff(TM, limits="s<=t") + 
                   ff(V,  limits="s<=t"), 
                 data = df_fda)

## 4. CALCUL DES MÉTRIQUES
Y_hat     <- predict(fit_full)
SST       <- sum((Y_mat - mean(Y_mat))^2)
R2_global <- 1 - sum((Y_mat - Y_hat)^2) / SST

cat("\n--- METRIQUES (BOSCHI ET AL.) ---")
cat("\nR2 Global :", round(R2_global, 4), "\n")

## 5. VISUALISATION

variables_meteo <- c("DH", "RR", "UM", "TM", "V")
par(mfrow = c(2, 3), mar = c(1, 1, 3, 1))

for(v in variables_meteo) {
  coef_surf <- coef(fit_full, n1 = 50, n2 = 50)$smterms[[paste0("ff(", v, ")")]]
  if(!is.null(coef_surf)) {
    s_grid <- unique(coef_surf$x)
    t_grid <- unique(coef_surf$y)
    beta_mat <- matrix(coef_surf$value, nrow = length(s_grid), ncol = length(t_grid))
    persp(s_grid, t_grid, beta_mat,
          theta = 40, phi = 25,
          ticktype = "detailed",
          col = "lightblue", shade = 0.5,
          main = paste("Effet historique :", v),
          xlab = "Temps météo (s)",
          ylab = "Temps maladie (t)",
          zlab = "Beta(s,t)")
  }
}

par(mfrow = c(1, 1))



##### Qlq verif

# 1. Vérifier les dimensions
dim(X_list$RR)    # doit être N x 121

# 2. Vérifier que pffr a bien créé des termes ff()
summary(fit_full)  # cherchez "ff(RR)" dans les smooth terms
# La méthode la plus directe sur tes données
levels(data_reg$variete)[1]

# Ou en regardant directement comment le modèle a codé la variable
contrasts(data_reg$variete)




