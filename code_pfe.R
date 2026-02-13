setwd("C:/Users/Republic Of Computer/OneDrive - GENES/Documents/pfeENSAI/FDA")

## ===== Régression fonctionnelle sur variétés fréquentes (présentes tous les ans) =====
## Hypothèse : courbe logistique
##   y(t) = a / (1 + exp((xmid - t)/scal))

rm(list = ls())

## Packages
pkgs <- c("dplyr", "readr", "lubridate", "fda", "ggplot2", "tidyr", "stringr")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if(length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

## 1) Charger les données 
data <- readr::read_csv2("don.csv")  

## 2) Filtrage de base 
data_clean <- data %>%
  dplyr::filter(annee >= 1995,
                status == "Temoin",
                fix_asymp == FALSE)


## 3) Garder seulement les variétés fréquentes
var_keep <- c("ARKA","BINTJE","DESIREE","ESCORT","HD R1","HD R2","HD R3","HD R4")

data_reg <- data_clean %>%
  dplyr::filter(variete %in% var_keep) %>%
  dplyr::mutate(variete = factor(variete, levels = var_keep))  # ordre contrôlé

# Vérif effectifs
print(table(data_reg$variete))

## 4) Reconstruction des courbes 
t <- seq(0, 120, 0.1)

f_disease <- function(a, xmid, scal, t){
  a / (1 + exp((xmid - t) / scal))
}

Ymat <- sapply(1:nrow(data_reg), function(i){
  f_disease(
    a    = data_reg$a[i],
    xmid = data_reg$b[i],
    scal = data_reg$c[i],
    t    = t
  )
})

## 5) Visualisation des courbes par variété 
var <- data_reg$variete
cols <- as.integer(var)

matplot(t, Ymat, type="l", col=cols, lwd=1,
        xlab="Temps (jours)",
        ylab="Sévérité maladie",
        main="Courbes reconstruites (variétés fréquentes)")

legend("topleft",
       legend=levels(var),
       col=1:length(levels(var)),
       lty=1,
       cex=0.7)

## 6) Conversion en objet fonctionnel (fd) + lissage 
range_t <- range(t)
nbasis  <- 15
basis   <- create.bspline.basis(range_t, nbasis)

Y_fd <- smooth.basis(t, Ymat, basis)$fd

plot(Y_fd, main="Courbes lissées (fd) — variétés fréquentes",
     xlab="Temps (jours)", ylab="Sévérité")

## 7) Régression fonctionnelle : Y(t) ~ variete 
# Fixer explicitement la référence (ici ARKA, car 1ère dans var_keep)
variete <- data_reg$variete

# Matrice design (intercept + dummies)
Xmat <- model.matrix(~ variete)

cat("\nDimensions Xmat:", dim(Xmat), "\n")
cat("Colonnes Xmat:", paste(colnames(Xmat), collapse=" | "), "\n")

# Base pour les coefficients beta(t)
nbasis_beta <- 10
basis_beta  <- create.bspline.basis(range_t, nbasis_beta)

# Covariables en fonctions constantes
xfdlist <- vector("list", ncol(Xmat))
const_basis <- create.constant.basis(range_t)

for(j in 1:ncol(Xmat)){
  xfdlist[[j]] <- fd(t(as.matrix(Xmat[, j])), basisobj = const_basis)
}

# Paramètres de lissage pour beta(t)
betalist <- vector("list", ncol(Xmat))
for(j in 1:ncol(Xmat)){
  betalist[[j]] <- fdPar(basis_beta, Lfdobj = 2, lambda = 1e-2)
}

# Fit
fit <- fRegress(Y_fd, xfdlist, betalist)

print(fit)

## 8) Tracer les coefficients beta(t) 
p <- length(fit$betaestlist)
nr <- ceiling(sqrt(p))
nc <- ceiling(p / nr)

par(mfrow=c(nr, nc), mar=c(4,4,2,1))

for(j in 1:p){
  plot(fit$betaestlist[[j]]$fd,
       main = colnames(Xmat)[j],
       xlab = "Temps (jours)",
       ylab = expression(beta(t)))
  abline(h = 0, lty = 2)
}

par(mfrow=c(1,1))
