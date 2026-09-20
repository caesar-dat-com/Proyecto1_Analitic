# Project 1 - Modelacion espacial (metodologia de clase: Script_Spatial_Analytics.R)
# Hallazgo EDA: NASA POWER tiene grilla nativa mas gruesa que CHIRPS.
#   temp completa: 123 pixeles | rad completa: 28 | interseccion: 2
# Decision documentada: modelar con precip + temp + altitud en los 123 pixeles
#   con cobertura completa (la radiacion se reporta como limitacion de cobertura,
#   exactamente el chequeo que pidio el profe en clase del 7-sep).
suppressMessages({library(terra)})
D <- "data/datos_proyecto_1/imagenes_semanales"
if (!dir.exists(D)) {
  stop("No existe ", D, "\n  Ejecuta primero: unzip data/datos_proyecto_1.zip -d data/")
}

## ---- 1. Dataset espacial: media semanal por punto, 2010-2025 ----
alt <- rast(file.path(D, "altitud_valle.tif")); alt[alt < 0] <- NA
av <- !is.na(values(alt)[, 1])

ok_t <- rep(TRUE, ncell(alt))
prec_sum <- 0; n_prec <- 0
for (y in 2010:2025) {
  rv <- values(rast(file.path(D, sprintf("chirps_semanal_valle_%d.tif", y))))
  pv <- values(rast(file.path(D, sprintf("power_semanal_valle_%d.tif", y))))
  rv[rv < 0] <- NA
  ok_t <- ok_t & (rowSums(is.na(pv[, 1:52])) == 0)
  prec_sum <- prec_sum + rowMeans(rv, na.rm = TRUE)
  n_prec <- n_prec + 1
}
cells <- which(av & ok_t)
xy <- xyFromCell(alt, cells)
df <- data.frame(
  x = xy[, 1], y = xy[, 2],
  prec = prec_sum[cells] / n_prec,
  altitud = values(alt)[cells, 1]
)
# temp media semanal (climatologia) en los mismos puntos
ct <- rast(file.path(D, "power_temp_climatologia_semanal_valle.tif"))
df$temp <- rowMeans(extract(ct, xy)[, -1], na.rm = TRUE)
cat("Dataset:", nrow(df), "puntos con precip+temp+altitud completas\n")

## ---- 2. Autocorrelacion espacial (metodo manual de clase) ----
n <- nrow(df)
dall <- as.matrix(dist(df[, c("x", "y")]))
umbral <- quantile(dall[dall > 0], 0.05)
W <- ifelse(dall < umbral & dall > 0, 1, 0)
rs <- rowSums(W); W <- W / ifelse(rs == 0, 1, rs)
zm <- mean(df$prec)
moran <- (n / sum(W)) * sum(W * outer(df$prec - zm, df$prec - zm)) / sum((df$prec - zm)^2)
geary <- ((n - 1) / (2 * sum(W))) * sum(W * (outer(df$prec, df$prec, "-"))^2) / sum((df$prec - zm)^2)
cat(sprintf("Moran I = %.4f | Geary C = %.4f (umbral vecindad %.3f grados)\n", moran, geary, umbral))

## ---- 3. GLM con enlace log (log(lambda) = X*beta, metodo de clase) ----
df$lprec <- log(df$prec)
m1 <- lm(lprec ~ temp + altitud, data = df)
print(summary(m1)$coefficients)
cat("R2 ajustado:", round(summary(m1)$adj.r.squared, 3), "\n")

## ---- 4. Validacion: LOO por punto ----
pred <- numeric(n)
for (i in 1:n) {
  mi <- lm(lprec ~ temp + altitud, data = df[-i, ])
  pred[i] <- predict(mi, df[i, ])
}
rmse_log <- sqrt(mean((df$lprec - pred)^2))
cat("RMSE (escala log, LOO):", round(rmse_log, 4), "\n")

## ---- 5. Residuos: queda autocorrelacion espacial? ----
df$res <- residuals(m1)
zr <- mean(df$res)
moran_res <- (n / sum(W)) * sum(W * outer(df$res - zr, df$res - zr)) / sum((df$res - zr)^2)
cat(sprintf("Moran I residuos = %.4f (cercano a 0 => covariables capturan la estructura)\n", moran_res))

## ---- 6. Semivariograma experimental de residuos (metodo manual de clase) ----
maxd <- max(dall[dall > 0])
bins <- seq(0.05, maxd, length.out = 20)
sv <- sapply(1:(length(bins) - 1), function(h) {
  idx <- which(dall >= bins[h] & dall < bins[h + 1] & dall > 0)
  if (length(idx) == 0) return(NA)
  ii <- idx %% n; jj <- idx %/% n + 1
  ok <- ii > 0 & jj <= n & ii < jj
  if (!any(ok)) return(NA)
  mean((df$res[ii[ok]] - df$res[jj[ok]])^2) / 2
})
semiv <- data.frame(h = head(bins, -1), gamma = sv)
print(semiv[!is.na(semiv$gamma), ])

## ---- 7. Prediccion: kriging del trend + residuo sobre la grilla ----
df$pred_trend <- predict(m1, df)
grilla <- as.data.frame(xyFromCell(alt, which(av)))
colnames(grilla) <- c("x", "y")
grilla$altitud <- values(alt)[which(av), 1]
grilla$temp <- rowMeans(extract(ct, grilla[, c("x","y")])[, -1], na.rm = TRUE)
grilla$lpred <- predict(m1, grilla)
cat("Prediccion generada para", sum(!is.na(grilla$lpred)), "pixeles validos\n")
cat("prec media predicha (mm/sem):", round(exp(mean(grilla$lpred, na.rm=TRUE)), 1),
    "| observada:", round(mean(df$prec), 1), "\n")

saveRDS(list(df = df, semiv = semiv, modelo = m1), "modelo.rds")
write.csv(df, "puntos_modelo.csv", row.names = FALSE)
cat("Guardado: modelo.rds, puntos_modelo.csv\n")
