# Project 1 - Analitica de Datos - UAO 2026-2
# EDA y chequeos previos (metodologia de clase, prof. Johann Ospina)
# Reproducible: Rscript R/01_eda.R desde la raiz del proyecto
#
# Datos: descomprimir primero el zip original del profesor:
#   unzip data/datos_proyecto_1.zip -d data/
suppressMessages(library(terra))
D <- "data/datos_proyecto_1/imagenes_semanales"
if (!dir.exists(D)) {
  stop("No existe ", D, "\n  Ejecuta primero: unzip data/datos_proyecto_1.zip -d data/")
}

## ---- 1. Chequeo de grilla: extension y resolucion identicas ----
r2015 <- rast(file.path(D, "chirps_semanal_valle_2015.tif"))
p2015 <- rast(file.path(D, "power_semanal_valle_2015.tif"))
stopifnot(compareGeom(r2015, p2015, res = TRUE, ext = TRUE))

## ---- 2. Limpieza de flags CHIRPS (valores negativos -9999 -> NA) ----
limpia_chirps <- function(r) { v <- values(r); v[v < 0] <- NA; values(r) <- v; r }
chirps2015 <- limpia_chirps(r2015)

## ---- 3. Mascara de pixeles validos (altitud + chirps completos) ----
alt <- rast(file.path(D, "altitud_valle.tif"))
alt[alt < 0] <- NA
mascara <- (!is.na(alt)) & (!is.na(subset(chirps2015, 1)))
celulas_validas <- which(values(mascara)[, 1])
cat("Pixeles validos en el Valle:", length(celulas_validas),
    "de", ncell(mascara), "\n")

## ---- 4. Muestreo aleatorio de puntos (indicacion de clase:
##      muestrear, NO remuestrear; covariables completas en cada punto) ----
set.seed(42)
sel <- sample(celulas_validas, 500)
xy <- xyFromCell(mascara, sel)
colnames(xy) <- c("x", "y")

## ---- 5. Extraccion de covariables en los puntos ----
prec <- extract(chirps2015, xy)[, -1]                # 52 semanas
temp <- extract(subset(p2015, 1:52), xy)[, -1]      # temp_semana_01..52
rad  <- extract(subset(p2015, 53:104), xy)[, -1]    # rad_semana_53..104
altv <- extract(alt, xy)[, 1]

completos <- complete.cases(prec, temp, rad, altv)
xy <- xy[completos, ]; prec <- prec[completos, ]; temp <- temp[completos, ]
rad <- rad[completos, ]; altv <- altv[completos]
cat("Puntos con covariables completas:", sum(completos), "/500\n")

## ---- 6. EDA ----
p_media <- rowMeans(as.matrix(prec))
cat("Precip semanal 2015 (mm): media", round(mean(p_media), 1),
    "| sd", round(sd(p_media), 1),
    "| min", round(min(as.matrix(prec)), 1),
    "| max", round(max(as.matrix(prec)), 1), "\n")
cat("Temp media (C):", round(mean(as.matrix(temp)), 1),
    "| Rad media (MJ/m2):", round(mean(as.matrix(rad)), 1),
    "| Altitud:", round(min(altv)), "-", round(max(altv)), "m\n")
cat("cor(prec, temp):", round(cor(p_media, rowMeans(as.matrix(temp))), 3),
    "| cor(prec, rad):", round(cor(p_media, rowMeans(as.matrix(rad))), 3),
    "| cor(prec, alt):", round(cor(p_media, altv), 3), "\n")

## ---- 7. Guardar dataset de puntos para modelacion ----
saveRDS(list(xy = xy, prec = prec, temp = temp, rad = rad, altitud = altv),
        "eda_puntos.rds")
cat("Guardado: eda_puntos.rds\n")
