# More than points: Insights on the potential distribution areas and niche of the giant anteater 

# Marisa de O. Novaes, Alessandra Bertassoni, Mônica Mafra Valença Montenegro, Renata Bocorny de Azevedo, Paulo De Marco Júnior

if (!"devtools"%in%installed.packages()){install.packages("devtools")}  
devtools::install_github("andrefaa/ENMTML") 

setwd("./Dados/Capitulo 2_ModeloClima")
library(ENMTML)
library(graphics)
library(ggplot2)
library(sf)
library(terra)
library(ape)
library(spdep)
library(spatialEco)
library(openxlsx)
library(sp)
library(ncf)

# ORIGINAL ----------------------------------------------------------------

ENMTML (pred_dir = "./Preditores/PCA",
        #proj_dir=,
        occ_file = "neoxen_cur_tb.txt", 
        sp = "species",
        x = "longitude",
        y = "latitude",
        min_occ = 10,
        thin_occ = NULL,
        #thin_occ=c(method='USER-DEFINED', distance='50'),
        #part = c(method = "BOOT", replicates = '1', proportion = "0.7"),
        part = c(method='BLOCK'),
        colin_var = NULL,
        sp_accessible_area = NULL,
        pseudoabs_method = c(method= 'ENV_CONST'),
        pres_abs_ratio = 1,
        algorithm = c("GLM", "RDF", "MLK", "MXS", "GAU", "SVM"),
        thr = c(type = "SORENSEN"),
        msdm = NULL, #PRES
        ensemble = c(method= 'SUP', metric='Sorensen'))

climap <- rast(".\\Result\\Ensemble\\SUP\\Myrmecophaga_tridactyla.tif")
climapa <- rast(".\\Preditores\\Tamandua\\Projection\\Chelsa_22kBP\\Ensemble\\SUP\\Myrmecophaga_tridactyla.tif")

pontos <- read.delim("G:/Meu Drive/Doutorado_Marisa PPGBAN/Dados/Capitulo 2_ModeloClima/neoxen_cur_tb.txt")
points = st_as_sf(pontos, coords = c("longitude", "latitude")) 


# Correlation -------------------------------------------------------------

# Upload the raster
climap <- rast(".\\Preditores\\Tamandua\\Ensemble\\SUP\\Myrmecophaga_tridactyla.tif")

# Upload occurrence data (latitude and longitude)
pontos <- read.delim("G:/Meu Drive/Doutorado_Marisa PPGBAN/Dados/Capitulo 2_ModeloClima/neoxen_cur_tb.txt")

library(readxl)
occtb <- read.delim("G:/Meu Drive/Doutorado_Marisa PPGBAN/Dados/Capitulo 2_ModeloClima/Result/Occurrences_Cleaned.txt")
View(occtb)

adeqvalues=terra::extract(climap, occtb[,2:3])

cur_tb <- cbind(occtb[, 2:3], adeqvalues$layer)

View(cur_tb)

contagem_coordenadas <- table(cur_tb$x, cur_tb$y)

media_adequabilidade <- tapply(cur_tb$layer, list(cur_tb$x, cur_tb$y), mean)


points(cur_tb)

n=nrow(cur_tb)
sptb=SpatialPoints(cur_tb[,1:2])
for (i in 1:n){
  print(i)
  b=buffer(sptb[i],100000)
  coisinha=sptb %over% b
  
  cur_tb$n[i]= length(coisinha[!is.na(coisinha)])
}
View(cur_tb)

x11()
plot(cur_tb$`adeqvalues$layer` ~cur_tb$n)


# Save as a text file
write.table(cur_tb, file = "dados1.txt", sep = "\t", quote = FALSE, row.names = FALSE)

setwd("G:/Meu Drive/Doutorado_Marisa PPGBAN/Dados/Capitulo 2_ModeloClima")
pt1 = read.table("Dados1.txt", header = T)
names(pt1)
#pt1=openxlsx::read.xlsx('pt1.xlsx')
pt=SpatialPoints(pt1)
pt=SpatialPoints(pt1[,1:2])
pt2=pt1[,c(1,2,3)]

x11()
plot(pt2[,1:2],col="red")
nrow(pt2)


nrow(pt2)
names(pt2)
coordinates(pt2) <- c("x", "y")

fig1=ncf::correlog(x=pt2$x,y=pt2$y,
                   z=pt2$adeqvalues.layer,increment=25,
                   latlon=TRUE)
summary(fig1)

View(resultados)

x11()
tiff("correlation1.tif", res=600, width = 4300, height = 3200)
plot(fig1,xlim=c(0,300),ylim=c(-0,1),
     xlab="Distance (km)",
     ylab="Spatial Correlation",
     main="")
axis(side=1,at=seq(0,300,by=25))
axis(side=2,at=seq(0,1,by=0.1))
dev.off()


# Cut ---------------------------------------------------------------------
mask <- rast(".\\Dados\\Capitulo 3_ModeloDinamico\\fechado.tif")
climap <- rast(".\\Preditores\\Tamandua\\Ensemble\\SUP\\Myrmecophaga_tridactyla.tif")
br <- vect(".\\shps\\BRlimite.shp")
raster_crop <- crop(climapa, mask)
raster_res <- resample(raster_crop, mask)
raster_mask <- mask(raster_res, mask)
writeRaster(raster_mask,filename= "climapaBR.tif",overwrite=TRUE)
x11();plot(raster_mask,  col=terrain.colors(256))


# Graphs ----------------------------------------------------------------
#past
temp <- rast(".\\Preditores\\Chelsa_22kBP\\bio01.tif")
pluv <- rast(".\\Dados\\Capitulo 2_ModeloClima\\Preditores\\Chelsa_22kBP\\bio12.tif")
pontos <- read.delim("./neoxen_cur_tb.txt")
tempvalue=terra::extract(temp, pontos[,2:3])
pluvalue=terra::extract(pluv, pontos[,2:3])
tempvalue=tempvalue$`CHELSA_TraCE21k_bio01_-200_V1.0`
pluvalue=pluvalue$`CHELSA_TraCE21k_bio12_-200_V1.0`
df=data.frame(temp=temp[!is.nan(temp)],pluv=pluv[!is.nan(temp)])
names(df)=c("temp","pluv")
df=df[!is.nan(df$temp),]
df1=data.frame(tempvalue=tempvalue,pluvalue=pluvalue)

tiff("templuv22k.tif", res=600, width = 4015, height =3800, bg="transparent")
x11();ggplot() +
  geom_point(data = df, aes(x = temp, y = pluv), color = "yellow") +
  geom_point(data = df1, aes(x = tempvalue, y = pluvalue), color = "black") +
  labs(x = "Mean temperature", y = "Precipitation") +
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 15), # Ajusta o tamanho da fonte do eixo x
    axis.title.y = element_text(size = 15))
dev.off()

#present
temp <- rast(".\\Preditores\\Presente\\bio01.tif")
pluv <- rast(".\\Preditores\\Presente\\bio12.tif")
caa <- vect("./shps/caatinga1.shp")
am <- vect("./shps/amazonia.shp")
sul <- vect("./shps/sul.shp")

extract_region_data <- function(temp, pluv, region_shp) {
  temp_masked <- mask(temp, region_shp)
  pluv_masked <- mask(pluv, region_shp)
  data <- data.frame(
    temp = as.vector(values(temp_masked)),
    pluv = as.vector(values(pluv_masked))
  )
  data[!is.nan(data$temp) & !is.nan(data$pluv), ]
}

df_caatinga <- extract_region_data(temp, pluv, caa)
df_amazonia <- extract_region_data(temp, pluv, am)
df_sul <- extract_region_data(temp, pluv, sul)

df <- data.frame(temp = as.vector(values(temp)), pluv = as.vector(values(pluv)))
df <- df[!is.nan(df$temp) & !is.nan(df$pluv), ]

pontos <- read.delim("./Dados/Capitulo 2_ModeloClima/ocorrencias.txt")
tempvalue <- terra::extract(temp, pontos[, 2:3])$CHELSA_bio10_01
pluvalue <- terra::extract(pluv, pontos[, 2:3])$CHELSA_bio10_12
df1 <- data.frame(tempvalue = tempvalue, pluvalue = pluvalue)

tiff("templuv1.tif", res=600, width = 4015, height =3800, bg="transparent")
x11(); ggplot() +
  geom_point(data = df, aes(x = temp, y = pluv), color = "yellow", alpha = 0.5) +
  geom_point(data = df_caatinga, aes(x = temp, y = pluv), color = "#40E0D0", alpha = 0.8) +
  geom_point(data = df_amazonia, aes(x = temp, y = pluv), color = "#32CD32", alpha = 0.8) +
  geom_point(data = df_sul, aes(x = temp, y = pluv), color = "magenta", alpha = 0.8) +
  geom_point(data = df1, aes(x = tempvalue, y = pluvalue), color = "black") +
  geom_vline(xintercept = 15, color = "red") +
  geom_vline(xintercept = 36, color = "red") +
  labs(x = "Mean Temperature (ºC)", y = expression("Precipitation "("kg m"^-2))) +
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 15), # Ajusta o tamanho da fonte do eixo x
    axis.title.y = element_text(size = 15))
dev.off()

#temp <- rast("G:\\Meu Drive\\Doutorado_Marisa PPGBAN\\Dados\\Capitulo 2_ModeloClima\\Preditores\\Presente\\bio10.tif")
#pluv <- rast("G:\\Meu Drive\\Doutorado_Marisa PPGBAN\\Dados\\Capitulo 2_ModeloClima\\Preditores\\Presente\\bio18.tif")
br <- shapefile(".\\BR_Pais_2022.shp")
library(raster)
plot(bioma)

mtemp <- rast(".\\Preditores\\Presente\\bio10.tif")
mpluv <- rast(".\\Preditores\\Presente\\bio18.tif")
pontos <- read.delim("./neoxen_cur_tb.txt")
tempvalue=terra::extract(mtemp, pontos[,2:3])
pluvalue=terra::extract(mpluv, pontos[,2:3])
tempvalue=tempvalue$CHELSA_bio10_10
pluvalue=pluvalue$CHELSA_bio10_18
df=data.frame(mtemp=mtemp[!is.nan(mtemp)],mpluv=mpluv[!is.nan(mtemp)])
names(df)=c("mtemp","mpluv")
df=df[!is.nan(df$mtemp),]
df1=data.frame(tempvalue=tempvalue,pluvalue=pluvalue)

tiff("maximuntempluv.tif", res=600, width = 4015, height =3800, bg="transparent")
x11();ggplot() +
  geom_point(data = df, aes(x = mtemp, y = mpluv), color = "yellow") +
  geom_point(data = df1, aes(x = tempvalue, y = pluvalue), color = "black") +
  geom_vline(xintercept=15, color="red")+
  geom_vline(xintercept=36, color="red")+
  labs(x = "Mean Temperature of Warmest Quarter (ºC)", y = expression("Precipitation of Warmest Quarter " ("kg m"^-2))) +
  theme_minimal()
dev.off()

# Vicio -------------------------------------------------------------------
library(terra)

# library(rgbif)
# 
# buffer=br
# 
# # Rasterize o polígono no objeto raster
# rasterized = rasterize(buffer, raster_layer)
# plot(rasterized)
# 
# # Extrair as coordenadas dos polígonos
# coords <- lapply(slot(buffer, "polygons"), function(x) slot(x, "Polygons")[[1]]@coords)
# 
# pontos<-coords[[1]]
# wkt_pontos <- paste(pontos[,1], " ", pontos[,2], sep = "", collapse = ",")
# 
# # Crie a representação WKT final do polígono
# wkt <- paste("POLYGON((", wkt_pontos, "))", sep = "")
# 
# # Recupere dados do GBIF dentro dos limites do rasterizado
# #occurrences =occ_search(taxonKey = 359, geometry = wkt)
# occurrences =occ_search(country="BR", hasCoordinate = T)
# #occurrences =lis(occ_search(taxonKey=40102,country="BR", hasCoordinate = T))
# x11();plot(pontos)
# 
# params <- list(taxonKey=40102,
#                country=105)
#                #hasCoordinate=TRUE)
# resul <- occ_search(taxonKey=359, country="BR", limit = 60000, hasCoordinate = TRUE)$data
# View(mammal)
# nrow(resul)

library(readr)
library(sf)
library(terra)
library(dplyr)
library(rgdal)
library(tidyverse)
library(sp)
library(devtools)
library(raster)

mammal <- vect("mammal_pts.shp")
mammal <- raster("mammalsbr.tif")
mammal <- aggregate(mammal, 8)
x11();plot(mammal)
br <- raster("climapBR.tif")
#hist(br)
tiff("histBR.tif", res=600, width = 3680, height = 2865)
hist(br,
     main = "",
     xlab = "Suitability value", ylab = "Frequency",
     col = "tomato")
dev.off()

##
setEPS()  # Standard settings for EPS (resolution, size, etc.)
postscript("histBR.eps", width = 3680/600, height = 2865/600)
hist(br,
     main = "",
     xlab = "Suitability value", 
     ylab = "Frequency",
     col = "tomato")

dev.off()
climapaulo=br<0.20
plot(maskclima)

climapaulo=resample(climapaulo,mammal)     
maskclima=climapaulo>0

writeRaster(maskclima, "maskclima.tif")

caa <- vect("./shps/caatinga1.shp")
am <- vect("./shps/amazonia.shp")
sul <- vect("./shps/sul.shp")
plot(am)

mamcaa=crop(mammal,caa)
mamcaa=mask(mamcaa,caa)
writeRaster(mamcaa, "mamcaa.tif")
plot(mamcaa)
hist(mamcaa)
tiff("histcaa.tif", res=600, width = 3900, height = 2865)
hist(mamcaa,
     main = "",
     xlab = "Sampling intensity", ylab = "Frequency",
     col = "tomato")
dev.off()


mamam=crop(mammal,am)
mamam=mask(mamam,am)
writeRaster(mamam, "mamam.tif")
plot(mamam)
hist(mamam)
tiff("histam.tif", res=600, width = 3580, height = 2865)
hist(mamam,
     main = "",
     ylim=c(0,35),
     xlab = "Sampling intensity", ylab = "Frequency",
     col = "tomato")
dev.off()

mamsul=crop(mammal,sul)
mamsul=mask(mamsul,sul)
writeRaster(mamsul, "mamsul.tif")
plot(mamsul)
hist(mamsul)
tiff("histsul.tif", res=600, width = 3880, height = 2865)
hist(mamsul,
     main = "",
     ylim=c(0,35),
     xlab = "Sampling intensity", ylab = "Frequency",
     col = "tomato")
dev.off()

# Rarefaction ---------------------------------------------------------------

ENMTML (pred_dir = "./Preditores/PCA",
        #proj_dir=,
        occ_file = "neoxen_cur_tb.txt", 
        sp = "species",
        x = "longitude",
        y = "latitude",
        min_occ = 10,
        #thin_occ = NULL,
        thin_occ=c(method='USER-DEFINED', distance='30'),
        #part = c(method = "BOOT", replicates = '1', proportion = "0.7"),
        part = c(method='BLOCK'),
        colin_var = NULL,
        sp_accessible_area = NULL,
        pseudoabs_method = c(method= 'ENV_CONST'),
        pres_abs_ratio = 1,
        algorithm = c("GLM", "RDF", "MLK", "MXS", "GAU", "SVM"),
        thr = c(type = "SORENSEN"),
        msdm = NULL,
        ensemble = c(method= 'SUP', metric='Sorensen'))


# Suitability vs. density graph ----------------------------------

adebr <- rast("res200_br.tif")
denbr <- rast("mammalsbr.tif")
adebr <- rast("climapBR.tif")

vals <- cbind(values(adebr), values(denbr))
colnames(vals) <- c("suitability", "intensity")

x_mean <- mean(vals[, "suitability"], na.rm = TRUE)
y_mean <- mean(vals[, "intensity"], na.rm = TRUE)

# Categorize into quadrants
quadrant <- ifelse(vals[, "suitability"] >= x_mean & vals[, "intensity"] >= y_mean, 1,
                   ifelse(vals[, "suitability"] < x_mean & vals[, "intensity"] >= y_mean, 2,
                          ifelse(vals[, "suitability"] < x_mean & vals[, "intensity"] < y_mean, 3, 4)))

# Count the number of points in each quadrant
quadrant_counts <- table(quadrant)

colors <- c("#d7191c", "#2b83ba", "#abdda4", "#fdae61")
color_quadrant <- colors[quadrant]

tiff("graph_with_histogram.tif", res = 600, width = 8000, height = 5000)
 layout(matrix(c(1, 2), nrow = 1), widths = c(2.5, 1))  

#par(mar = c(5, 5, 2, 2))  
plot(vals[, 1], vals[, 2], 
     col = color_quadrant, 
     pch = 20, cex = 0.5, 
     xlab = "Climatic suitability", 
     ylab = "Sampling intensity",
     cex.lab = 1.4, 
     cex.axis = 1.2,)
abline(v = x_mean, col = "black", lty = 2, lwd = 2) 
abline(h = y_mean, col = "black", lty = 2, lwd = 2)

#par(mar = c(5, 1, 2, 2))  
barplot(quadrant_counts, 
        col = colors, 
        names.arg = c("Q1", "Q2", "Q3", "Q4"),
        #main = "Points per Quadrant",
        ylab = "Number of Points",
        xlab = "Quadrants",
        las = 1,
        cex.lab = 1.4,
        cex.axis = 1.2,
        cex.names = 1.4)

dev.off()

