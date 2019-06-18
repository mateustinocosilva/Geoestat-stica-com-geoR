                              ###Script para  Análises Geoestatístiscas###
                                    ###Autor: Mateus Tinôco Silva##
                                  #email: mateus-tinoco@hotmail.com
#Referências: 
#Spatial Data Management and Mapping with R: Preston Sorenson 
    #https://www.linkedin.com/pulse/spatial-data-management-mapping-r-preston-sorenson/
# The geoR package: Paulo J. Ribeiro Jr e Peter J. Diggle 
      #http://www.leg.ufpr.br/geoR/geoRdoc/geoR.pdf


# Carregando os pacotes utilizados #
library("geoR")
library("sf")
library("sp")
library("raster")
library("rgeos")
library("rgdal")

# Definindo o diretório de trabalho  #
getwd()
#setwd('diretorio')

# Carregando os dados_amostra #
dados_amostra <- read.table('Dados_Campo.csv', h=T, sep=";")
limite_interp <- st_read("Talhao_Interpolacao.shp")

# Criando o arquivo geodata #
geoadata <- as.geodata(dados_amostra, coords.col=2:3, data.col=4)

#Análise exploratoria dos dados #
summary(geoadata)
plot(geoadata)

#Criando o variograma empírico #
variog_empir=variog(geoadata, estimator.type='modulus', trend='cte')

#Plotando o variograma empírico #
plot(variog_empir)

#Ajustando o semivariograma com diferentes funções de correlação
variog.wls.linear=variofit(variog_empir,cov.model="linear",weights="cressie")
variog.wls.exp=variofit(variog_empir,cov.model="exponential",weights="cressie")
variog.wls.spherical=variofit(variog_empir,cov.model="spherical",weights="cressie")
variog.wls.gaussian=variofit(variog_empir,cov.model="gaussian",weights="cressie")
variog.wls.cubic=variofit(variog_empir,cov.model="cubic",weights="cressie")
variog.wls.circular=variofit(variog_empir,cov.model="circular",weights="cressie")


#Plotando os modelos no variograma #
plot(variog_empir)
lines(variog.wls.cubic)

#Determinando modelo com menor SSE #
lista_modelos <-list(variog.wls.spherical$cov.model, variog.wls.linear$cov.model,variog.wls.gaussian$cov.model,
                variog.wls.exp$cov.model,variog.wls.cubic$cov.model,variog.wls.circular$cov.model)
lista_sse <-list(variog.wls.exp$value,variog.wls.spherical$value,variog.wls.gaussian$value,variog.wls.linear$value,
              variog.wls.cubic$value,variog.wls.circular$value)
model_erros <- do.call(rbind,Map(data.frame,Mod=lista_modelos,SSE=lista_sse))

min_sse <- min(model_erros$SSE)
model <- model_erros[which(model_erros$SSE==min_sse),1]

#Ajustando o variograma com menor SSE #
variog_ajust <- variofit(variog_empir,cov.model=toString(model),weights="cressie")

#Plotando os modelos no variograma #
plot(variog_empir)
lines(variog_ajust)

#Ajustando o variograma com menor SSE #
valid_cruz <- xvalid(geoadata, model = variog_ajust)
plot(valid_cruz)

# Criando pontos para krigagem #
grid_krig <- st_make_grid(limite_interp, cellsize = 1, what = "centers")
grid_krig <- st_intersection(grid_krig,limite_interp)
grid_krig <- matrix(unlist(grid_krig),nrow=length(grid_krig), byrow = T)

# Realizando a krigagem #
krig_control <- krige.control(cov.model=toString(model),cov.pars=c(variog_ajust$cov.par[1], variog_ajust$cov.par[2]))
pH_plot <- krige.conv(geoadata, loc=grid_krig, krige=krig_control)

# Salvando resultados da krigagem em um spatial data frame #
pH_plot <- data.frame(cbind(pH_plot$pedict, grid_krig[,1], grid_krig[,2]))
colnames(pH_plot)=c("pH", "x", "y")
coordinates(pH_plot) <- ~x+y

# Convertendo em raster e definindo a projeção #
pH_plot <- rasterFromXYZ(pH_plot)
crs(pH_plot)<- "+init=epsg:31984"

#Plotando o resultado #
image(pH_plot, asp=1, xlab="Easting", ylab="Northing", main="pH", legend.only = T)
image.plot( legend.only = T, breaks=sl, lab.breaks=sll)
contour(pH_plot, add=T)

# Adicionando imagem da área #
rgb_image <- brick("Orto_Talhao_Interpol.tif")

# Convertendo o raster em shapefile #
pH_shape <- pH_plot

pH_shape@data@values[pH_shape@data@values <5.8  ]=1
pH_shape@data@values[pH_shape@data@values >=5.8]=2

pH_shape <- rasterToPolygons(pH_shape, dissolve=TRUE)

spplot(pH_shape,xlab="Easting", ylab="Northing", main="pH", legend.only = T,gridded=T, asp=1)

#Criando plotagem final #

jpeg(filename = "pH_plot.jpg", width = 1000, height = 727)
par(xpd=FALSE, bty="l")

pal <-colorRampPalette(c("red","yellow", "green"))
plot(pH_plot, main="pH",xlab="Easting", ylab="Northing", cex=1.5, col= pal(10) )
plotRGB(rgb_image, add=T, stretch="lin", colNA='white', main="pH",xlab="Easting (m)", ylab="Northing (m)")
plot(pH_plot, add=T, main="pH",xlab="Easting", ylab="Northing", legend=T, cex=1.5, alpha=0.75)
grid(lty=1)
par(xpd=TRUE)

contour(pH_plot, add=T, col="white")

dev.off()

# Exportando resultados #
writeRaster(pH_plot, "pH_plot.tiff", format="GTiff")
writeOGR(pH_shape, dsn= getwd(),layer = "pH_shape", driver = 'ESRI Shapefile')

