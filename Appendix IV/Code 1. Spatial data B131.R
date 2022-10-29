library(raster)
library(gstat)
library(viridisLite)
library(automap)

#.............Preparation..............

# Get data
DataB131<- read.csv("Table 5. Spatial data B.131.csv")
coordinates(DataB131) <- ~X+Y

# grid to interpolate values to
r <- raster(extent(DataB131))
res(r) <- 0.05
r.grd <- as(r, "SpatialGrid")

# shape of the sampled area
shp <- shapefile("MaskB131.shp")


#..............Spatial analysis........

#16) Starch Concentration 131 variogram

var16 <- variogram(Total.Starch~1, DataB131, cutoff=5)
plot(var16)
fve16 <- fit.variogram(var16, vgm(139, "Exp", .50, .01))
plot(var16, fve16)
Scon131<- plot(var16, fve16, main= "Starch concentration", cex.main=1)
plot(Scon131)

# Krigign Starach concentration

krg16 <- gstat(formula=Total.Starch~1, locations=DataB131, model=fve16)
pred16 <-predict(krg16, r.grd131)
pred16<- mask(raster(pred16), shp131)
plot(pred16, main= "Starch concentration", cex.main= 1,col= turbo (100),axes= F)


#17) Triticeae starch variogram

var17 <- variogram(TriticeaeST~1, DataB131, cutoff=2.5)
plot(var17 )
fve17<- fit.variogram(var17 , vgm(100, "Exp", 1, 0))
plot(var17 ,fve17)
Triti131<- plot(var17 , fve17, main= "Triticeae starch", cex.main=1)
plot(Triti131)

# Kriging Triticeae

krg17 <- gstat(formula= TriticeaeST~1, locations= DataB131, model=fve17)
pred17 <- predict(krg17 , r.grd131)
pred17<- mask(raster(pred17), shp131)
plot(pred17, main= "Triticeae starch", cex.main= 1,col= turbo (100),axes= F)


#18) Faboideae starch IDW

gs18<- gstat(formula=FaboideaeST~1, location= DataB131,nmax = 5)
idw18<- interpolate(r131,gs18)
idw18<- mask(idw18, shp131)
plot(idw18, main= "Faboideae starch", cex.main= 1,col= turbo (100),axes= F)


#19) Panicoideae starch variogram

var19 <- variogram(PanicoideaeST~1, DataB131, cutoff=3)
plot(var19)
fve19<- fit.variogram(var19, vgm(.7, "Exp", .90, .4))
plot(var19,fve19)
Panist131<-plot(var19, fve19, main= "Panicoideae starch", cex.main=1)
plot(Panist131)

# Kriging Panicoideae

krg19 <- gstat(formula= PanicoideaeST~1, locations= DataB131, model=fve19)
pred19 <- predict(krg19, r.grd131)
pred19<- mask(raster(pred19), shp131)
plot(pred19, main= "Panicoideae starch", cex.main= 1,col= turbo (100),axes= F)

#20) USO starch variogram

var20 <- variogram(USO~1, DataB131, cutoff=2)
plot(var20)
fve20<- fit.variogram(var20, vgm(15, "Exp", 1.2, 9))
plot(var20,fve20)
USO131<-plot(var20, fve20, main= "USO starch", cex.main=1)
plot(USO131)

#Kriging USO

krg20 <- gstat(formula=USO~1, locations=DataB131, model=fve20)
pred20 <- predict(krg20 , r.grd131)
pred20 <- mask(raster(pred20 ), shp131)
plot(pred20 , main= "USO starch", cex.main= 1,col= turbo (100), axes= F)

#......Export kriging prediction as rasters.....

writeRaster(pred16, "Starch concentration 131.asc")
writeRaster(pred17, "Triticeae Starch 131.asc")
writeRaster(idw18, "Faboideae Starch 131.asc")
writeRaster(pred19, "Panicoideae Starch 131.asc")
writeRaster(pred20, "USO Starch 131.asc")




#.................PHYTOLITHS..............


# Phyto Concentration IDW 

gs21<- gstat(formula=PxGram~1, location=DataB131,nmax = 5)
idw21<- interpolate(r131,gs21)
idw21<- mask(idw21, shp131)
plot(idw21,main= "Phytoliths concentration", cex.main= 1, col= turbo (100),axes= F)

# 22) Grass inflorescence IDW

Inf22<- gstat(formula=Grass.inflorescence~1, location=DataB131,nmax = 5)
idw22<- interpolate(r131,Inf22)
idw22<- mask(idw22, shp131)
plot(idw22,main= "Inflorescence", cex.main= 1,col= turbo (100),axes= F)

#23) Grass culm/leaves  IDW

gs23<- gstat(formula=Grass.Culm.Leaves~1, location=DataB131,nmax = 5)
idw23<- interpolate(r131,gs23)
idw23<- mask(idw23, shp131)
plot(idw23,main= "Grass culm.leaves", cex.main= 1,col= turbo (100),axes= F)

#24) Panicoideae phytolith variogram

var24<- variogram(Panicoideae.grass~1, DataB131, cutoff=3)
plot(var24)
fve24<- fit.variogram(var24, vgm(2.8, "Exp", 1.1,.5))
plot(var24,fve24)
Paniph131<-plot(var24, fve24, main= "Panicoideae phytoliths", cex.main=1)
plot(Paniph131)

# Kriging Panicoideae phytoltihs

krg24 <- gstat(formula= Panicoideae.grass~1, locations=DataB131, model=fve24)
pred24 <- predict(krg24, r.grd131)
pred24<- mask(raster(pred24), shp131)
plot(pred24, main= "Panicodeae phytoliths", cex.main= 1,col= turbo (100),axes= F)

#25) Cyperaceae phytolith variogram

var25 <- variogram(Cyp.papillae~1, DataB131, cutoff=3)
plot(var25)
fve25<- fit.variogram(var25, vgm(2.5, "Exp", 1.0, 1))
plot(var25,fve25)
CyperP131<-plot(var25, fve25, main= "Cyperaceae phytoltins", cex.main=1)
plot(CyperP131)

#Kriging Cyperaceae

krg25 <- gstat(formula=Cyp.papillae~1, locations=DataB131, model=fve25)
pred25 <- predict(krg25, r.grd131)
pred25<- mask(raster(pred25), shp131)
plot(pred25, main= "Cyperaceae phytoltihs", cex.main= 1,col= turbo (100), axes= F)

#26) Arecaceae phytolith variogram

var26 <- variogram(Arecaceae~1, DataB131, cutoff=3)
plot(var26)
fve26<- fit.variogram(var26, vgm(0.14, "Exp", 4.10, 0.09))
plot(var26,fve26)
Arec131<-plot(var26, fve26, main= "Arecaceae phytoltihs", cex.main=1)
plot(Arec131)

# Kriging Palm

krg26 <- gstat(formula=Arecaceae~1, locations=DataB131, model=fve26)
pred26 <- predict(krg26, r.grd131)
pred26<- mask(raster(pred26), shp131)
plot(pred26, main= "Arecaceae phytoliths", cex.main= 1,col= turbo (100), axes= F)


#......Export kriging prediction as rasters.....

writeRaster(idw21, "Phytololiths concentration 131.asc")
writeRaster(pred22, "Inflorescence 131.asc")
writeRaster(pred23, "Grass.culm 131.asc")
writeRaster(pred24, "Panicoideae Phyto 131.asc")
writeRaster(pred25, "Cyperaceae 131.asc")
writeRaster(pred26, "Arecaceae 131.asc")