library(raster)
library(gstat)
library(viridisLite)
library(automap)

#.............Preparation..............


# Get data
DataB80<- read.csv("Table 5. Spatial data B.80.csv")
coordinates(DataB80) <- ~X+Y

# grid to interpolate values
r <- raster(extent(DataB80))
res(r) <- 0.05
r.grd <- as(r, "SpatialGrid")

# shape of the sampled area
shp <- shapefile("MaskB80.shp")


#..............Spatial analysis........

#1) Total Starch Concentration automatically fitted variogram

a<-autofitVariogram(Total.Starch~1, DataB80)
plot(a)

var1 <- variogram(Total.Starch~1, DataB80, cutoff=10)
plot(var1)
fve1 <- fit.variogram(var1, vgm(108474869, "Exp", 1294, 3026))
plot(var1, fve1)
totalstB80<- plot(var1, fve1, main= "Starch concentration", cex.main=1)
plot(totalstB80)

# Krigign Total tarch concentration

krg1 <- gstat(formula=Total.Starch~1, locations=DataB80, model=fve1)
pred1 <- predict(krg1, r.grd)
pred1<- mask(raster(pred1), shp)
plot(pred1, main= "Starch concentration", cex.main=1,col= turbo (100),axes= F)

#2) Triticeae Starch variogram

var2 <- variogram(TriticeaeST~1, DataB80, cutoff=5)
plot(var2)
fve2 <- fit.variogram(var2, vgm(1.1, "Exp", 1.5, .04))
plot(var2, fve2)
tritiB80<- plot(var2, fve2, main= "Triticeae Starch", cex.main=1)
plot(tritiB80)

#Krigign Triticeae Starch

krg2 <- gstat(formula=TriticeaeST~1, locations=DataB80, model=fve2)
pred2 <- predict(krg2, r.grd)
pred2 <- mask(raster(pred2), shp)
plot(pred2, main= "Triticeae Starch", cex.main= 1,col= turbo (100),axes= F)


#3) Faboideae starch variogram

var3 <- variogram(FaboideaeST~1, DataB80, cutoff=10)
plot(var3)
fve3 <- fit.variogram(var3, vgm(0.30, "Exp", 4.1, 0.17))
plot(var3,fve3)
faboB80<-plot(var3, fve3, main= "Faboideae Starch", cex.main=1)
plot(faboB80)

# Kriging Faboideae starch

krg3 <- gstat(formula= FaboideaeST~1, locations=DataB80, model=fve3)
pred3 <- predict(krg3, r.grd)
pred3<- mask(raster(pred3), shp)
plot(pred3, main= "Faboideae Starch", cex.main= 1,col= turbo (100),axes= F)


#4) Panicoideae starch variogram

var4 <- variogram(PanicoideaeST~1, DataB80, cutoff=5)
plot(var4)
fve4<- fit.variogram(var4 , vgm(150, "Exp", 1, 0.01))
plot(var4,fve4)
panistB80<- plot(var4,fve4, main= "Panicoideae Starch ", cex.main=1)
plot(panistB80)

# Kriging Panicoideae Starch

krg4 <- gstat(formula= PanicoideaeST ~1, locations= DataB80, model=fve4)
pred4 <- predict(krg4, r.grd)
pred4<- mask(raster(pred4), shp)
plot(pred4, main= "Panicoideae Starch", cex.main= 1,col= turbo (100),axes= F)


#5) Type 5 starch variogram

var5 <- variogram(Type5~1, DataB80, cutoff=3)
plot(var5)
fve5<- fit.variogram(var5, vgm(97, "Exp", 0.90, 52))
plot(var5,fve5)
Type5<-plot(var5, fve5, main= "Type 5", cex.main=1)
plot(Type5)

# Kriging Type 5

krg5 <- gstat(formula=Type5~1, locations= DataB80, model=fve5)
pred5 <- predict(krg5, r.grd)
pred5<- mask(raster(pred5), shp)
plot(pred5, main= "Type 5", cex.main= 1,col= turbo (100), axes= F)

#6) USO starch variogram

var6 <- variogram(USO~1, DataB80, cutoff=3)
plot(var6)
fve6<- fit.variogram(var6, vgm(5, "Exp", 0.60, 2.03))
plot(var6,fve6)
USO80<-plot(var6, fve6, main= "USO Starch", cex.main=1)
plot(USO80)

# Kriging USO

krg6 <- gstat(formula=USO~1, locations=DataB80, model=fve6)
pred6 <- predict(krg6, r.grd)
pred6<- mask(raster(pred6), shp)
plot(pred6, main= "USO Starch", cex.main= 1,col= turbo (100), axes= F)

#7) USO Group A starch variogram

var7 <- variogram(USOA~1, DataB80, cutoff=3)
plot(var7)
fve7<- fit.variogram(var7, vgm(.75, "Exp", 1.1, .01))
plot(var7,fve7)
USOA<-plot(var7, fve7, main= "USO Group A", cex.main=1)
plot(USOA)

# Kriging USO Group A

krg7 <- gstat(formula=USOA~1, locations=DataB80, model=fve7)
pred7 <- predict(krg7, r.grd)
pred7<- mask(raster(pred7), shp)
plot(pred7, main= "F. USO Group A", cex.main= 1,col= turbo (100), axes= F)


#8) USO Group B starch variogram

var8 <- variogram(USOB~1, DataB80, cutoff=3)
plot(var8)
fve8<- fit.variogram(var8, vgm(.10, "Exp", 1.1, .02))
plot(var8,fve8)
USOB<-plot(var8, fve8, main= "USO Group B", cex.main=1)
plot(USOB)

# Kriging USO Group B

krg8 <- gstat(formula=USOB~1, locations= DataB80, model=fve8)
pred8 <- predict(krg8, r.grd)
pred8<- mask(raster(pred8), shp)
plot(pred8, main= "USO Group B", cex.main= 1,col= turbo (100), axes= F)

#9) USO Group C starch variogram

var9 <- variogram(USOC~1, DataB80, cutoff=15)
plot(var9)
fve9<- fit.variogram(var9, vgm(2.50, "Exp", 3.80, 1.55))
plot(var8s,fve9)
USOC<-plot(var9, fve9, main= "USO Group C", cex.main=1)
plot(USOC)

# Kriging USO Group C

krg9 <- gstat(formula=USOC~1, locations= DataB80, model=fve9)
pred9 <- predict(krg9, r.grd)
pred9<- mask(raster(pred9), shp)
plot(pred9, main= "USO Group C", cex.main= 1,col= turbo (100), axes= F)


#......Export kriging prediction as rasters.....

writeRaster(pred1, "Starch concentration 80.asc")
writeRaster(pred2, "Triticeae Starch 80.asc")
writeRaster(pred3, "Faboideae Starch 80.asc")
writeRaster(pred4, "Panicoideae Starch 80.asc")
writeRaster(pred5, "Type 5 80.asc")
writeRaster(pred6, "USO Starch 80.asc")
writeRaster(pred7, "USO Group A.asc")
writeRaster(pred8, "USO Group B.asc")
writeRaster(pred9, "USO Group C.asc")



#.................PHYTOLTIHS..............


#10) Phyto Concentration B80 variogram

b<- autofitVariogram(PxGram~1, DataB80)
plot(b)

var10 <- variogram(PxGram~1, DataB80, cutoff=17)
plot(var10)
fve10 <- fit.variogram(var10, vgm(83616470014, "Exp", 2.8, 56702186070))
plot(varP, fve10)
Pcon80<- plot(var10, fve10, main= "Phytoliths concentration", cex.main=1)
plot(Pcon80)

# Krigign Phyto concentratio B80
krg10 <- gstat(formula=PxGram~1, locations=DataB80, model=fve10)
pred10 <- predict(krg10, r.grd)
pred10<- mask(raster(pred10), shp)
plot(pred10, cex.main= 1,col= turbo (100),axes= F)


#11) Infloresence phytoltiths variogram

var11 <- variogram(Grass.inflorescence~1, DataB80, cutoff=5)
plot(var11)
fve11 <- fit.variogram(var11, vgm(650, "Exp", 1, 380))
plot(var11, fve11)
Inflo80<- plot(var11, fve11, main= "Inflorescences", cex.main=1)
plot(Inflo80)

#  Kriging Inflorescence

krg11 <- gstat(formula=Grass.inflorescence~1, locations= DataB80, model=fve11)
pred11 <- predict(krg11, r.grd)
pred11<- mask(raster(pred11), shp)
plot(pred11, main= "Grass Inflorescences", cex.main= 1,col= turbo (100),axes= F)


#12) Grass.Culm.leaf phytoltiths variogram

var12 <- variogram(Grass.Culm.Leaves~1, DataB80, cutoff=5)
plot(var12)
fve12<- fit.variogram(var12, vgm(900, "Exp", 1.40, 380))
plot(var12,fve12)
Grassculm80<- plot(var12, fve12, main= "Culms.leaves", cex.main=1)
plot(Grassculm80)


# kriging Culms.leaves

krg12 <- gstat(formula= Grass.Culm.Leaves~1, locations=DataB80, model=fve12)
pred12 <- predict(krg12, r.grd)
pred12<- mask(raster(pred12), shp)
plot(pred12, main= "Grass culms.leaves", cex.main= 1,col= turbo (100),axes= F)


#13) Panicoideae phytoltiths variogram

var13 <- variogram(Panicoideae.grass~1, DataB80, cutoff=6)
plot(var13)
fve13<- fit.variogram(var13, vgm(29, "Exp", 4.1, 19))
plot(var13,fve13)
PanicoPh80<-plot(var13, fve13, main= "Panicoideae phytoliths", cex.main=1)
plot(PanicoPh80)

#Kriging Panicoideae Phyto

krg13 <- gstat(formula=Panicoideae.grass~1, locations=DataB80, model=fve13)
pred13 <- predict(krg13, r.grd)
pred13<- mask(raster(pred13), shp)
plot(pred13, main= "Panicoideae Phytoliths", cex.main= 1,col= turbo (100), axes= F)

#14) Cyperaceae phytoltiths variogram

var14 <- variogram(Cyp.papillae~1, DataB80, cutoff=3)
plot(var14)
fve14<- fit.variogram(var14, vgm(15, "Exp", 1.8, 0.01))
plot(var14,fve14)
Cyperaceae80<-plot(var14, fve14, main= "Cyperaceae phytoliths", cex.main=1)
plot(Cyperaceae80)

#  Kriging Cyperaceae Phyto

krg14 <- gstat(formula=Cyp.papillae~1, locations= DataB80, model=fve14)
pred14 <- predict(krg14, r.grd)
pred14<- mask(raster(pred14), shp)
plot(pred14, main= "Cyperaceae phytoliths", cex.main= 1,col= turbo (100), axes= F)

#15) Arecaceae phytoltiths variogram

var15 <- variogram(Arecaceae~1, DataB80, cutoff=5)
plot(var15)
fve15<- fit.variogram(var15, vgm(0.28, "Exp", 1.5, 0.15))
plot(var15,fve15)
Palm80<-plot(var15, fve15, main= "Arecaceae", cex.main=1)
plot(Palm80)

# Kriging Arecaceae 

krg15 <- gstat(formula=Arecaceae~1, locations=DataB80, model=fve15)
pred15 <- predict(krg15, r.grd)
pred15<- mask(raster(pred15), shp)
plot(pred15, main= "Arecaceae ", cex.main= 1,col= turbo (100), axes= F)


#......Export kriging prediction as rasters.....

writeRaster(pred10, "Phytololiths concentration 80.asc")
writeRaster(pred11, "Inflorescence 80.asc")
writeRaster(pred12, "Grass.culm 80.asc")
writeRaster(pred13, "Panicoideae Phyto 80.asc")
writeRaster(pred14, "Cyperaceae 80.asc")
writeRaster(pred15, "Arecaceae 80.asc")