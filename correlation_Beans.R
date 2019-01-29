library(raster)
# analysis Beans
###########camara micasen data
bands_M_10_09_I=read.table("clipboard", sep="\t", header = T)
id=c()
pos=which(base10$V1==6)
for(i in 1:length(base10$V1)){
  if(length(which(pos_f==i))>0){
    id=c(id,rep("c",base10$V1[i]))
  }else{id=c(id,rep(which(pos==i),base10$V1[i]))}
}
bands_M_10_09_I$ID=id
balance_data=bands_M_10_09_I[-which(bands_M_10_09_I$ID=="c"),]
##data
spect_M_10_09_I=read.table("clipboard", sep="\t", header = T)
means_reflec= function(reflectance_data_all){
  # reflectance_data_all=data_M_5_09_NI
  #################spectral data
  blue=c(465,485)
  green=c(550,570)
  red=c(663,673)
  reg=c(712,722)
  nir=c(820,880)
  bands=list(blue,green,red,reg,nir);names(bands)=c("blue","green","red","red-edge","nir")
  names(reflectance_data_all)=c("Date","Time","Treatment","Genotype","Leaf.side","Spectra",seq(325,1075,1))
reflec=list()
for (i in 1:length(bands)) {
  LW1=which(names(reflectance_data_all)==as.character(bands[[i]][1]))
  LW2=which(names(reflectance_data_all)==as.character(bands[[i]][2]))
  lw=seq(LW1,LW2,1)
  reflec[[i]]=reflectance_data_all[-c(which(reflectance_data_all$Genotype=="WR"),which(reflectance_data_all$Genotype=="BAN")),c(5,lw)]
}
names(reflec)=c("blue","green","red","red-edge","nir")
E=lapply(reflec,function(x){apply(as.matrix(x[-1]),1,mean)})
return(E)
}
M10_9=means_reflec(spect_M_10_09_I)
################Regression linear all
regre=lm(M10_9$nir~balance_data$NIR_MEAN)
correlation=cor(M10_9$nir,balance_data$NIR_MEAN)
summary(regre)
plot(balance_data$NIR_MEAN,M10_9$nir,xlab = "Camera Multiespectral",ylab = "Spectroradiometer",main = "Data Dispersion NIR")
correlation
################Regression linear per leaf
############leaf AD
AD=seq(1,280,2)
regre=lm(means_reflec$green[AD]~bands_M_5_09_NI$GREEN_MEAN[AD]) 
correlation=cor(means_reflec$blue[AD],bands_M_5_09_NI$BLUE_MEAN[AD])
summary(regre)
plot(bands_M_5_09_NI$GREEN_MEAN[AD],means_reflec$green[AD])
############leaf AB
regre=lm(means_reflec$green[-AD]~bands_M_5_09_NI$GREEN_MEAN[-AD]) 
correlation=cor(means_reflec$blue[-AD],bands_M_5_09_NI$BLUE_MEAN[-AD])
summary(regre)
plot(bands_M_5_09_NI$GREEN_MEAN[-AD],means_reflec$green[-AD])

#################per genotype
means_reflec_geno= function(reflectance_data_all){
  # reflectance_data_all=data_M_5_09_NI
  #################spectral data
  blue=c(465,485)
  green=c(550,570)
  red=c(663,673)
  reg=c(712,722)
  nir=c(820,880)
  bands=list(blue,green,red,reg,nir);names(bands)=c("blue","green","red","red-edge","nir")
  names(reflectance_data_all)=c("Date","Time","Treatment","Genotype","Leaf.side","Spectra",seq(325,1075,1))
  reflec=list()
  for (i in 1:length(bands)) {
    LW1=which(names(reflectance_data_all)==as.character(bands[[i]][1]))
    LW2=which(names(reflectance_data_all)==as.character(bands[[i]][2]))
    lw=seq(LW1,LW2,1)
    reflec[[i]]=reflectance_data_all[-c(which(reflectance_data_all$Genotype=="WR"),which(reflectance_data_all$Genotype=="BAN")),c(4,5,lw)]
  }
  names(reflec)=c("blue","green","red","red-edge","nir")
band=list()
for(i in 1:length(reflec)){
  genotype=list()
  for(j in 1:42){
    genotype[[j]]=mean(apply(reflec[[i]][which(reflec[[i]][,1]==j),-c(1:2)],1,mean))
  }
  names(genotype)=seq(1,42)
  band[[i]]=genotype
}
names(band)=c("blue","green","red","red-edge","nir")
band
}
M10_9_geno=means_reflec_geno(spect_M_10_09_I)

base10=read.table("clipboard", sep="\t", header = F)
id=c()
pos=which(base10$V1==6)
for(i in 1:length(base10$V1)){
  if(length(which(pos_f==i))>0){
    id=c(id,rep("c",base10$V1[i]))
  }else{id=c(id,rep(which(pos==i),base10$V1[i]))}
}
bands_M_10_09_I_geno=data.frame()
for(i in 1:42){
geno=cbind(t(apply(bands_M_10_09_I[which(bands_M_10_09_I$ID==as.character(i)),-11],2,mean)),ID=i)
bands_M_10_09_I_geno=rbind(bands_M_10_09_I_geno,geno)
}
##############################
regre=plot(bands_M_10_09_I_geno$BLUE_MEAN,unlist(M10_9_geno$blue))
summary(regre)
################################correlation NDVI
NDVI_CAMERA=(balance_data$NIR_MEAN-balance_data$RED_MEAN)/(balance_data$NIR_MEAN+balance_data$RED_MEAN)
NDVI_SPECTRO=(M10_9$nir-M10_9$red)/(M10_9$nir+M10_9$red)
regre=lm(NDVI_CAMERA~NDVI_SPECTRO)
correlation=cor(NDVI_CAMERA,NDVI_SPECTRO)
summary(regre)
plot(NDVI_CAMERA,NDVI_SPECTRO,xlab = "NDVI Camera",ylab = "NDVI Spectroradiometer",main = "Data Dispersion NDVI")
correlation
