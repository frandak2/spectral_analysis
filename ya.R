if(!require(ggplot2)){install.packages('ggplot2');library(ggplot2)}else{library(ggplot2)}
if(!require(stringr)){install.packages('stringr');library(stringr)}else{library(stringr)}
if(!require(plotly)){install.packages('plotly');library(plotly)}else{library(plotly)}
if(!require(utils)){install.packages('utils');library(utils)}else{library(utils)}
if(!require(signal)){install.packages('signal');library(signal)}else{library(signal)}
if(!require(reshape2)){install.packages('reshape2');library(reshape2)}else{library(reshape2)}
if(!require(scales)){install.packages('scales');library(scales)}else{library(scales)}
if(!require(MASS)){install.packages('MASS'); library(MASS)} else {library(MASS)}
if(!require(reshape2)){install.packages('reshape2'); library(reshape2)} else {library(reshape2)}
if(!require(car)){install.packages('car'); library(car)} else {library(car)}
if(!require(nortest)){install.packages('nortest'); library(nortest)} else {library(nortest)}
if(!require(PMCMR)){install.packages('PMCMR'); library(PMCMR)} else {library(PMCMR)}
if(!require(DescTools)){install.packages('DescTools'); library(DescTools)} else {library(DescTools)}
#########################lectura de firmas espectrales
lee_firmas=function(dir_signature,Treatment,Levels,Replicates,Date_start,extension,decimal,Separations,delte_line,Wavelength_start,Wavelength_end,Col_Wavelength,Col_Reflectance){
  setwd(dir_signature)
  Treatment=as.character(Treatment)
  Levels=as.character(Levels)
  foldersF=list.dirs(getwd(),full.names = T)
  if(!missing(Replicates)){
    Replicates=as.character(Replicates)
    e=list()
    for(i in 1:length(Treatment)){
       # i=2
      foldersT=list.dirs(foldersF,full.names = T)
      posT=str_which(foldersF,Treatment[i])
      matrix.fir=list()
      for(j in 1:length(Levels)){
        # j=1
        x=paste0("/",Levels[j])
        foldersN=foldersF[posT]
        posN=str_which(foldersN,x)
        if(length(posN)>0){
          repli=lapply(Replicates,function(y){
             # y=Replicates[1]
            g=paste0("/",y)
            foldersR=foldersN[posN]
            posR=str_which(foldersR,g)
            if(length(posR)>0){
              r=foldersR[posR]
              fir_1=str_subset(list.files(r,full.names = T),extension)
              longi=as.data.frame(read.delim(fir_1[1],header = F,sep = Separations,dec = decimal,skip =delte_line))[Col_Wavelength]
              lon=c(which(longi==Wavelength_start):which(longi==Wavelength_end))
              fir_1=t(sapply(fir_1, function(k){as.data.frame(read.delim(k,header = F,sep = Separations,dec = decimal,skip =delte_line))[lon,Col_Reflectance]}))
              dates=na.omit(unique(str_locate(row.names(fir_1),Date_start)))
              dates=str_sub(row.names(fir_1),dates[1,1],dates[1,2])
              len=max(unlist(str_locate_all(row.names(fir_1)[1],"/")))
              id=paste0(dates,"_",Treatment[i],"_",Levels[j],"_",y,"_",str_remove(str_sub(row.names(fir_1),len+1,len+30),extension))
              row.names(fir_1)=id;colnames(fir_1, do.NULL = FALSE);colnames(fir_1)=longi[lon,]
              if(max(fir_1)>=40){max_wave=100}else{max_wave=1}
              outliers=ifelse(fir_1>=0 & fir_1<=max_wave,0,1);outliers=apply(outliers,1,sum)
              if(length(which(outliers>0))==0){
                fir_1
              }else{
                fir_1=fir_1[-which(outliers>0),]
              }
              # fix(fir_1)
            }
          })
          names(repli)=Replicates
          co=c()
          for (q in 1:length(Replicates)) {
            if(is.null(repli[[q]])){
              co=c(co,j)
            }
          }
          if(length(co)>0){
            repli=repli[-co]
          }
          
        }
        matrix.fir[[j]]=repli
      }
      names(matrix.fir)=Levels
      print(paste0("Detected Replicates-",i))
      co=c()
      for (q in 1:length(Levels)) {
        if(is.null(matrix.fir[[q]])){
          co=c(co,j)
        }
      }
      if(length(co)>0){
        matrix.fir=matrix.fir[-co]
      }
      e[[i]]=matrix.fir
    }
    names(e)=Treatment
    dir=paste0(getwd(),"/Matrices_Firma")
    if(!file.exists(dir)){
      dir.create(dir)
    }
    dir=paste0(dir,"/",names(e))
    if(!file.exists(dir)){
      lapply(dir,function(x){dir.create(x)})
    }
    for (q in 1:length(Treatment)) {
      for (w in 1:length(Levels)) {
        for (r in 1:length(Replicates)) {
          dir1=paste0(dir[q],"/",names(e)[q],"_",names(e[[q]])[w],"_",names(e[[q]][[w]][r]),".csv")
          db=e[[q]][[w]][r]
          write.csv(db,dir1,row.names = TRUE)
        }
      }
    }
  }else{e=list()
  for(i in 1:length(Treatment)){
    # i=2
    # foldersT=list.dirs(foldersF,full.names = T)
    posT=which(str_extract(foldersF,Treatment[i])==Treatment[i])
    matrix.fir=lapply(Levels,function(x){
      # x=Levels[2]
      foldersN=foldersF[posT]
      posN=which(str_extract(foldersN,x)==x)
      if(length(posN)>0){
        r=foldersN[posN]
        fir_1=str_subset(list.files(r,full.names = T),extension)
        longi=as.data.frame(read.delim(fir_1[1],header = F,sep = Separations,dec = decimal,skip =delte_line))[Col_Wavelength]
        lon=c(which(longi==Wavelength_start):which(longi==Wavelength_end))
        fir_1=t(sapply(fir_1, function(t){as.data.frame(read.delim(t,header = F,sep = Separations,dec = decimal,skip =delte_line))[lon,Col_Reflectance]}))
        dates=na.omit(unique(str_locate(row.names(fir_1),Date_start)))
        dates=str_sub(row.names(fir_1),dates[1,1],dates[1,2])
        len=max(unlist(str_locate_all(row.names(fir_1)[1],"/")))
        id=paste0(dates,"_",Treatment[i],"_",x,"_",str_remove(str_sub(row.names(fir_1),len+1,len+30),extension))
        row.names(fir_1)=id;colnames(fir_1, do.NULL = FALSE);colnames(fir_1)=longi[lon,]
        if(max(fir_1)>=40){max_wave=100}else{max_wave=1}
        outliers=ifelse(fir_1>=0 & fir_1<=max_wave,0,1);outliers=apply(outliers,1,sum)
        if(length(which(outliers>0))==0){
          fir_1
        }else{
          fir_1=fir_1[-which(outliers>0),]
        }
        fir_1
      }
      
    })
    names(matrix.fir)=Levels
    co=c()
    for (j in 1:length(Levels)) {
      if(is.null(matrix.fir[[j]])){
        co=c(co,j)
      }
    }
    print(paste0("no detected Replicates-",i))
    if(length(co)>0){
      matrix.fir=matrix.fir[-co]
    }
    e[[i]]=matrix.fir
  }
  names(e)=Treatment
  dir=paste0(getwd(),"/Matrices_Firma")
  if(!file.exists(dir)){
    dir.create(dir)
  }
  dir=paste0(dir,"/",names(e))
  if(!file.exists(dir)){
    lapply(dir,function(x){dir.create(x)})
  }
  for (i in 1:length(Treatment)) {
    dir1=paste0(dir[i],"/",names(e)[i],"_",names(e[[i]]),".csv")
    for (j in 1:length(dir1)) {
      db=e[[i]][[j]]
      write.csv(db,dir1[j],row.names = TRUE)
      
    }
  }
  }
  return(e)
  print("Done and saved")
}
###############Promedio firmas espectrales
stats_fir=function(.fir,Replicate=F){
  # .fir=firmasSpectrales;Replicate=T
  if(Replicate){
    G=list()
    for(m in 1:length(.fir)){
      # m=1
      x=.fir[[m]]
      for (h in 1:length(x)){
        # h=1
        d=x[[h]]
        for (k in 1:length(d)) {
          mean=apply(d[[k]],2,mean)
          sd=apply(d[[k]],2,sd)
          cv=(sd/mean)*100
          db=rbind(mean,sd,cv)
          d[[k]]=db
        }
        x[[h]]=d
      }
      G[[m]]=x
    }
    names(G)=names(.fir)
    dir=paste0(getwd(),"/Matrices_Firma")
    dir=list.dirs(dir,full.names = T)[-1]
    for (q in 1:length(G)) {
      for (w in 1:length(G[[q]])) {
        for (r in 1:length(G[[q]][[w]])) {
          dir1=paste0(dir[q],"/",names(G)[q],"_",names(G[[q]])[w],"_",names(G[[q]][[w]][r]),".csv")
          db=G[[q]][[w]][r]
          write.csv(db,dir1,row.names = TRUE)
        }
      }
    }
  }else{
  G=list()
  for(m in 1:length(.fir)){
    # m=1
    x=.fir[[m]]
    for (h in 1:length(x)){
      # h=1
      d=x[[h]]
      mean=apply(d,2,mean)
      sd=apply(d,2,sd)
      cv=(sd/mean)*100
      db=rbind(mean,sd,cv)
      x[[h]]=db
    }
    G[[m]]=x
  }
  names(G)=names(.fir)
  dir=paste0(getwd(),"/Matrices_Firma")
  dir=list.dirs(dir,full.names = T)[-1]
  for (i in 1:length(dir)) {
    # i=1
    dir1=paste0(dir[i],"/",names(G)[i],"_",names(G[[i]]),"_stats.csv")
    for (j in 1:length(dir1)) {
      # j=1
      db=G[[i]][[j]]
      write.csv(db,dir1[j])
    }
  }}
  return(G)
  print("Done and saved")
}
#######################calculo indices de vegetacion
cal_IV=function(.fir,Replicate=F){
  # .fir=firmasSpectrales;Replicate=T
  if(Replicate){
    G=list()
    for(m in 1:length(.fir)){
      # m=1
      x=.fir[[m]]
      for (h in 1:length(x)){
        # h=1
        d=x[[h]]
        for (t in 1:length(d)) {
          p=melt(d[[t]]);p$id=str_sub(p$Var1,12,nchar(as.character(p$Var1[1])));p$Var1=str_sub(p$Var1,1,10);names(p)=c("fecha","wave","reflec","id");p$fecha=as.Date(p$fecha)
          ######W1
          r670=p$reflec[which(p$wave==670)]
          r700=p$reflec[which(p$wave==700)]
          r740=p$reflec[which(p$wave==740)]
          r780=p$reflec[which(p$wave==780)]
          re=(r670+r780)/2
          REP_li=700+40*((re-r700)/(r740-r700))
          
          #normalised difference red edge
          # NDRE=(790-720)/(790+720)
          r790=p$reflec[which(p$wave==790)]
          r720=p$reflec[which(p$wave==720)]
          NDRE=(r790-r720)/(r790+r720)
          ##NDVI705 Indice modificado Red Edge Ratio simple
          r750=p$reflec[which(p$wave==750)]
          r705=p$reflec[which(p$wave==705)]
          NDVI705=(r750-r705)/(r750+r705)
          # Chl_I, chlorophyll indep
          r740=p$reflec[which(p$wave==740)]
          r720=p$reflec[which(p$wave==720)]
          chl_I=r740/r720
          # Indice modificado Red Edge Ratio simple
          r445=p$reflec[which(p$wave==445)]
          mSR=(r705-r445)/(r705+r445)
          #MCARI Indice modificado de relaciI?n en absorciI?n de clorofila
          r550=p$reflec[which(p$wave==500)]
          MCARI=((r700-r670)-0.2*(r700-r550))*(r700/r670)
          # TCARI Indice de reflectancia de absorciOn transformada de clorofila
          TCARI=3*(((r700-r670)-0.2*(r700-r550))*(r700/r670))
          #OSAVI Indice optimizado de vegetaciOn del suelo ajustado
          r800=p$reflec[which(p$wave==800)]
          OSAVI=1.16*((r800-r670)/(r800+r670+0.16))
          #division entre TCARI Y OSAVI
          div=TCARI/OSAVI
          #formula clorofilometro
          r940=p$reflec[which(p$wave==940)]
          r650=p$reflec[which(p$wave==650)]
          Spad=log(r940/r650)
          f=p$fecha[which(p$wave==650)]
          id=p$id[which(p$wave==650)]
          ####correlaciones entre clorofila y indices de vegetacion
          df=data.frame(fecha=f,id=id,REP=REP_li,NDVI=NDVI705,NDRE=NDRE,CHL=chl_I,OSAVI=OSAVI,
                        MSR=mSR,MCARI=MCARI,TCARI=TCARI,DIV=div,SPAD=Spad)
          d[[t]]=df
        }
        x[[h]]=d
      }
      G[[m]]=x
    }
    names(G)=names(.fir)
    dir=paste0(getwd(),"/Matrices_Firma")
    dir=list.dirs(dir,full.names = T)[-1]
    for (q in 1:length(G)) {
      for (w in 1:length(G[[q]])) {
        for (r in 1:length(G[[q]][[w]])) {
          dir1=paste0(dir[q],"/",names(G)[q],"_",names(G[[q]])[w],"_",names(G[[q]][[w]][r]),".csv")
          db=G[[q]][[w]][r]
          write.csv(db,dir1,row.names = TRUE)
        }
      }
    }
  }else{G=list()
  for(m in 1:length(.fir)){
    # m=1
    x=.fir[[m]]
    for (h in 1:length(x)){
      # h=1
      d=x[[h]]
      p=melt(d);p$id=str_sub(p$Var1,12,nchar(as.character(p$Var1[1])));p$Var1=str_sub(p$Var1,1,10);names(p)=c("fecha","wave","reflec","id");p$fecha=as.Date(p$fecha)
      ######W1
      r670=p$reflec[which(p$wave==670)]
      r700=p$reflec[which(p$wave==700)]
      r740=p$reflec[which(p$wave==740)]
      r780=p$reflec[which(p$wave==780)]
      re=(r670+r780)/2
      REP_li=700+40*((re-r700)/(r740-r700))
      
      #normalised difference red edge
      # NDRE=(790-720)/(790+720)
      r790=p$reflec[which(p$wave==790)]
      r720=p$reflec[which(p$wave==720)]
      NDRE=(r790-r720)/(r790+r720)
      ##NDVI705 Indice modificado Red Edge Ratio simple
      r750=p$reflec[which(p$wave==750)]
      r705=p$reflec[which(p$wave==705)]
      NDVI705=(r750-r705)/(r750+r705)
      # Chl_I, chlorophyll indep
      r740=p$reflec[which(p$wave==740)]
      r720=p$reflec[which(p$wave==720)]
      chl_I=r740/r720
      # Indice modificado Red Edge Ratio simple
      r445=p$reflec[which(p$wave==445)]
      mSR=(r705-r445)/(r705+r445)
      #MCARI Indice modificado de relaciI?n en absorciI?n de clorofila
      r550=p$reflec[which(p$wave==500)]
      MCARI=((r700-r670)-0.2*(r700-r550))*(r700/r670)
      # TCARI Indice de reflectancia de absorciOn transformada de clorofila
      TCARI=3*(((r700-r670)-0.2*(r700-r550))*(r700/r670))
      #OSAVI Indice optimizado de vegetaciOn del suelo ajustado
      r800=p$reflec[which(p$wave==800)]
      OSAVI=1.16*((r800-r670)/(r800+r670+0.16))
      #division entre TCARI Y OSAVI
      div=TCARI/OSAVI
      #formula clorofilometro
      r940=p$reflec[which(p$wave==940)]
      r650=p$reflec[which(p$wave==650)]
      Spad=log(r940/r650)
      f=p$fecha[which(p$wave==650)]
      id=p$id[which(p$wave==650)]
      ####correlaciones entre clorofila y indices de vegetacion
      df=data.frame(fecha=f,id=id,REP=REP_li,NDVI=NDVI705,NDRE=NDRE,CHL=chl_I,OSAVI=OSAVI,
                    MSR=mSR,MCARI=MCARI,TCARI=TCARI,DIV=div,SPAD=Spad)
      x[[h]]=df
    }
    G[[m]]=x
  }
  names(G)=names(.fir)
  dir=paste0(getwd(),"/Matrices_Firma")
  dir=list.dirs(dir,full.names = T)[-1]
  for (i in 1:length(dir)) {
    # i=1
    dir1=paste0(dir[i],"/",names(G)[i],"_",names(G[[i]]),"_IV.csv")
    for (j in 1:length(dir1)) {
      # j=1
      db=G[[i]][[j]]
      write.csv(db,dir1[j])
      
    }
  }}
  return(G)
  print("Done and saved")# assign(paste0(unique(str_sub(row.names(.fir[[1]][[1]]),12,12)),".IV's"),G,envir = .GlobalEnv)
}
############PLOT iv
plot_IV=function(.IV,nameIV){
  .IV=IV
  nameIV="NDVI"
  p=list()
  for (j in 1:length(.cor)) {
    j=1
    x=.IV[[j]]
    all=data.frame()
    for (i in 1:length(x)) {
      i=1
      a=x[[i]]
      all=rbind(all,a)
    }
    db=cbind(fecha=all$fecha,IV=all[which(names(all)==nameIV)],variable=names(x))
    plot=ggplot(db, aes(fecha,IV))+ geom_point(aes(color=Level),size=0.72)+
      theme_bw()+
      labs(title = paste0("Correlacion entre ",unique(all$variable)," y Refelctancia"))+
      labs(x = "Fecha", y = "IV")+
      geom_hline(yintercept =0,linetype="solid",color = "red")+
      theme(legend.text=element_text(size=10),plot.title =element_text(hjust=0.5),plot.subtitle=element_text(hjust=0.5))+
      scale_colour_manual(values = 1:length(unique(all$Level)))
    # name="C:/Users/frand/Google Drive/la U/ultimo semestre/tg2/presentacion/analisis/indices/CORlw9-in.png"
    # png(filename=name,width = 5260, height = 469)
    p[[j]]=plot
  }
  p
}
#########################plot stats###############
plots_stats=function(stats_fir){
##########create dir
dir=paste0(getwd(),"/Plots")
if(!file.exists(dir)){
  dir.create(dir)
}
dir1=paste0(dir,"/Treatment")
if(!file.exists(dir1)){
  dir.create(dir1)
}
dir2=paste0(dir,"/Levels")
if(!file.exists(dir2)){
  dir.create(dir2)
}
stats_fir=melt(stats_fir);names(stats_fir)=c("stats","wave","reflec","Levels","Treatment")
for (i in 1:length(unique(stats_fir$Treatment))) {
  for (j in 1:3) {
    Tra=stats_fir[which(stats_fir$Treatment==unique(stats_fir$Treatment)[i]),]
    mean=Tra[which(Tra$stats==unique(Tra$stats)[j]),]
    name=paste0(dir1,"/",unique(Tra$stats)[j],"_",unique(stats_fir$Treatment)[i],".png")
    png(filename=name,width = 726, height = 469)
    ###########plot Treatment
    p=ggplot(mean, aes(wave,reflec))+ geom_line(aes(color=Levels),size=0.72)+
      theme_bw()+
      geom_vline(xintercept = c(450,510,530,590,620,670,740,820),
                 linetype="dashed",color = "black")+
      annotate("text", x = c(480,560,645,780,700), y = max(mean$reflec)/2,
               label = c("B","G","R","NIR","RE"))+
      theme(legend.text=element_text(size=10))
    if(j==1){
      p=p+labs(title = paste0("Mean Feature Spectral"),
               subtitle = paste0("Levels ",unique(Tra$Treatment)),
               x = "Wavelength", y = "Reflectance (%)") 
    }
    if(j==2){
      p=p+labs(title = paste0("SD Feature Spectral"),
               subtitle = paste0("Levels ",unique(Tra$Treatment)),
               x = "Wavelength", y = "Standard Deviation") 
    }
    if(j==3){
      p=p+labs(title = paste0("CV Feature Spectral"),
               subtitle = paste0("Levels ",unique(Tra$Treatment)),
               x = "Wavelength", y = "Coefficient Variation (%)") 
    }
    print(p)
    dev.off() 
  }
}
####################plot Levels
for (i in 1:length(unique(stats_fir$Treatment))) {
  for (j in 1:3) {
    niv=stats_fir[which(stats_fir$Levels==unique(stats_fir$Levels)[i]),]
    mean=niv[which(niv$stats==unique(niv$stats)[j]),]
    name=paste0(dir2,"/",unique(niv$stats)[j],"_",unique(stats_fir$Levels)[i],".png")
    png(filename=name,width = 726, height = 469)
    p=ggplot(mean, aes(wave,reflec))+ geom_line(aes(color=Treatment),size=0.72)+
      theme_bw()+
      geom_vline(xintercept = c(450,510,530,590,620,670,740,820),
                 linetype="dashed",color = "black")+
      annotate("text", x = c(480,560,645,780,700), y = max(mean$reflec)/2,
               label = c("B","G","R","NIR","RE"))+
      theme(legend.text=element_text(size=10))
    if(j==1){
      p=p+labs(title = paste0("Mean Feature Spectral"),
               subtitle = paste0("Levels ",unique(niv$Levels)),
               x = "Wavelength", y = "Reflectance (%)") 
    }
    if(j==2){
      p=p+labs(title = paste0("SD Feature Spectral"),
               subtitle = paste0("Levels ",unique(niv$Levels)),
               x = "Wavelength", y = "Standard Deviation (%)") 
    }
    if(j==3){
      p=p+labs(title = paste0("CV Feature Spectral"),
               subtitle = paste0("Levels ",unique(niv$Levels)),
               x = "Wavelength", y = "Coefficient Variation (%)") 
    }
    print(p)
    dev.off() 
    }
  }
}
