#########CONSTRUCCION DE MATRIX DE FIRMAS###########
# 
# dir_signature = "D:/OneDrive-CGIAR/Workspace/analys_spectro/SPECTROMETER/analisis_espectro";Treatment = c("T1","T2","T3","T4");Levels = c("N1","N3","N5","N7","N9");
# Date_start = "2016-12-28";extension = ".txt";decimal=".";Separations=" ";delte_line=1;Wavelength_start=400;Wavelength_end=950;Col_Wavelength=3;Col_Reflectance=4  
# 
# dir_signature = "D:/OneDrive-CGIAR/Workspace/app_analisis_spectral/demo2";Treatment = c("K","N");Levels = c("100","50","0");Replicates = c(9,12,3)
# Date_start = "2017-05-04";extension = ".TRM";decimal=".";Separations=" ";delte_line=2;Wavelength_start=400;Wavelength_end=950;Col_Wavelength=2;Col_Reflectance=4

dir_signature="D:/OneDrive-CGIAR/Workspace/analys_spectro/SPECTROMETER/beans/spectro";
Treatment=c("NI","IR");Levels=c("AB","AD");Replicates=c("M","A");
Date_start="2018-09-05";extension=".txt";decimal=".";Separations=",";delte_line=2;Wavelength_start=400;Wavelength_end=950;Col_Wavelength=2;Col_Reflectance=3



firmasSpectrales=lee_firmas(dir_signature = "D:/OneDrive-CGIAR/Workspace/app_analisis_spectral/demo2",Treatment = c("K","N"),Levels = c("100","50","0"),Replicates = c(9,12,3),
                            Date_start = "2017-05-04",extension = ".TRM",decimal=".",Separations=" ",delte_line=2,Wavelength_start=400,Wavelength_end=950,Col_Wavelength=2,Col_Reflectance=4)
########estadistica firmas#############
firm_stats=stats_fir(firmasSpectrales,Replicate = T)
##################CAL VEGETATION INDEX###############
IV=cal_IV(firmasSpectrales,Replicate = T)
############plot stats
plots_stats(firm_stats)
