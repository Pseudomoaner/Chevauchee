library(tidyverse)
library(imager)

#this script takes all the PIV measurement outputs from each experiment and computes median velocities for each timepoint

setwd("D:/Sean/SurfaceColonyPIV")

files <- Sys.glob("*new_v-mags.csv")

setwd("D:/Sean/SurfaceColonyPIV/Fluorescence_Blocks")
f_dirs<-list.dirs(recursive=F)[c(1:12)]

setwd("D:/Sean/SurfaceColonyPIV")

alldat <- list()

for (i in 1:length(files)) {
    print(paste0("i is: ",i))
    print(paste0("Loading file: ",files[i]))
    dat <- data.table::fread(files[i]) %>% as_tibble()  #use fread for speed with large files
    colnames(dat) <- c("Time","x","y","mag_vector", "mean-vector_mag")
    dat <- dat %>% mutate(`mag_vector`= `mag_vector`*(0.227/2), `mean-vector_mag` = `mean-vector_mag`*(0.227/2)) 
    
    dat$Density <- substr(files[i],gregexpr("_",files[i])[[1]][2]+1, gregexpr("_",files[i])[[1]][3]-1)
    dat$Experiment <- substr(files[i],1, gregexpr("_",files[i],length(files[i]))[[1]][1]-1)

    mask <- load.image(paste0(sub("\\.","D:/Sean/SurfaceColonyPIV/Fluorescence_Blocks",f_dirs[which(substr(files[i],1,regexpr("_",files[i])-1)==substr(f_dirs,3,regexpr("_",f_dirs)-1))]),"/","C1_segment.tif"))
    if(grepl("Feb22",files[i])){
      mask <-imsub(mask,x%inr%c(4,2751),y%inr%c(4,2211)) #account for Feb22 weird sized mask
    }
    #make the PIV data and fluorescence mask the same size if they aren't
    sc_mask <- mask %>% imsub(x%inr%c((dim(mask)[1]-max(dat$x)*4)/2+1,dim(mask)[1]-(dim(mask)[1]-max(dat$x)*4)/2),
          y%inr%c((dim(mask)[2]-max(dat$y)*4)/2+1,dim(mask)[2]-(dim(mask)[2]-max(dat$y)*4)/2)) %>% resize(size_x=max(dat$x),size_y=max(dat$y))

    maskdf <- as.data.frame(sc_mask) %>% mutate(Time=z/2,mask=ifelse(value>0,1,NA)) %>% select(-value,-z) %>% as_tibble()
    masked <- left_join(dat,maskdf)
    
    out <-masked %>% as_tibble() %>% mutate(mag = `mean-vector_mag`*mask) %>% select(-x,-y,-`mean-vector_mag`,-mag_vector,-mask) %>%
      group_by(Time,Density,Experiment) %>% 
      mutate(coverage = sum(!is.na(mag))/n(),
             covered = (coverage>0.99)) %>% ungroup() %>%
      mutate(mag = ifelse(is.na(mag),ifelse(Time>Time[min(which(covered))],0,NA),mag)) 
    
    alldat[[i]] <-out %>% group_by(Time,Density,Experiment)  %>%summarise(coverage = sum(!is.na(mag))/n(),
                                   mean=mean(mag,na.rm=T),
                                   sd=sd(mag,na.rm=T),
                                   median=median(mag,na.rm=T),
                                   log_mean=mean(log10(ifelse(mag==0,10^-5,mag)),na.rm=T),
                                   log_sd=sd(log10(ifelse(mag==0,10^-5,mag)),na.rm=T),
                                   log_median=median(log10(ifelse(mag==0,10^-5,mag)),na.rm=T)
    ) %>% 
      ungroup()
    print(ggplot(out,aes(x=log10(ifelse(mag==0,10^-5,mag)),color=as.factor(Time)))+geom_density(size=1.5)+scale_color_viridis_d(option="C")+facet_wrap(.~Time)+labs(title=files[i])+theme_bw())
    
    print(ggplot(out,aes(x=log10(ifelse(mag==0,10^-5,mag)),color=as.factor(Time)))+geom_density(size=1.5)+scale_color_viridis_d(option="C")+labs(title=files[i])+theme_bw())
    
}


sdat <- do.call("rbind",alldat)
sdat <- sdat %>% mutate(Density = case_when(Density=="OD0001" ~ 0.001,Density=="OD001" ~ 0.01,Density=="OD01" ~ 0.1,Density=="OD1" ~ 1))
sdat$Density <- factor(sdat$Density,levels=c(1,0.1,0.01,0.001))
sdat <- sdat %>% group_by(Experiment,Density) %>% 
  mutate(max_time = Time[which.max(mean)],
         covered=(coverage>0.99),
         covered_time=Time[(min(which(covered)))]) %>% ungroup() %>%
  group_by(Density) %>% mutate(group_min_time = min(Time-covered_time),group_min_v_time = min(Time-max_time)) %>% ungroup()

write_csv(sdat,"PIV_dat_V7_f-masked.csv")