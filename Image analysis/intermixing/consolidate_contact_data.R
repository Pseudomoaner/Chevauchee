library(tidyverse)

root <- "D:/Sean/SurfaceColonyPIV/Fluorescence_Blocks"
setwd(root)
folders <- Sys.glob("*")
folders <- folders[!grepl("tif",folders)]
folders <- folders[!grepl("csv",folders)]
folders <- folders[!grepl("Flu_Piv_combo_analyze",folders)]
alldat <- list()

for (i in 1:length(folders)){
  setwd(paste0(root,.Platform$file.sep,folders[i]))
  print(paste0(root,.Platform$file.sep,folders[i]))
  dat <- read_csv("output.csv") %>% rename(contacts=pKs,flu_coverage=packFracs)
  dat$Experiment <- folders[i]
  dat$Time <- seq(0.5,dim(dat)[1]/2,0.5)
  alldat[[i]] <- dat
}

dat <- do.call("rbind",alldat)
dat <- dat %>% mutate(Density = case_when(grepl("OD1",Experiment) ~ "1",
                                          grepl("OD01",Experiment) ~ "0.1",
                                          grepl("OD001",Experiment) ~ "0.01",
                                          grepl("OD0001",Experiment) ~ "0.001"
                                          ),
                      Density = factor(Density,levels=c("1","0.1","0.01","0.001")),
                      Experiment = substr(Experiment,1,regexpr("_",Experiment)-1)
) 

ggplot(dat, aes(x=Time,y=contacts,color=Density))+geom_path(size=2,aes(group=Experiment))+theme_bw()

setwd(root)
write_csv(dat,"flu_dat.csv")