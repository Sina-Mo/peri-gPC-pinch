# Clears any previous data
rm(list=ls())

# Sets working directory
setwd("/Users/Bosque/IBAMR/peri-gPC/newcode_runs/curvefiles")

# Displays the current working directory.
getwd()
s1<- "newcode"
n<-681
num<-seq(1,n)
num<-sprintf("%03d", num)

#allocate space
allofit_UmP<-matrix(0,n,11)
allofit_Ux<-matrix(0,n,8)
allofit_means<-matrix(0,n,8)

for (i in 1:n){

  # Constructs filenames
  filename_aorta_P<-paste(s1,num[i],"_aorta_P",sep="")
  filename_aorta_Um<-paste(s1,num[i],"_aorta_Um",sep="")
  filename_aorta_Ux<-paste(s1,num[i],"_aorta_Ux",sep="")
  
  filename_topcross_Um<-paste(s1,num[i],"_topcross_Um",sep="")
  filename_topcross_Ux<-paste(s1,num[i],"_topcross_Ux",sep="")
  
  filename_connect_Um<-paste(s1,num[i],"_connect_Um",sep="")
  filename_connect_Uy<-paste(s1,num[i],"_connect_Uy",sep="")
  
  filename_vena_P<-paste(s1,num[i],"_vena_P",sep="")
  filename_vena_Um<-paste(s1,num[i],"_vena_Um",sep="")
  filename_vena_Ux<-paste(s1,num[i],"_vena_Ux",sep="")
  
  #Loads data
  data_aorta_P_avg <- read.table(paste(filename_aorta_P,"_avg.curve",sep=""), header=FALSE, sep="")
  data_aorta_P_max <- read.table(paste(filename_aorta_P,"_max.curve",sep=""), header=FALSE, sep="")
  
  data_aorta_Um_avg <- read.table(paste(filename_aorta_Um,"_avg.curve",sep=""), header=FALSE, sep="")
  data_aorta_Um_max <- read.table(paste(filename_aorta_Um,"_max.curve",sep=""), header=FALSE, sep="")
  
  data_aorta_Ux_avg <- read.table(paste(filename_aorta_Ux,"_avg.curve",sep=""), header=FALSE, sep="")
  data_aorta_Ux_max <- read.table(paste(filename_aorta_Ux,"_max.curve",sep=""), header=FALSE, sep="")
  
  data_topcross_Um_avg <- read.table(paste(filename_topcross_Um,"_avg.curve",sep=""), header=FALSE, sep="")
  data_topcross_Um_max <- read.table(paste(filename_topcross_Um,"_max.curve",sep=""), header=FALSE, sep="")
  
  data_topcross_Ux_avg <- read.table(paste(filename_topcross_Ux,"_avg.curve",sep=""), header=FALSE, sep="")
  data_topcross_Ux_max <- read.table(paste(filename_topcross_Ux,"_max.curve",sep=""), header=FALSE, sep="")
  
  data_connect_Um_avg <- read.table(paste(filename_connect_Um,"_avg.curve",sep=""), header=FALSE, sep="")
  data_connect_Um_max <- read.table(paste(filename_connect_Um,"_max.curve",sep=""), header=FALSE, sep="")
  
  data_connect_Uy_avg <- read.table(paste(filename_connect_Uy,"_avg.curve",sep=""), header=FALSE, sep="")
  data_connect_Uy_max <- read.table(paste(filename_connect_Uy,"_max.curve",sep=""), header=FALSE, sep="")
  
  data_vena_P_avg <- read.table(paste(filename_vena_P,"_avg.curve",sep=""), header=FALSE, sep="")
  data_vena_P_max <- read.table(paste(filename_vena_P,"_max.curve",sep=""), header=FALSE, sep="")
  
  data_vena_Um_avg <- read.table(paste(filename_vena_Um,"_avg.curve",sep=""), header=FALSE, sep="")
  data_vena_Um_max <- read.table(paste(filename_vena_Um,"_max.curve",sep=""), header=FALSE, sep="")
  
  data_vena_Ux_avg <- read.table(paste(filename_vena_Ux,"_avg.curve",sep=""), header=FALSE, sep="")
  data_vena_Ux_max <- read.table(paste(filename_vena_Ux,"_max.curve",sep=""), header=FALSE, sep="")
  
  #Plot things
  #pdf(paste(filename_topcross_Um,".pdf",sep=""))
  #plot(V2~V1,data=data_topcross_Um_avg,type="l",col="black",
  #     xlab=list("Time (s)",cex=1),ylab=list("data",cex=1),ylim=c(min(data_topcross_Um_avg$V2),max(data_topcross_Um_avg$V2)))
  #dev.off()
  
  #Calculating values
  allofit_UmP[i,1]<-max(data_topcross_Um_avg$V2)
  allofit_UmP[i,2]<-max(data_topcross_Um_max$V2)
  allofit_UmP[i,3]<-max(data_connect_Um_avg$V2)
  allofit_UmP[i,4]<-max(data_connect_Um_max$V2)
  allofit_UmP[i,5]<-max(data_aorta_P_max$V2)
  allofit_UmP[i,6]<-max(data_aorta_Um_avg$V2)
  allofit_UmP[i,7]<-max(data_aorta_Um_max$V2)
  allofit_UmP[i,8]<-max(data_vena_P_max$V2)
  allofit_UmP[i,9]<-max(data_vena_Um_avg$V2)
  allofit_UmP[i,10]<-max(data_vena_Um_max$V2)

  Pdiff_max<-data_aorta_P_max$V2 - data_vena_P_max$V2
  Pdiff_max[ is.infinite(Pdiff_max) ]<-0
  allofit_UmP[i,11]<-max(Pdiff_max)
  

  #Calculating values
  allofit_Ux[i,1]<-max(data_topcross_Ux_avg$V2)
  allofit_Ux[i,2]<-max(data_topcross_Ux_max$V2)
  allofit_Ux[i,3]<-max(data_connect_Uy_avg$V2)
  allofit_Ux[i,4]<-max(data_connect_Uy_max$V2)
  allofit_Ux[i,5]<-max(data_aorta_Ux_avg$V2)
  allofit_Ux[i,6]<-max(data_aorta_Ux_max$V2)
  allofit_Ux[i,7]<-max(data_vena_Ux_avg$V2)
  allofit_Ux[i,8]<-max(data_vena_Ux_max$V2)

  
  #reshape data for mean calculation, takes off the first 0.5 second of data, averages the remainder. 
  allofit_means[i,1]<-mean(data_topcross_Um_avg[12:length(data_topcross_Um_avg$V1),2])
  allofit_means[i,2]<-mean(data_topcross_Ux_avg[12:length(data_topcross_Ux_avg$V1),2])
  allofit_means[i,3]<-mean(data_connect_Um_avg[12:length(data_connect_Um_avg$V1),2])
  allofit_means[i,4]<-mean(data_connect_Uy_avg[12:length(data_connect_Uy_avg$V1),2])
  allofit_means[i,5]<-mean(data_aorta_Um_avg[12:length(data_aorta_Um_avg$V1),2])
  allofit_means[i,6]<-mean(data_aorta_Ux_avg[12:length(data_aorta_Ux_avg$V1),2])
  allofit_means[i,7]<-mean(data_vena_Um_avg[12:length(data_vena_Um_avg$V1),2])
  allofit_means[i,8]<-mean(data_vena_Ux_avg[12:length(data_vena_Ux_avg$V1),2])
}

colnames(allofit_UmP)<-c("upper avg","upper max","connect avg","connect max","aorta p","aorta avg","aorta max","vena p","vena avg","vena max","delta p max")

colnames(allofit_Ux)<-c("upper Ux avg","upper Ux max","connect Ux avg","connect Ux max","aorta Ux avg","aorta Ux max","vena Ux avg","vena Ux max")

colnames(allofit_means)<-c("upper Um avg","upper Ux avg","connect Um avg","connect Ux avg","aorta Um avg","aorta Ux avg","vena Um avg","vena Ux avg")

#Saves it all into a CSV file
write.csv(allofit_UmP,file=paste(s1,"_alldata_UmP.csv",sep=""))
write.csv(allofit_Ux,file=paste(s1,"_alldata_Ux.csv",sep=""))
write.csv(allofit_means,file=paste(s1,"_alldata_means.csv",sep=""))
