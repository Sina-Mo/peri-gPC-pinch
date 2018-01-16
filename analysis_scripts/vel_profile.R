# Clears any previous data
rm(list=ls())

# Sets working directory
setwd("/Users/Bosque/IBAMR/peri-gPC/newcode_runs/curvefiles")

#Setting up numbers for name convention
s1<- "newcode"
n<-681
num<-seq(1,n)
#num<-sprintf("%03d", num)

#Define parameters
sub<-12		# Number of initial time steps to ignore
d0=0.1		# Diameter of tube
Qall<-matrix(0,n,1)

for (k in 1:n){
	
	# Sets working directory
	wd<-paste("/Users/Bosque/IBAMR/peri-gPC/newcode_runs/curvefiles/",s1,num[k],"_profiles/",sep="")
	setwd(wd)
	# Displays the current working directory.
	getwd()
	# Gets file list.
	f<-list.files(path=".",all.files=FALSE)
	# Sets number of time steps based on number of files.
	ts<-length(f)
	
	# allocates space for flow rates
	Qts<-matrix(0,ts,1)
		
	for (i in 1:length(f)){
			
		data <- read.table(f[i]) 	# Loads data for time step
		s2<-dim(data)[1]			# Sets size of data file
		dr<-d0/s2					# Sets distance between data points
		q<-matrix(0,s2-1,1)			# allocates space for individual time step
		for (j in 1:(s2-1)){
				A<-0.5*pi*((d0-dr*(j-1))-(d0-dr*j)) # calculates concentric circle
				q[j]<-A*data[j,2]					# calculates volume flow for each circle
		}
		
		Qts[i]<-sum(q)				# sums individual volume flows for time step
		
	}

	Qall[k]<-mean(Qts[sub:ts])			# Calculates mean across simulation
	rm(Qts)
}
	plot(Qall)						# Plots data for each simulation
	setwd("/Users/Bosque/IBAMR/peri-gPC/newcode_runs/curvefiles")
	write.csv(Qall,file=paste(s1,"_Um_Q.csv",sep=""))	#Saves data for each simulation



