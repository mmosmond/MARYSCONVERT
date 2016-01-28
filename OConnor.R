############### Import Patrick's data:

load("/Users/mscott0106/Library/Mail Downloads/Simulation abundances.RData")
df1<-Abundance.df
### Patrick's data is given as abundances in each patch of a 20x20 grid. The R code below (based on the Mathematica package from Azaele, 2015) does not use this type of abundance data - the number of rows for each species in a "plot" (grid cell) is used as the number of counts in that plot. We can replicate each row according to the abundance column to make the data fit the format required for the package. 
library(splitstackshape)
df<-expandRows(df1,"Abundance")
# some rows are removed because the species is absent. These species are removed as factors. 
df$Species<-droplevels(df$Species)
###Note, Patrick's data gives the y coordinate in the first column and x coordinate in the second column. The grid is 20*20 so this doesn't matter in this case. 


############### Import Robin's data:
df1<-read.csv("~/Downloads/sbc_fish_2010.csv")
head(df1)
max(df1$Latitude)
min(df1$Latitude)
max(df1$Longitude)
min(df1$Longitude)

library(splitstackshape)
df2<-expandRows(df1,"Abundance")

df3<-df2[,6:7]
df3$SpeciesName<-as.factor(with(df2,paste(Genus,Species)))
nrow(df3)
df3$Latitude<-(df3$Latitude-min(df1$Latitude))*10+0.01
df3$Longitude<-(df3$Longitude-min(df1$Longitude))*10+0.01
df<-df3

plot(df[,1:2])


############### Plot Dimensions
#here the dimensions of the ecosystem area are inputted.
Lx<-4.5
Ly<-8
#Here the area of the plots is inputted. This number gives the level of resolution of the analysis. You need to choose a resolution so to have enough statistics and at the same time reduces the noise. The output of the PCF gives you an idea of this trade off: you need to have a n>30 points and not too noisy
Asub<-0.5
###NOTE THAT Lx/Asub should be an integer, otherwise some plots will extend over the edge of the area or be lost. 

################ summary information
#species names:
USpecies<-levels(df[,3])
#number of different species
NSpecies<-length(USpecies)

################ Plot information 
### plots are indexed (1,1) (1,2) (1,3) etc. and will then be converted back to real distances by multiplying by Asub

#number of plots that fit into the area, n.b., should be an integer
numx<-Lx/Asub
numy<-Ly/Asub

#gives the positions (from 1 to the total possible number)in X and Y of all possible plots (allpos)
allX<-rep(1:numx,each=numy)
allY<-rep(1:numy,numx)
allpos<-data.frame(allX,allY)
#the central position of each plot (used for calculating distances)
xyCpos<-allpos-0.5
#the total number of plots:
nplot<-nrow(allpos)

################ Which species are in each plot (possub)
#For each species point (row of the original data), what plot is it in?
positions<-ceiling(df[,1:2]/Asub)
#XYpos will be used to go through all the plots in allpos and return a list with the row number in the original data if present in that plot 
XYpos<-function(i){intersect(which(positions[,1]==allX[i]),which(positions[,2]==allY[i]))}
allpos$possub<-sapply(1:nplot,function(x)(XYpos(x)))

################ Pairwise comparisons between plots (will also copy across possub)
#all possible combinations of plot positions
allpair<-data.frame(rep(allX,each=nplot),rep(allY,each=nplot),rep(allX,nplot),rep(allY,nplot))
colnames(allpair)<-c("Xpair1","Ypair1","Xpair2","Ypair2")
#I will also transport across the possub so that the original data can be returned for the first and second parts of the pair
allpair$possub1<-rep(allpos$possub,each=nplot)
allpair$possub2<-rep(rep(allpos$possub,nplot))


################ Calculate distances between each plot (alldist)
# Rather than using the indexes to each plot, the central position of each plot will be used to calculate the distances between particular pairs of plots (note that, this is not necessary, distances are just measured relative to one another so you could use the indexes)
allpairCpos<-cbind(rep(allX,each=nplot),rep(allY,each=nplot),rep(allX,nplot),rep(allY,nplot))-0.5
#dist takes the euclidian distance using pythagoras' rule
dist<-function(i){
sqrt((allpairCpos[i,3]-allpairCpos[i,1])^2+(allpairCpos[i,4]-allpairCpos[i,2])^2)
}
npair<-nrow(allpair)
allpair$alldist<-as.factor(sapply(1:npair,function(x)(dist(x))))
# rvec gives the catagories for the distances between plots. Below the average covariance for each distance is calculated.
rvec<-levels(allpair$alldist)

################ Calculate covariances between each pair (allcov)
#SabA gets the counts of the species in each plot using poslist
SabA<-function(i,poslist){
	table(df[poslist[[i]],3])
	}
#CovFunc will use these counts to get the covariance between pairs 
CovFunc<-function(i,poslist1,poslist2){
	mean(SabA(i,poslist1)*SabA(i,poslist2))/(mean(SabA(i,poslist1))*mean(SabA(i,poslist2)))
}
allpair$allcov<-sapply(1:nrow(allpair),function(x)(CovFunc(x,allpair$possub1,allpair$possub2)))
sapply(1:50,function(x)(CovFunc(x,allpair$possub1,allpair$possub2)))

### because there are some grid cells with no data, I simply omit these rows from allpair. 
allpair<-na.omit(allpair)

###Take the mean covariance for each distance between plots
allpair$alldist<-droplevels(allpair$alldist)
rvec<-levels(allpair$alldist)
empPCF<-tapply(allpair$allcov,allpair$alldist,mean)
#put the PCF values and their distances in a data frame with appropriate columns
empPCFdf<-data.frame(as.numeric(rvec)*Asub,as.numeric(empPCF))
rownames(empPCFdf)<-1:nrow(empPCFdf)
colnames(empPCFdf)<-c("distance","PCF")

plot(empPCFdf)

#We will fit a model that is always Infinite when distance=0. So we remove the first row to help fit
x<-empPCFdf$distance
y<-empPCFdf$PCF
dataset<-data.frame(x,y)
dataset<-dataset[-1,]

#model often doesn't converge well enough. 
modelfit<-nls(y ~ 1+(1/(2*pi)) * ((ro/lambda)^2) * besselK(x/lambda,0), data=dataset, start=list(ro=50,lambda=1),control=list(minFactor = 1/1024,warnOnly=FALSE,maxiter=200))

foo<-summary(modelfit)
params<-c(foo$parameters[1,1],foo$parameters[2,1])
roEst<-params[1]
lamEst<-params[2]

gFunc<-function(x,ro,lambda){
	 1+(1/(2*pi)) * ((ro/lambda)^2) * besselK(x/lambda,0)
}
plot(empPCFdf)
lines(seq(0,20,by=0.01),sapply(seq(0,20,by=0.01),function(x)(gFunc(x, roEst, lamEst))))

empPCFdf

#number of individuals in each plot
ind<-sapply(1:nplot,function(x)(sum(SabA(x,allpos$possub)))) 
#total number of individuals
totalind<-sum(ind)

###functions taken from Mathematica file, used to generate SAR and SAD curved
alpha2p<-function(L,ro,lambda){
	pi * ((L/ro)^2)*(1-(2*lambda/L) * (besselK(L/lambda,1) * besselI(L/lambda,1)) / (besselK(L/lambda,1) * besselI(L/lambda,0) + besselK(L/lambda,0) * besselI(L/lambda,1) ))^(-1)
}
beta2<-function(L,ro,lambda){
	(totalind/(NSpecies * Lx * Ly))* ro^2 *(1-(2*lambda/L) * (besselK(L/lambda,1) * besselI(L/lambda,1)) / (besselK(L/lambda,1) * besselI(L/lambda,0) + besselK(L/lambda,0) * besselI(L/lambda,1) ))
}
rsapres2medp<-function(n,L,ro,lambda){
	pgamma(2^(n+1),alpha2p(L,ro,lambda),scale=beta2(L,ro,lambda))-pgamma(2^n,alpha2p(L,ro,lambda),scale=beta2(L,ro,lambda))
}
sarmedp<-function(L,ro,lambda){
	1-pgamma(1,alpha2p(L,ro,lambda),scale=beta2(L,ro,lambda))
}

#SAR used for downscaled SAR at radius r when all the information in the study region is available. The largest scale is Lx*Ly. Here, r<sqrt(Lx*Ly/pi)
SAR<-function(r,ro,lambda){
	NSpecies*(sarmedp(r,ro,lambda))/(sarmedp(sqrt((Lx*Ly)/pi),ro,lambda))
}

plot(SAR(1:14,roEst,lamEst),xlab="Radius (0.1 lat-long coordinates, not corrected) of sampled area",ylab="Mean number of species")

upSAR<-function(r,ro,lambda){
	NSpecies*(sarmedp(r,ro,lambda))/(sarmedp(sqrt((20*0.1)/pi),ro,lambda))
}

plot(upSAR(1:14,roEst,lamEst),xlab="Radius (0.1 lat-long coordinates, not corrected) of sampled area",ylab="Mean number of species")


#Downscaled SAD at radius r when all the information in the study region is available. The largest scale is Lx*Ly. Here r<sqrt((Lx*Ly)/pi)
SAD<-function(n,r,ro,lambda){
	NSpecies*(rsapres2medp(n,r,ro,lambda))/(sarmedp(sqrt((Lx*Ly)/pi),ro,lambda))
}

plot(SAD(1:20, 5,roEst,lamEst),xlab="Preston Classes",ylab="Number of species")
