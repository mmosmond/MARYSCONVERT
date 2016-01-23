#here the dimensions of the ecosystem area are inputted. These numbers refer to the 1000m x 500m of the BCI forest
Lx<-1000
Ly<-500

#Here the area of the plots is inputted. This number gives the level of resolution of the analysis. You need to choose a resolution so to have enough statistics and at the same time reduces the noise. The output of the PCF gives you an idea of this trade off: you need to have a n>30 points and not too noisy
Asub<-100

###NOTE THAT Lx/Asub should be an integer, otherwise some plots will extend over the edge of the area or be lost. 

###################################################
################## Generate data ##################
set.seed(1234)

#number of samples
nsamp<-1000
#number of species
nspecies=3

#position of each sample will be drawn from a normal distribution with the following parameters (in the middle of the study area)
meanx<-Lx/2
meany<-Ly/2
sdx<-Lx/5
sdy<-Ly/5
xcoor<-rnorm(nsamp,mean=meanx,sd=sdx)
ycoor<-rnorm(nsamp,mean=meany,sd=sdy)
#Species names are given by sequential letters
specieslabels<-sample(letters[1:nspecies],nsamp,replace=TRUE)
#Abundances at each site are drawn from a uniform distribution between 10 and 100
abundances<-sample(10:100,nsamp, replace=TRUE)

################ Where data should be imported

#df is the data that will be imported, which should be a data frame with columns of xcoordinate, ycoordinate, species labels and then abundances.
df<-data.frame(xcoor,ycoor,specieslabels,abundances)

################ summary information
#species names:
USpecies<-levels(df[,3])
#number of different species
NSpecies<-length(USpecies)
#total number of individuals in the ecosystem
TotAbundance=sum(df[,4])

################ Plot information 
### plots are indexed 1,1 1,2 1,3 etc. and will then be converted back to real distances by multiplying by Asub

#number of plots that fit into the area, n.b., should be an integer
numx<-Lx/Asub
numy<-Ly/Asub

#gives the positions (from 1 to the total possible number)in X and Y of all possible plots (allpos)
allX<-rep(1:numx,each=numy)
allY<-rep(1:numy,numx)
allpos<-data.frame(allX,allY)
#the central position of each plot:
xyCpos<-allpos-0.5
#the total number of plots:
nplot<-nrow(allpos)

################ Which species are in each plot (possub)
#For each species point (row of the original data), what plot is it in?
positions<-ceiling(df[,1:2]/Asub)
#XYpos will be used to go through all the plots in allpos and return a list with the row number in the original data if present in that plot 
XYpos<-function(i){intersect(which(positions$xcoor==allX[i]),which(positions$ycoor==allY[i]))}
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
	table(df$specieslabels[poslist[[i]]])
	}
#CovFunc will use these counts to get the covariance between pairs 
CovFunc<-function(i,poslist1,poslist2){
	mean(SabA(i,poslist1)*SabA(i,poslist2))-mean(SabA(i,poslist1))*mean(SabA(i,poslist2))
}
allpair$allcov<-sapply(1:nrow(allpair),function(x)(CovFunc(x,allpair$possub1,allpair$possub2)))

###Take the mean covariance for each distance between plots
empPCF<-tapply(allpair$allcov,allpair$alldist,mean)
#convert row names back to real distances.
rownames(empPCF)<-as.numeric(rvec)*Asub

empPCFdf<-data.frame(as.numeric(rvec)*Asub,as.numeric(empPCF))
rownames(empPCFdf)<-1:nrow(empPCFdf)
colnames(empPCFdf)<-c("distance","PCF")

#model fitting
empPCFdf[empPCFdf$distance == 0,] <- NA #remove rows with distance of zero because this breaks the bessel function
m <- nls(PCF ~ 1 + 1/(2*pi) * (ro/lambda)^2 * besselK(abs(distance/lambda), 0), data = empPCFdf[-1,], start = list(ro = 10000, lambda = 50000), trace = T) #take absolute value because with made up data nls looks for negative lambda, which breaks bessel (set lower bounds somehow?)
#m <- nls(y ~ 1 + 1/(2*pi) * (ro/lambda) * besselK(x/lambda, 0), data = df, start = list(ro = 10000, lambda = 50000), trace = T, lower = list(0, 0))

summary(m) #parameter estimates
plot(empPCFdf) #data
lines(empPCFdf$distance ,predict(m)) #fitted curve