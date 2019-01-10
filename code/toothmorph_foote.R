########################################################################
#
#     Toothmorph Evolutionary Rates: Foote
#
#
#
######################################################################

## Notes for Foote, MARK analyses 
#  in/out:  morphdat - full dataset of all teeth to include in analysis
#           sc.xxx - strat column object - xxx corresponds to appropriate dataset
#           ranges - counts of each morphotype at each age; can be used for make.inp; also range.sc$counts
#           ad - table of FAD and LAD
#           taxa - list of taxa
#           AgeID.unique - list of dates used in plotting, etc. 

## things that might be able to be deleted: 
# (morphdat.con, morphdat.lib, morphdat.morph, morphdat.orig)

##### Load/save datasets #####

# load('RData/toothmorph_rates.RData')
# save.image('RData/toothmorph_rates.RData')
source('code/toothmorph/toothmorph_functions.R')  #



# useful variables that should be around already but may not be for whatever reason
# need to specify dataset for these to work.
AgeID<-morphdat$AgeID
AgeID.unique<-unique(AgeID)
taxa<-levels(morphdat$Morphotype_name)[2:length(levels(morphdat$Morphotype_name))] #excludes the '' blank; should have 136 items

##### 1. Foote Boundary Crossers - more notes in foote.R #####
# output:   pp - per capita origination
#           qq - per capita extinction
#           dt - time interval (youngest to oldest)
#           dt.rev - time intervals (oldest to youngest)
###############
# a) count up crossers, Using AgeID.unique and the dataset 'ad' 
# to calculate foote rates
# code from above section for 'ad' object:
# ad<-a.datums(range.sc, depths=range.sc$absolute.ages)  #gives matrix of fads and lads for each type

dt<-c()                           # change in time (interval)
for(i in 1:(length(AgeID.unique))-1) {
   dt[i]<-(AgeID.unique[i+1]-AgeID.unique[i]) }
dt.rev<-dt[length(dt):1]          #in reverse


## select a dataset
ad <- ad.con

## calculate singletons, crossers present in each time bin - yay for for-loops?
Nfl<-c()  # osingletons (no boundary crossers)
Nbl<-c()  # bottom crossers only
Nft<-c()  # top crossers only
Nbt<-c()  # both bottom and top crossers
for (i in 1:length(AgeID.unique)) {
   #subset of all species present or assumed present at the time point [i]
   ad.sub<-subset(ad, fads >= AgeID.unique[i] & lads <= AgeID.unique[i])
   #how many of each type of occurrance are present? 
   Nfl[i] <- length(subset(ad.sub, fads==lads)[,1])  # "singletons" (no boundary cross)
   Nbl[i] <- length(subset(ad.sub, fads>=AgeID.unique[i] & lads==AgeID.unique[i] & fads!=lads)[,1]) # bottom crossers only
   Nft[i] <- length(subset(ad.sub, fads==AgeID.unique[i] & lads<=AgeID.unique[i] & fads!=lads)[,1]) # top crossers only 
   Nbt[i] <- length(subset(ad.sub, fads> AgeID.unique[i] & lads< AgeID.unique[i])[,1]) # both top and bottom  crossers
}

# b) calculate relevant metrics
Ntot <- apply(cbind(Nfl, Nbl, Nft, Nbt), 1, sum) # total Diversity observed
Nb   <- apply(cbind(Nbl, Nbt), 1, sum)           # All bottom boundary crossers
Nt   <- apply(cbind(Nft, Nbt), 1, sum)           # all top boundary crossers
No   <- apply(cbind(Nfl, Nft), 1, sum)           # number of originations
Ne   <- apply(cbind(Nfl, Nbl), 1, sum)           # number of extinctions
Ndiv <- apply(cbind(Nb, Nt), 1, function(x) sum(x)/2)     # estimated mean standing diversity
pp   <- apply(cbind(Nbt, Nt, c(dt,0)), 1, function(x) {-log(x[1]/x[2]) / x[3]}) # per capita origination
qq   <- apply(cbind(Nbt, Nb, c(0,dt)), 1, function(x) {-log(x[1]/x[2]) / x[3]}) # per capita extinction

# c) cleanup
rm(Nfl, Nft, Nbl, Nbt, ad.sub, Nb, Nt,i)

# d) summary plot
plot(AgeID.unique, pp, col='blue', pch=16, type='o', xlim=c(73, 42),
   main='Foote 2000 boundary crosser extinction and origination')  #AgeID.unique
points(AgeID.unique, qq, col='red', pch=16, type='o')
legend ('topright', legend=c('origination', 'extinction'), col=c('blue', 'red'), lty=1, pch=16)

# add diversity estimated
par(new=T)
plot(AgeID.unique, Ndiv, axes=F, xlim=c(73, 42), ylim=c(0, 85))
axis(4)

# add diversity total
points(AgeID.unique, Ntot, pch=16)


