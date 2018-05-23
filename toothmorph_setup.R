########################################
# This document includes distance calculations for a 
# generic dataset as well as individual 
# "original", "liberal" and "conservative" datasets
# including use of dist.clust function
### Nameing protocol: 
# distpairs.xx   is data frame from distances.clust()
# distmat.xx     is distance matrix
# nmds#.xx       is #-dimensional ordination for the relevant distmat.xx 


##### Setup: Working Directory and Packages #####
# setwd("toothmorph")  #set your working directory wherever you want it to be
library(doParallel) #cluster
library(stratigraph) #for making range chart setup 

# save.image('RData/toothmorph_setup.RData')
# load('RData/toothmorph_setup.RData')

source('code/toothmorph/toothmorph_functions.R')  #


##### Functions required #####
#these are all in the toothmorph_functions.R file
# call.traits.csvs()
# toothdat.cleanup(toothdat, fix_dat = FALSE, sortby = 'age-obj')
# distances.clust(morph, traits, weights)
# distmat.fn(distpairs, type='avg')

##### Naming Protocol #####
# distpairs.xx   is data frame from distances.clust()
# distmat.xx     is distance matrix
# nmds#.xx       is #-dimensional ordination for the relevant distmat.xx 
# ad.xx          is the FAD and LAD for each tooth morphotype
# sc.xx          is the strat column object used as input for rangechart3 function

############### Morphometric Analysis Setup ###################
##### EXAMPLE: calculate pairwise distances between individually described teeth - generic code (uses morphotypes dataset as example) #####

## 1. Define traits, weights, and morphdat 

traits<-call.traits.csvs()

# assign weights for each trait (if no weights, distances.clust() will assign equal weight to each trait.)
weights<-c(1,1,1,1,1, 0.5, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 1, 1, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 1)
names(weights)<-c(names(traits), 'AR', 'LEN', 'WID')

toothdat<-read.csv('csv/Morphotypes.csv')
# toothdat<-read.csv('csv/Morph3.0_596_Compiled_orig.csv')
# toothdat<-read.csv('csv/Morph3.0_596_Compiled_con.csv')
# toothdat<-read.csv('csv/Morph3.0_596_Compiled_lib.csv')

morphdat<-toothdat.cleanup(toothdat, fix_dat = TRUE, sortby = 'morph')  #sortby = 'age-obj' is standard but doesn't work for morphotypes

## 2. Calculate pairwise Distances (and save output as .csv files)

distpairs<-distances.clust(morph=morphdat, traits=traits, weights=weights)
write.csv(distpairs, 'csv/pairwisedist_morph.csv')
distmat<-distmat.fn(distpairs)
write.csv(distmat, 'csv/distmat_morph.csv')

## 2 (alternate). Read in pairwise distances data from written .csv files

#distpairs
distpairs.read<-read.csv('csv/pairwisedist_morph.csv', header=T)
distpairs.read<-distpairs.read[,2:length(distpairs.read)] #remove extra column of row list

# distmat (needs some post-processing to coerce properly)
distmat.read<-read.csv('csv/distmat_morph.csv', header=T)
tooth.ID<-distmat.read[,1]
distmat.read<-distmat.read[,2:length(distmat.read)]
names(distmat.read)<-tooth.ID
rownames(distmat.read)<-tooth.ID
distmat.read<-as.matrix(distmat.read)
rm(tooth.ID)

# give them the right names and clean up
distmat <- distmat.read
distpairs <- distpairs.read
rm(distmat.read, distpairs.read) 


##### ACTUAL DATA: calculate pairwise distances between individually described teeth - Code/variable names for specific datasets #####
#### 0. pre-set information #####
traits<-call.traits.csvs()
weights<-c(1,1,1,1,1, 0.5, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 1, 1, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 1)
names(weights)<-c(names(traits), 'AR', 'LEN', 'WID')


#### 1. morphotypes dataset #####
toothdat.morph<-read.csv('csv/Morphotypes.csv')
morphdat.morph<-toothdat.cleanup(toothdat.morph, fix_dat = TRUE, sortby = 'morph')
rm(toothdat.morph) #clean up

## calculate distances and save output (!)
distpairs.morph<-distances.clust(morph=morphdat.morph, traits = traits, weights=weights)
write.csv(distpairs.morph, 'csv/pairwisedist_morph.csv')
distmat.morph<-distmat.fn(distpairs.morph, type='avg')
write.csv(distmat.morph, 'csv/distmat_morph.csv')

## alternatively, read in .csv files written previously...
# distpairs
distpairs.morph<-read.csv('csv/pairwisedist_morph.csv', header=T)
distpairs.morph<-distpairs.morph[,2:length(distpairs.morph)] #remove extra column of row list
# distmat 
distmat.morph<-read.csv('csv/distmat_morph.csv', header=T)
tooth.ID<-distmat.morph[,1]
distmat.morph<-distmat.morph[,2:length(distmat.morph)]
names(distmat.morph)<-tooth.ID
rownames(distmat.morph)<-tooth.ID
distmat.morph<-as.matrix(distmat.morph)
rm(tooth.ID) #clean up


#### 2. original dataset #####

## call the data
toothdat.orig<-read.csv('csv/Morph3.0_596_Compiled_orig.csv')
morphdat.orig<-toothdat.cleanup(toothdat.orig, fix_dat = TRUE, sortby = 'age-obj')
rm(toothdat.orig) #clean up

## calculate distances and save output (!)
distpairs.orig<-distances.clust(morph=morphdat.orig, traits = traits, weights=weights)
write.csv(distpairs.orig, 'csv/pairwisedist_orig.csv')
distmat.orig<-distmat.fn(distpairs.orig, type='avg')
write.csv(distmat.orig, 'csv/distmat_orig.csv')

## alternatively, read in .csv files you've written previously...
# distpairs
distpairs.orig<-read.csv('csv/pairwisedist_orig.csv', header=T)
distpairs.orig<-distpairs.orig[,2:length(distpairs.orig)] #remove extra column of row list
# distmat 
distmat.orig<-read.csv('csv/distmat_orig.csv', header=T)
tooth.ID<-distmat.orig[,1]
distmat.orig<-distmat.orig[,2:length(distmat.orig)]
names(distmat.orig)<-tooth.ID
rownames(distmat.orig)<-tooth.ID
distmat.orig<-as.matrix(distmat.orig)
rm(tooth.ID) #clean up


#### 3. Conservative Dataset #####
toothdat.con<-read.csv('csv/Morph3.0_596_Compiled_con.csv')
morphdat.con<-toothdat.cleanup(toothdat.con, fix_dat = TRUE, sortby = 'age-obj')
rm(toothdat.con) #clean up
## calculate distances and save output (!)
distpairs.con<-distances.clust(morph=morphdat.con, traits = traits, weights=weights)
write.csv(distpairs.con, 'csv/pairwisedist_con.csv')
distmat.con<-distmat.fn(distpairs.con, type='avg')
write.csv(distmat.con, 'csv/distmat_con.csv')

## alternatively, read in .csv files you've written previously...
# distpairs
distpairs.con<-read.csv('csv/pairwisedist_con.csv', header=T)
distpairs.con<-distpairs.con[,2:length(distpairs.con)] #remove extra column of row list
# distmat 
distmat.con<-read.csv('csv/distmat_con.csv', header=T)
tooth.ID<-distmat.con[,1]
distmat.con<-distmat.con[,2:length(distmat.con)]
names(distmat.con)<-tooth.ID
rownames(distmat.con)<-tooth.ID
distmat.con<-as.matrix(distmat.con)
rm(tooth.ID) #clean up


#### 4. liberal dataset  #####
toothdat.lib<-read.csv('csv/Morph3.0_596_Compiled_lib.csv')
morphdat.lib<-toothdat.cleanup(toothdat.lib, fix_dat = TRUE, sortby = 'age-obj')
rm(toothdat.lib) #clean up
## calculate distances and save output (!)
distpairs.lib<-distances.clust(morph=morphdat.lib, traits = traits, weights=weights)
write.csv(distpairs.lib, 'csv/pairwisedist_lib.csv')
distmat.lib<-distmat.fn(distpairs.lib, type='avg')
write.csv(distmat.lib, 'csv/distmat_lib.csv')

## alternatively, read in .csv files you've written previously...
# distpairs
distpairs.lib<-read.csv('csv/pairwisedist_lib.csv', header=T)
distpairs.lib<-distpairs.lib[,2:length(distpairs.lib)] #remove extra column of row list
# distmat 
distmat.lib<-read.csv('csv/distmat_lib.csv', header=T)
tooth.ID<-distmat.lib[,1]
distmat.lib<-distmat.lib[,2:length(distmat.lib)]
names(distmat.lib)<-tooth.ID
rownames(distmat.lib)<-tooth.ID
distmat.lib<-as.matrix(distmat.lib)
rm(tooth.ID) #clean up





#################### Range Chart Setup #################
##### EXAMPLE CODE: Create range chart objects for calculating evolutionary rates (uses orig dataset, as it has age data) #####

### Generic Code 
# Call in morphdat object
morphdat<-morphdat.orig

#build strat column object for use with package stratigraph
sc<-build.strat.obj(morphdat)

#create FAD/LAD object  
#####NOTE IN OLD CODE, INCREASING.DOWN = FALSE - NEED TO DOUBLE CHECK ALL THE THINGS AS YOU GO ALONG!!!
ad<-a.datums(sc, depths=sc$absolute.ages, increasing.down = TRUE)  #gives matrix of fads and lads for each type

### (optional step) Combine samples to increase number of objects in each age bin 
# inputs: (1) morphdat; (2) combines - list of vectors defining samples in each bin 
# Note that the combines vectors are currently defined in the function for clarity of this code.

# combine samples
morphdat<-combine.samples(morphdat)

#reset sc and ad objects
sc<-build.strat.obj(morphdat)
ad<-a.datums(sc, depths=sc$absolute.ages, increasing.down = TRUE)


##### CODE RUN: dataset specific code ##### 
# Original
toothdat.orig<-read.csv('csv/Morph3.0_596_Compiled_orig.csv')
morphdat.orig<-toothdat.cleanup(toothdat.orig, fix_dat = TRUE, sortby = 'age-obj')
rm(toothdat.orig) #clean up
morphdat.orig<-combine.samples(morphdat.orig)
sc.orig<-build.strat.obj(morphdat.orig)
ad.orig<-a.datums(sc.orig, depths=sc.orig$absolute.ages, increasing.down = TRUE)

# Conservative
toothdat.con<-read.csv('csv/Morph3.0_596_Compiled_con.csv')
morphdat.con<-toothdat.cleanup(toothdat.con, fix_dat = TRUE, sortby = 'age-obj')
rm(toothdat.con) #clean up
morphdat.con<-combine.samples(morphdat.con)
sc.con<-build.strat.obj(morphdat.con)
ad.con<-a.datums(sc.con, depths=sc.con$absolute.ages, increasing.down = TRUE)

# Liberal
toothdat.lib<-read.csv('csv/Morph3.0_596_Compiled_lib.csv')
morphdat.lib<-toothdat.cleanup(toothdat.lib, fix_dat = TRUE, sortby = 'age-obj')
rm(toothdat.lib) #clean up
morphdat.lib<-combine.samples(morphdat.lib) 
sc.lib<-build.strat.obj(morphdat.lib)
ad.lib<-a.datums(sc.lib, depths=sc.lib$absolute.ages, increasing.down = TRUE)


