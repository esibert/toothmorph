#####################################
#  Code for calculating NMDS of tooth assemblages
# Dependencies:   toothmorph_dist_setup.RData
#                 toothmorph_functions.R
#                 toothmorph_rangechart_setup.RData
#########################################

# setwd("")
source('code/toothmorph/toothmorph_functions.R')
library(doParallel) #cluster
library(vegan)
# save.image('RData/toothmorph_nmds.RData')
# load('RData/toothmorph_nmds.RData')

##### ABOUT THIS DOCUMENT #####
# Dependent datasets (combined to above toothmorph_nmds.RData) call in both distances and range chart datasets: 
# load('RData/toothmorph_dist_setup.RData') #includes distpairs, distmat, and nmds.list 
# load('RData/toothmorph_rangechart_setup.RData')  #ad, sc, etc. for plotting

# Nameing protocol: 
# distpairs.xx   is data frame from distances.clust()
# distmat.xx     is distance matrix
# nmds#.xx       is #-dimensional ordination for the relevant distmat.xx 

# Time bins: 
# Cretaceous: >66
# Pulse 1: 66-60
# Pulse 2: 60-55
# Eocene: <55
time.splits<-c(66,60,55,43)  #K/Pg, P1/P2, P2/E, Eoc pre-end time point


##### Step 0: NMDS Run and setup #####
load('RData/toothmorph_dist_setup.RData') 
#generate NMDS ordination objects: 

nmds.nam<-c(paste('NMDS', 2:5, '.orig', sep=''), paste('NMDS', 2:5, '.con', sep=''), paste('NMDS', 2:5, '.lib', sep=''), paste('NMDS', 2:5, '.morph', sep=''))

# make cluster 
cores<-detectCores()
cl<-makeCluster(cores-2) #detect the number of cores and leave one open
registerDoParallel(cl)

# run nmds's (note this is all done in parallel) 
nmds.list<-foreach(i=1:length(nmds.nam), .export=c('distmat.orig', 'distmat.con', 'distmat.lib', 'distmat.morph'), .packages=c('vegan')) %dopar% {
   nam<-paste(nmds.nam[i])
   dims<-as.numeric(substr(nam, 5,5))   #gives number of NMDS dimensions
   dataset<-unlist(strsplit(nam, '.', fixed=T))[2]  #gives "orig", "con", "lib"
   dataset<-get(paste('distmat.', dataset, sep=''))  #call the distance matrix
   nmds<-metaMDS(dataset, k=dims)
   nmds
}

# kill cluster
stopCluster(cl)
rm(cl, cores)

# gve nmds.list the appropriate names
names(nmds.list)<-nmds.nam


##### (a) Make NMDS Plots (stressplot, ordiplot, timeseries, grouped) #####

nmds.plots.all<-function(nmds.list, time.splits, output.type = "pdf") {
   for(i in 1:length(nmds.list)) { 
      # call appropriate NMDS results and tooth ageID's for plotting
      nmds<-nmds.list[[i]]   
      # dataset is a character vector defining which of the 3 datasets (orig/con/lib) are being manipulated
      dataset<-unlist(strsplit(nmds.nam[i], '.', fixed=T))[2]  #gives "orig", "con", "lib"
      # ageID vector is the ages of each individual tooth in the dataset
      ageID<-get(paste('morphdat.', dataset, sep=''))$AgeID
      
      # Call NMDS points
      nmds.df<-data.frame(cbind(ageID, nmds$points)) #Make data frame for plot
      dims<-grep("MDS", names(nmds.df))         # call all MDS columns
      combs<-t(combn(dims,2))                   # combinations for each
      
      # set wd to put these plots (orig, con, lib)
      # setwd(paste('figures/nmds_all/nmds_', dataset, sep='')) 
   
      # stressplot figure
      if(output.type == "pdf") {pdf(file=paste('figures/nmds_all/', nmds.nam[i], '_stressplot.pdf', sep=''), height=7, width=7, useDingbats = FALSE)}
      if(output.type == "jpg") {jpeg(filename=paste('figures/nmds_all/', nmds.nam[i], '_stressplot.jpg', sep='')) }
      stressplot(nmds, main=paste(nmds.nam[i], '_stressplot', sep=''))
      dev.off()
      
      for(j in 1:length(combs[,1])) {
         
         # call relevant MDS values
         c<-combs[j,]
         
         # Make text-ordination plot
         cc<-c-1
         if(output.type == "pdf") {pdf(file=paste('figures/nmds_all/', nmds.nam[i], '_', cc[1], '-', cc[2],'_ordiplot.pdf', sep=''), height=10, width=10, useDingbats=FALSE)}
         if(output.type == "jpg") {jpeg(filename=paste('figures/nmds_all/', nmds.nam[i], '_', cc[1], '-', cc[2],'_ordiplot.jpg', sep=''), height=1600, width=1600) }
         plot(nmds, choices=cc, type='t', main=paste(nmds.nam[i], cc[1], cc[2],'ordiplot', sep='_'))
         dev.off()
         
         #Make timeseries figure
         if(output.type == "pdf") {pdf(file=paste('figures/nmds_all/', nmds.nam[i], '_', cc[1], '-', cc[2], '_timeseries.pdf', sep='')) }
         if(output.type == "jpg") {jpeg(filename=paste('figures/nmds_all/', nmds.nam[i], '_', cc[1], '-', cc[2], '_timeseries.jpg', sep='')) } ########Change filename here!
         plot(nmds.df[,c[1]], nmds.df[,c[2]], type='n', xlab=paste('MDS', cc[1], sep=''), ylab=paste('MDS', cc[2], sep=''), main=paste(nmds.nam[i], cc[1], cc[2], 'timeseries', sep='_'))
         cols<-rainbow(length(AgeID.unique))
         # cols<-colfunc(length(AgeID.unique))
         for(k in 1:length(AgeID.unique)) {
            sub<-subset(nmds.df, ageID==AgeID.unique[k])
            points(sub[,c[1]], sub[,c[2]], pch=k, col=cols[k])
            Plot_ConvexHull(sub[,c[1]], sub[,c[2]], lcol=cols[k], lwd=2)
         }
         dev.off() 
         
         #Make clumped figure
         pal<-subset(nmds.df, ageID>=time.splits[3] & ageID<=time.splits[1])
         eoc<-subset(nmds.df, ageID<=time.splits[3])
         cret<-subset(nmds.df, ageID>=time.splits[1])
         ages.grouped<-list(eoc, pal, cret)
         
         if(output.type == "pdf") {pdf(file=paste('figures/nmds_all/', nmds.nam[i], '_', cc[1], '-', cc[2], '_grouped.pdf', sep=''), height=10, width=10, useDingbats = FALSE)}
         if(output.type == "jpg") {jpeg(filename=paste('figures/nmds_all/', nmds.nam[i], '_', cc[1], '-', cc[2], '_grouped.jpg', sep=''))  }
         plot(nmds.df[,c[1]], nmds.df[,c[2]], type='n', xlab=paste('MDS', cc[1], sep=''), ylab=paste('MDS', cc[2], sep=''), main=paste(nmds.nam[i], cc[1], cc[2], 'grouped', sep='_'))
         cols<-rainbow(length(ages.grouped))
         # Step 2: loop through the dataset to add to plot
         for(k in 1:length(ages.grouped)) {
            sub<-data.frame(ages.grouped[k])
            points(sub[,c[1]], sub[,c[2]], pch=k, col=cols[k])
            Plot_ConvexHull(sub[,c[1]], sub[,c[2]], lcol=cols[k], lwd=2)
         }
         legend('topleft', legend=c("Eoc", "Pal", "Cret"), lty=1, col=cols)
         dev.off()
         
         # Make 2-pulses clumped figure
         eoc<-subset(nmds.df, ageID<=time.splits[3])
         p2<-subset(nmds.df, ageID>=time.splits[3] & ageID <=time.splits[2])
         p1<-subset(nmds.df, ageID>=time.splits[2] & ageID <=time.splits[1])
         cret<-subset(nmds.df, ageID>=time.splits[1])
         ages.grouped<-list(eoc, p2, p1, cret)
         
         if(output.type == "pdf") {pdf(file=paste('figures/nmds_all/', nmds.nam[i], '_', cc[1], '-', cc[2], '_pulses.pdf', sep=''), height=10, width=10, useDingbats = FALSE)}
         if(output.type == "jpg") {jpeg(filename=paste('figures/nmds_all/', nmds.nam[i], '_', cc[1], '-', cc[2], '_pulses.jpg', sep='')) }
         plot(nmds.df[,c[1]], nmds.df[,c[2]], type='n', xlab=paste('MDS', cc[1], sep=''), ylab=paste('MDS', cc[2], sep=''), main=paste(nmds.nam[i], cc[1], cc[2], 'grouped', sep='_'))
         cols<-rainbow(length(ages.grouped))
         # Step 2: loop through the dataset to add to plot
         for(k in 1:length(ages.grouped)) {
            sub<-data.frame(ages.grouped[k])
            points(sub[,c[1]], sub[,c[2]], pch=k, col=cols[k])
            Plot_ConvexHull(sub[,c[1]], sub[,c[2]], lcol=cols[k], lwd=2)
         }
         legend('bottomleft', legend=c("Eoc", "Pulse 2", "Pulse 1", "Cret"), lty=1, col=cols)
         dev.off()
         
      }
   }
}

#make the plots:
nmds.plots.all(nmds.list=nmds.list, time.splits=time.splits, output.type = "pdf")
nmds.plots.all(nmds.list=nmds.list, time.splits=time.splits, output.type = "jpg")


### ****OLD**** Create data frame for individual teeth that will plot NMDS by origination, extinction, extant, etc.

make.nmds.plotting.dfs<-function(time.splits) {
   
   df.list<-list()
   
   for(i in 1:12) { #for the first 12 entries in the nmds.list [ignores nmds.morph]
      #call the correct NMDS to run
      nmds<-nmds.list[[i]] 
      dataset<-unlist(strsplit(nmds.nam[i], '.', fixed=T))[2]
      
      # get ageID for each individual tooth and taxon name for each individual tooth
      ageID<-get(paste('morphdat.', dataset, sep=''))$AgeID
      taxa<-as.character(get(paste('morphdat.', dataset, sep=''))$Morphotype_name) #is as.character necessary - yes!
      toothID<-as.character(get(paste('morphdat.', dataset, sep=''))$ID)
      
      # Create dataframe to populate. Has columns: toothID (from morphdat$ID), ageID (from morphdat$AgeID)
      
      df<-data.frame(toothID=toothID, ageID=ageID, taxa=taxa, 
         c.exist=0, c.ext=0, c.orig=0, 
         p1.exist=0, p1.ext=0, p1.orig=0,
         p2.exist=0, p2.ext=0, p2.orig=0,
         e.exist=0, e.ext=0, e.orig=0)
      df<-cbind(df, nmds$points)
      
      #Give existance values to each tooth based on occurrance (not morphotype range)
      df$c.exist<-ifelse(df$ageID >= time.splits[1], 1, 0)
      df$p1.exist<-ifelse(df$ageID <=time.splits[1] & df$ageID >=time.splits[2], 1, 0)
      df$p2.exist<-ifelse(df$ageID <=time.splits[2] & df$ageID >=time.splits[3], 1, 0)
      df$e.exist<-ifelse(df$ageID <=time.splits[3], 1, 0)
      
      # get FAD and LAD for each taxon
      ad<-get(paste('ad.', dataset, sep=''))
      
      # Give morphotype FAD/LAD (origination and extinction) intervals. 
      for(j in 1:length(ad[,1])) {
         tax<-rownames(ad)[j]  #pick taxon to work on
         time.range<-ad[j,]  #look up the FAD/LAD for that particular taxon
         index<-which(taxa==tax)  #find all occurrances with that specific taxon's name. 
         #Cretaceous is >66 Ma
         df$c.ext[index] <- ifelse(time.range$lads >=time.splits[1], 1, 0)
         df$c.orig[index] <- ifelse(time.range$fads >=time.splits[1], 1, 0)
            #p1 is 66 to 62
         df$p1.ext[index] <- ifelse(time.range$lads <=time.splits[1] & time.range$lads >=time.splits[2], 1, 0)
         df$p1.orig[index] <- ifelse(time.range$fads <=time.splits[1] & time.range$fads >=time.splits[2], 1, 0)
         #p2 is from 62 to 55
         df$p2.ext[index] <- ifelse(time.range$lads <=time.splits[2] & time.range$lads >=time.splits[3], 1, 0)
         df$p2.orig[index] <- ifelse(time.range$fads <=time.splits[2] & time.range$fads >=time.splits[3], 1, 0)
         #eocene is <55 Ma
         df$e.ext[index] <- ifelse(time.range$lads <=time.splits[3] & time.range$lads >=time.splits[4], 1, 0)
         df$e.orig[index] <- ifelse(time.range$fads <=time.splits[3], 1, 0)
      }
      
      #add df to the df list for output
      df.list[[i]]<-df
   }
   
   names(df.list)<-nmds.nam[1:12]
   return(df.list)
}

nmds.plot.list<-make.nmds.plotting.dfs(time.splits=time.splits)


#check my work... 
df<-nmds.plot.list[[1]]
sum(ifelse(apply(df[4:15], 1, sum)==2 | apply(df[4:15], 1, sum)==3,0,1))  #should all be 2 or 3, and this should equal 0 (yay!)

### *** NEW *** Create data frame for individual teeth for plotting NMDS values ###
### Make new nmds.plot.list for figure that includes: triangle down, triangle up, diamond, and square... ###

# triangle down: originated in an earlier time bin, went extinct during this interaval. For Cretaceous, define all earliest appearing taxa as this. Others as square
# triangle up: originated in this bin, went extinct later
# diamond: originated earlier and went extinct later; for Eocene, define any that appear in latest sample as range through, others as squares
# square: originated and went extinct in this bin. 

### new plotting data frame function
time.splits<-c(72.5, 66, 60, 55, 43) #updated

### new plotting data frame function ###
time.splits<-c(72.5, 66, 60, 55, 43) #updated

new.nmds.plotting.dfs <- function(time.splits) {
   df.list <- list() #empty list to add the data frames to
   
   for(i in 1:12) {#for the first 12 entries in the nmds.list [ignores nmds.morph]
      #call the correct NMDS to look at
      nmds<-nmds.list[[i]] 
      dataset<-unlist(strsplit(nmds.nam[i], '.', fixed=T))[2]
      
      # get ageID for each individual tooth and taxon name for each individual tooth
      ageID<-get(paste('morphdat.', dataset, sep=''))$AgeID
      taxa<-as.character(get(paste('morphdat.', dataset, sep=''))$Morphotype_name) #is as.character necessary - yes!
      toothID<-as.character(get(paste('morphdat.', dataset, sep=''))$ID)
      
      #create data frame
      df<-data.frame(toothID=toothID, ageID=ageID, taxa=taxa, 
         c.up = 0, c.down = 0, c.diamond = 0, c.square = 0, 
         p1.up = 0, p1.down = 0, p1.diamond = 0, p1.square = 0,
         p2.up = 0, p2.down = 0, p2.diamond = 0, p2.square = 0,
         e.up = 0, e.down = 0, e.diamond = 0, e.square = 0)
      df<-cbind(df, nmds$points)
      
      # get FAD and LAD for each taxon
      ad<-get(paste('ad.', dataset, sep=''))
      
      # for each morphotype, evaluate which range values are true for each time slice
      
      for(j in 1:length(ad[,1])) {
         tax<-rownames(ad)[j]  #pick taxon to work on
         time.range<-ad[j,]  #look up the FAD/LAD for that particular taxon
         index<-which(taxa==tax)  #find all occurrances with that specific taxon's name. 
         
         #Cretaceous is >66 Ma (between time splits 1, 2)
         df$c.up[index] <- ifelse(time.range$fads >= time.splits[2] & time.range$fads <= time.splits[1] &      # originated in Cretaceous
               time.range$lads <= time.splits[2], 1, 0)                                                        # went extinct later
         df$c.down[index] <- ifelse(time.range$fads >= time.splits[1] &                                        # Originated before Cretaceous
               time.range$lads <= time.splits[1] & time.range$lads >= time.splits[2], 1, 0)                    # went extinct in Cretaceous
         df$c.diamond[index] <- ifelse(time.range$fads >= time.splits[1] &                                     # originated before time period
               time.range$lads <= time.splits[2], 1, 0)                                                        # went extinct after time period
         df$c.square[index] <- ifelse(time.range$fads >= time.splits[2] & time.range$fads <= time.splits[1] &  #originated in Cretaceous
               time.range$lads <= time.splits[1] & time.range$lads >= time.splits[2], 1, 0)                    # went extinct in Cretaceous
         
         #p1 is 66 to 60 (between time splits 2, 3)
         df$p1.up[index] <- ifelse(time.range$fads >= time.splits[3] & time.range$fads <= time.splits[2] &      # originated in P1
               time.range$lads <= time.splits[3], 1, 0)                                                         # went extinct later
         df$p1.down[index] <- ifelse(time.range$fads >= time.splits[2] &                                        # Originated before P1
               time.range$lads <= time.splits[2] & time.range$lads >= time.splits[3], 1, 0)                     # went extinct in P1
         df$p1.diamond[index] <- ifelse(time.range$fads >= time.splits[2] &                                     # originated before time period
               time.range$lads <= time.splits[3], 1, 0)                                                         # went extinct after time period
         df$p1.square[index] <- ifelse(time.range$fads >= time.splits[3] & time.range$fads <= time.splits[2] &  # originated in P1
               time.range$lads <= time.splits[2] & time.range$lads >= time.splits[3], 1, 0)                     # went extinct in P1
         
         
         #p2 is from 60 to 55 (between time splits 3, 4)
         df$p2.up[index] <- ifelse(time.range$fads >= time.splits[4] & time.range$fads <= time.splits[3] &      # originated in P2
               time.range$lads <= time.splits[4], 1, 0)                                                         # went extinct later
         df$p2.down[index] <- ifelse(time.range$fads >= time.splits[3] &                                        # Originated before P2
               time.range$lads <= time.splits[3] & time.range$lads >= time.splits[4], 1, 0)                     # went extinct in P2
         df$p2.diamond[index] <- ifelse(time.range$fads >= time.splits[3] &                                     # originated before time period
               time.range$lads <= time.splits[4], 1, 0)                                                         # went extinct after time period
         df$p2.square[index] <- ifelse(time.range$fads >= time.splits[4] & time.range$fads <= time.splits[3] &  # originated in P2
               time.range$lads <= time.splits[3] & time.range$lads >= time.splits[4], 1, 0)                     # went extinct in P2
         
         
         #eocene is <55 Ma (between time splits 4, 5)
         
         df$e.up[index] <- ifelse(time.range$fads >= time.splits[5] & time.range$fads <= time.splits[4] &      # originated in eocene
               time.range$lads <= time.splits[5], 1, 0)                                                         # went extinct later
         df$e.down[index] <- ifelse(time.range$fads >= time.splits[4] &                                        # Originated before eocene
               time.range$lads <= time.splits[4] & time.range$lads >= time.splits[5], 1, 0)                     # went extinct in eocene
         df$e.diamond[index] <- ifelse(time.range$fads >= time.splits[4] &                                     # originated before time period
               time.range$lads <= time.splits[5], 1, 0)                                                         # went extinct after time period
         df$e.square[index] <- ifelse(time.range$fads >= time.splits[5] & time.range$fads <= time.splits[4] &  # originated in eocene
               time.range$lads <= time.splits[4] & time.range$lads >= time.splits[5], 1, 0)                     # went extinct in eocene
         
         # special case of originating and going extinct in youngest time bin, force e.up = 1
         if(time.range$fads == time.range$lads & time.range$fads <= time.splits[5]) {df$e.up[index] <- 1} 
         
      }
      
      ### Clean up the df output so that only values of 1 are for within time range that the individual fossil tooth occurred
      ## Cretaceous cleanup: Make all teeth with ageID < 66 c.xx values = 0
      non.cret <- which(df$ageID < time.splits[2])
      df$c.up[non.cret] <- 0
      df$c.down[non.cret] <- 0
      df$c.diamond[non.cret] <- 0
      df$c.square[non.cret] <- 0
      
      ## P1 cleanup - make all teeth with ageID > 66 or < 60 go away 
      non.p1 <- which(df$ageID > time.splits[2] | df$ageID < time.splits[3])
      df$p1.up[non.p1] <- 0
      df$p1.down[non.p1] <- 0
      df$p1.diamond[non.p1] <- 0
      df$p1.square[non.p1] <- 0
      
      ## P2 cleanup - make all teeth with ageID > 60 or < 55 go away 
      non.p2 <- which(df$ageID > time.splits[3] | df$ageID < time.splits[4])
      df$p2.up[non.p2] <- 0
      df$p2.down[non.p2] <- 0
      df$p2.diamond[non.p2] <- 0
      df$p2.square[non.p2] <- 0
      
      ## Eocene cleanup: Make all teeth with ageID > 55 values = 0
      non.eoc <- which(df$ageID > time.splits[4])
      df$e.up[non.eoc] <- 0
      df$e.down[non.eoc] <- 0
      df$e.diamond[non.eoc] <- 0
      df$e.square[non.eoc] <- 0
      
      #add df to the df.list for output
      df.list[[i]]<-df
      
   }
   #give df.list names and return the function
   names(df.list)<-nmds.nam[1:12]
   return(df.list)
   
}


nmds.plot.list.new <- new.nmds.plotting.dfs(time.splits)

# check my work
df <- nmds.plot.list.new$NMDS3.con

sum(apply(df[, 16:19],1, sum) > 1) #check eocene, seems to be good (should be = 0)
sum(apply(df[, 12:15],1, sum) >1 ) #check P2, seems to be good (should be = 0)
sum(apply(df[, 8:11],1, sum) > 1) #check P1
sum(apply(df[, 4:7],1, sum) > 1) #check cretaceous
sum(apply(df[, 4:19], 1, sum) == 0) #check if there are any which do not appear at all (0 values) - this should also be 0 at the end

# which(apply(df[, 4:19], 1, sum) ==0) #find the ones that have no spot in this. 
# # the function misses one where fad and lad are time 5 (aka originated after last point in the Eocene. 
#in a perfect world, I think we're actually set here..... 



##### (b) Make NMDS 4-panel figure for each of the 3 datasets by 4 NMDS options (12 figures total once I loop...) #####

orig.ext.nmds.plots<-function(orig.color = 'darkgreen', ext.color = 'orange', 
   exist.color = 'gray40', prior.color = 'gray40', 
   output.to.file=TRUE) {
   for(i in 1:length(nmds.plot.list)) {
   
      df<-nmds.plot.list[[i]]
      nam.df<-names(nmds.plot.list)[i]
      dims.df<-dim(df)
      nmds.pts.all<-df[,16:dims.df[2]]  #just the MDS points
      
      dims<-grep("MDS", names(nmds.pts.all))         # call all MDS columns
      combs<-t(combn(dims,2))                    # combinations for each
      
      for(j in 1:length(combs[,1])) {
         
         dims.to.plot<-combs[j,]
         nmds.pts<-nmds.pts.all[,dims.to.plot]
         plotname<-paste(nam.df, '_MDS', dims.to.plot[1], '-MDS', dims.to.plot[2], sep='')
         
         #start output to file 
         if(output.to.file==TRUE) {
            jpeg(filename=paste('figures/orig-ext-nmds/', plotname, '.jpg', sep=''), 
               width=10, height=10, units='in', res=150)
         }
         
         par(mfrow=c(2,2)) #2 rows, 2 columns of figures
         
         ### 1. Cretaceous plot: 
         plottitle<-paste(plotname, '_Cretaceous', sep='')
         plot(nmds.pts, type='n', main=plottitle)  #make empty plot
         
         #Previous time bin points (n/a for Cretaceous)

         # add convex hulls for originated and extinct teeth before plotting points: 
         Plot_ConvexHull(xcoord=nmds.pts[which(df$c.orig==1),], lcolor=orig.color, lwd=1)
         Plot_ConvexHull(xcoord=nmds.pts[which(df$c.ext==1),], lcolor=ext.color, lwd=1)
         
         # Current time bin points
         #teeth that went extinct, and were also extant during time P1 (outlines black points)
         points(nmds.pts[which(df$c.ext==1&df$c.exist==1),], pch=16, cex=1.8, col=ext.color)
         #add extant teeth
         points(nmds.pts[which(df$c.exist==1),], pch=16, col=exist.color)
         #add origination teeth present during the time bin
         points(nmds.pts[which(df$c.orig==1&df$c.exist==1),], pch=16, col=orig.color)
         

         ### 2. plot P1:
         plottitle<-paste(plotname, '_Pulse1', sep='')
         plot(nmds.pts, type='n', main=plottitle)  #make empty plot
         
         # Previous time bin
         points(nmds.pts[which(df$c.exist==1),], col='lightgray', pch=16, cex=0.5)
         Plot_ConvexHull(xcoord=nmds.pts[which(df$c.exist==1),], border.line = F, lcolor=prior.color, lwd=0, 
            shade=T, scolor=adjustcolor(prior.color,alpha.f = 0.25))
         
         # add convex hulls for origination and extinction prior to points
         Plot_ConvexHull(xcoord=nmds.pts[which(df$p1.orig==1),], lcolor=orig.color, lwd=1)
         Plot_ConvexHull(xcoord=nmds.pts[which(df$p1.ext==1),], lcolor=ext.color, lwd=1)
         
         # Current time bin points
         #teeth that went extinct, and were also extant during time P1 (outlines black points)
         points(nmds.pts[which(df$p1.ext==1&df$p1.exist==1),], pch=16, cex=1.8, col=ext.color)
         #add extant points
         points(nmds.pts[which(df$p1.exist==1),], pch=16, col=exist.color)
         #add origination points present during the time bin
         points(nmds.pts[which(df$p1.orig==1&df$p1.exist==1),], pch=16, col=orig.color)
         
         
         ### 3. plot P2: 
         plottitle<-paste(plotname, '_Pulse2', sep='')
         plot(nmds.pts, type='n', main=plottitle) 
         
         # Previous time bin
         points(nmds.pts[which(df$p1.exist==1),], col='lightgray', pch=16, cex=0.5)
         Plot_ConvexHull(xcoord=nmds.pts[which(df$p1.exist==1),], border.line = F, lcolor=prior.color, lwd=0, 
            shade=T, scolor=adjustcolor(prior.color,alpha.f = 0.25))
         
         # plot convex hulls prior to points
         Plot_ConvexHull(xcoord=nmds.pts[which(df$p2.orig==1),], lcolor=orig.color, lwd=1)
         Plot_ConvexHull(xcoord=nmds.pts[which(df$p2.ext==1),], lcolor=ext.color, lwd=1)
         
         # Current time bin points
         #teeth that went extinct, and were also extant during time P1 (outlines black points)
         points(nmds.pts[which(df$p2.ext==1&df$p2.exist==1),], pch=16, cex=1.8, col=ext.color)
         #add extant points
         points(nmds.pts[which(df$p2.exist==1),], pch=16, col=exist.color)
         #add origination points present during the time bin
         points(nmds.pts[which(df$p2.orig==1&df$p2.exist==1),], pch=16, col=orig.color)
         
         ### 4. plot eocene:
         plottitle<-paste(plotname, '_Eocene', sep='')
         plot(nmds.pts, type='n', main=plottitle) 
         
         # Previous time bin
         points(nmds.pts[which(df$p2.exist==1),], col='lightgray', pch=16, cex=0.5)
         Plot_ConvexHull(xcoord=nmds.pts[which(df$p2.exist==1),], border.line = F, lcolor=prior.color, lwd=0, 
            shade=T, scolor=adjustcolor(prior.color,alpha.f = 0.25))
         
         # plot convex hulls prior to points
         Plot_ConvexHull(xcoord=nmds.pts[which(df$e.orig==1),], lcolor=orig.color, lwd=1)
         Plot_ConvexHull(xcoord=nmds.pts[which(df$e.ext==1),], lcolor=ext.color, lwd=1)
         
         # Current time bin points
         #teeth that went extinct, and were also extant during time P1 (outlines black points)
         points(nmds.pts[which(df$e.ext==1&df$e.exist==1),], pch=16, cex=1.8, col=ext.color)
         #add extant points
         points(nmds.pts[which(df$e.exist==1),], pch=16, col=exist.color)
         #add origination points present during the time bin
         points(nmds.pts[which(df$e.orig==1&df$e.exist==1),], pch=16, col=orig.color)
      
         # end output to file
         if(output.to.file==TRUE) {
            dev.off()
         }
         
      } #close J loop
      
   } #close i loop
} #end function

# generate the figures
orig.ext.nmds.plots(output.to.file=T, exist.color = 'gray60') 


##### (c) Make NMDS Orig/Ext Figure for papers/talks: ##### 

### 1. Origination/extinction plots ###

# Step 1: Define colors and objects to plot
orig.color = 'darkgreen'
ext.color = 'orange'
exist.color = 'gray40'
prior.color = 'gray40'

df<-nmds.plot.list$NMDS3.con  #pull the appropriate NMDS plot

#points to plot: MDS1 and MDS2 (cols 16 and 17)
nmds.pts<-df[,16:17]

# Step 2: Define output location 
# jpeg(filename='figures/NMDS_Figure_nmds2con.jpg', 
#    width=10, height=10, units='in', res=150)

pdf(file='figures/NMDS_Figure_nmds3con.pdf', 
   width=10, height=10, useDingbats = FALSE)

# Step 3: Construct the plots
par(mfrow=c(2,2)) #2 rows, 2 columns of figures

### 1. Cretaceous plot: 
plot(nmds.pts, type='n', main='Cretaceous (72-66 Ma)')  #make empty plot

#Previous time bin points (n/a for Cretaceous)

# add convex hulls for originated and extinct teeth before plotting points: 
Plot_ConvexHull(xcoord=nmds.pts[which(df$c.orig==1),], lcolor=orig.color, lwd=1)
Plot_ConvexHull(xcoord=nmds.pts[which(df$c.ext==1),], lcolor=ext.color, lwd=1)

# Current time bin points
#teeth that went extinct, and were also extant during time P1 (outlines black points)
points(nmds.pts[which(df$c.ext==1&df$c.exist==1),], pch=16, cex=1.8, col=ext.color)
#add extant teeth
points(nmds.pts[which(df$c.exist==1),], pch=16, col=exist.color)
#add origination teeth present during the time bin
points(nmds.pts[which(df$c.orig==1&df$c.exist==1),], pch=16, col=orig.color)


### 2. plot P1:
plot(nmds.pts, type='n', main='Paleocene Pulse 1 (66-60 Ma)')  #make empty plot

# Previous time bin
points(nmds.pts[which(df$c.exist==1),], col='lightgray', pch=16, cex=0.5)
Plot_ConvexHull(xcoord=nmds.pts[which(df$c.exist==1),], border.line = F, lcolor=prior.color, lwd=0, 
   shade=T, scolor=adjustcolor(prior.color,alpha.f = 0.25))

# add convex hulls for origination and extinction prior to points
Plot_ConvexHull(xcoord=nmds.pts[which(df$p1.orig==1),], lcolor=orig.color, lwd=1)
Plot_ConvexHull(xcoord=nmds.pts[which(df$p1.ext==1),], lcolor=ext.color, lwd=1)

# Current time bin points
#teeth that went extinct, and were also extant during time P1 (outlines black points)
points(nmds.pts[which(df$p1.ext==1&df$p1.exist==1),], pch=16, cex=1.8, col=ext.color)
#add extant points
points(nmds.pts[which(df$p1.exist==1),], pch=16, col=exist.color)
#add origination points present during the time bin
points(nmds.pts[which(df$p1.orig==1&df$p1.exist==1),], pch=16, col=orig.color)


### 3. plot P2: 
plot(nmds.pts, type='n', main='Paleocene Pulse 2 (60-55 Ma)') 

# Previous time bin
points(nmds.pts[which(df$p1.exist==1),], col='lightgray', pch=16, cex=0.5)
Plot_ConvexHull(xcoord=nmds.pts[which(df$p1.exist==1),], border.line = F, lcolor=prior.color, lwd=0, 
   shade=T, scolor=adjustcolor(prior.color,alpha.f = 0.25))

# plot convex hulls prior to points
Plot_ConvexHull(xcoord=nmds.pts[which(df$p2.orig==1),], lcolor=orig.color, lwd=1)
Plot_ConvexHull(xcoord=nmds.pts[which(df$p2.ext==1),], lcolor=ext.color, lwd=1)

# Current time bin points
#teeth that went extinct, and were also extant during time P1 (outlines black points)
points(nmds.pts[which(df$p2.ext==1&df$p2.exist==1),], pch=16, cex=1.8, col=ext.color)
#add extant points
points(nmds.pts[which(df$p2.exist==1),], pch=16, col=exist.color)
#add origination points present during the time bin
points(nmds.pts[which(df$p2.orig==1&df$p2.exist==1),], pch=16, col=orig.color)

### 4. plot eocene:
plot(nmds.pts, type='n', main='Eocene (55-42 Ma)') 

# Previous time bin
points(nmds.pts[which(df$p2.exist==1),], col='lightgray', pch=16, cex=0.5)
Plot_ConvexHull(xcoord=nmds.pts[which(df$p2.exist==1),], border.line = F, lcolor=prior.color, lwd=0, 
   shade=T, scolor=adjustcolor(prior.color,alpha.f = 0.25))

# plot convex hulls prior to points
Plot_ConvexHull(xcoord=nmds.pts[which(df$e.orig==1),], lcolor=orig.color, lwd=1)
Plot_ConvexHull(xcoord=nmds.pts[which(df$e.ext==1),], lcolor=ext.color, lwd=1)

# Current time bin points
#teeth that went extinct, and were also extant during time P1 (outlines black points)
points(nmds.pts[which(df$e.ext==1&df$e.exist==1),], pch=16, cex=1.8, col=ext.color)
#add extant points
points(nmds.pts[which(df$e.exist==1),], pch=16, col=exist.color)
#add origination points present during the time bin
points(nmds.pts[which(df$e.orig==1&df$e.exist==1),], pch=16, col=orig.color)

# end output to file
dev.off()






##### (d) Make NMDS extant figure #####

cret.color <- 'royalblue4'
p1.color <- 'springgreen4'
p2.color <- 'sandybrown'
e.color <- 'firebrick'

### Setup the workspace:
pdf(file='figures/NMDS_Figure_nmds3con_extant_horiz.pdf', 
   width=16, height=4, useDingbats = FALSE)
par(mfrow=c(1,4)) #1 row, 4 columns for

### 1. Cretaceous plot: 
plot(nmds.pts, type='n', main='Cretaceous (72-66 Ma)')  #make empty plot
#add extant teeth
points(nmds.pts[which(df$c.exist==1),], pch=16, col=cret.color)
#plot convex hull
Plot_ConvexHull(xcoord=nmds.pts[which(df$c.exist==1),], border.line = F, lcolor=cret.color, lwd=0, 
   shade=T, scolor=adjustcolor(cret.color,alpha.f = 0.8))

### 2. P1 plot: 
plot(nmds.pts, type='n', main='Paleocene Pulse 1 (66-60 Ma)')  #make empty plot
#add extant teeth
points(nmds.pts[which(df$p1.exist==1),], pch=16, col=p1.color)
#plot convex hull
Plot_ConvexHull(xcoord=nmds.pts[which(df$p1.exist==1),], border.line = F, lcolor=p1.color, lwd=0, 
   shade=T, scolor=adjustcolor(p1.color,alpha.f = 0.8))

### 3. P2 plot:
plot(nmds.pts, type='n', main='Paleocene Pulse 2 (60-55 Ma)') 
#add extant teeth
points(nmds.pts[which(df$p2.exist==1),], pch=16, col=p2.color)
#plot convex hull
Plot_ConvexHull(xcoord=nmds.pts[which(df$p2.exist==1),], border.line = F, lcolor=p2.color, lwd=0, 
   shade=T, scolor=adjustcolor(p2.color,alpha.f = 0.8))

### 4. Eocene plot:
plot(nmds.pts, type='n', main='Eocene (55-42 Ma)') 
#add extant teeth
points(nmds.pts[which(df$e.exist==1),], pch=16, col=e.color)
#plot convex hull
Plot_ConvexHull(xcoord=nmds.pts[which(df$e.exist==1),], border.line = F, lcolor=e.color, lwd=0, 
   shade=T, scolor=adjustcolor(e.color,alpha.f = 0.8))

# end output to file
dev.off()

##### Calculate percent of novelty that persists past time period  #####


### Calculate percent of novelty that goes extinct in the same time period

# Make count DF from age datums (counting morphotypes not individuals)
# - originated during time period
# - originated and also went extinct during time period
# - originated but did not go extinct during time period

ad<-ad.con

## Cretaceous
# originated during time period: 48
sum(ifelse(which(ad$fads >= time.splits[1]),1,0))
ad[which(ad$fads >= time.splits[1]),]
# originated and went extinct during time period: 2
sum(ifelse(which(ad$fads >= time.splits[1] & ad$lads >= time.splits[1]),1,0))
ad[which(ad$fads >= time.splits[1] & ad$lads >= time.splits[1]),]
# originated and did not go extinct during time period: 46 
sum(ifelse(which(ad$fads >= time.splits[1] & ad$lads <= time.splits[1]),1,0))
ad[which(ad$fads >= time.splits[1] & ad$lads <= time.splits[1]),]

## Pulse 1
# originated during time period: 34
sum(ifelse(which(ad$fads <= time.splits[1] & ad$fads >= time.splits[2]),1,0))
ad[which(ad$fads <= time.splits[1] & ad$fads >= time.splits[2]),]
# originated and went extinct during time period: 6
sum(ifelse(which(ad$fads <= time.splits[1] & ad$fads >= time.splits[2] & ad$lads >= time.splits[2]),1,0))
ad[which(ad$fads <= time.splits[1] & ad$fads >= time.splits[2] & ad$lads >= time.splits[2]),]
# originated and did not go extinct during time period: 28
sum(ifelse(which(ad$fads <= time.splits[1] & ad$fads >= time.splits[2] & ad$lads <= time.splits[2]),1,0))
ad[which(ad$fads <= time.splits[1] & ad$fads >= time.splits[2] & ad$lads <= time.splits[2]),]

## Pulse 2
# originated during time period: 39
sum(ifelse(which(ad$fads <= time.splits[2] & ad$fads >= time.splits[3]),1,0))
ad[which(ad$fads <= time.splits[2] & ad$fads >= time.splits[3]),]
# originated and went extinct during time period: 11
sum(ifelse(which(ad$fads <= time.splits[2] & ad$fads >= time.splits[3] & ad$lads >= time.splits[3]),1,0))
ad[which(ad$fads <= time.splits[2] & ad$fads >= time.splits[3] & ad$lads >= time.splits[3]),]
# originated and did not go extinct during time period: 28
sum(ifelse(which(ad$fads <= time.splits[2] & ad$fads >= time.splits[3] & ad$lads <= time.splits[3]),1,0))
ad[which(ad$fads <= time.splits[2] & ad$fads >= time.splits[3] & ad$lads <= time.splits[3]),]

## Eocene
# originated during time period: 14
sum(ifelse(which(ad$fads <= time.splits[3] & ad$fads >= time.splits[4]),1,0))
ad[which(ad$fads <= time.splits[3] & ad$fads >= time.splits[4]),]
# originated and went extinct during time period: 6
sum(ifelse(which(ad$fads <= time.splits[3] & ad$fads >= time.splits[4] & ad$lads >= time.splits[4]),1,0))
ad[which(ad$fads <= time.splits[3] & ad$fads >= time.splits[4] & ad$lads >= time.splits[4]),]
# originated and did not go extinct during time period: 8
sum(ifelse(which(ad$fads <= time.splits[3] & ad$fads >= time.splits[4] & ad$lads <= time.splits[4]),1,0))
ad[which(ad$fads <= time.splits[3] & ad$fads >= time.splits[4] & ad$lads <= time.splits[4]),]


##### Calculate and report stress values #####
lapply(nmds.list, function(x) return(x$stress))
# $NMDS2.orig
# [1] 0.1508316
# 
# $NMDS3.orig
# [1] 0.109862
# 
# $NMDS4.orig
# [1] 0.08433382
# 
# $NMDS5.orig
# [1] 0.07049754
# 
# $NMDS2.con
# [1] 0.1505748
# 
# $NMDS3.con
# [1] 0.110008
# 
# $NMDS4.con
# [1] 0.08434291
# 
# $NMDS5.con
# [1] 0.07051625
# 
# $NMDS2.lib
# [1] 0.1506772
# 
# $NMDS3.lib
# [1] 0.1096905
# 
# $NMDS4.lib
# [1] 0.08397646
# 
# $NMDS5.lib
# [1] 0.07047626
# 
# $NMDS2.morph
# [1] 0.1731968
# 
# $NMDS3.morph
# [1] 0.1264346
# 
# $NMDS4.morph
# [1] 0.09802952
# 
# $NMDS5.morph
# [1] 0.07903778
# 
