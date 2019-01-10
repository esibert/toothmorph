#################################################
#   functions written for toothmorph analyses   
# 
#  1. call.traits.csvs()
#     function crawls through the working directory to find all csv files with 
#     the filename Trait* (capital specified... this makes things easier... 
#  2. toothdat.cleanup(toothdat, fix_dat = FALSE, sortby = 'age-obj')
#     coerces the morphological dataset into an object that can be used 
#     in the distances.clust function, below. 
#  3. distances.clust(morph, traits, weights) 
#     calculates the pairwise distance between individual teeth, by timeslice.
#     requires package doParallel
#     input: morph matrix, list of traits distance matrices, vector of weights
#     output: data frame with row IDs, tooth ID's, distance sums, distance averages, etc.
#  4. distmat.fn(distpairs, type='avg')
#     input: distpairs data frame from distances.clust function
#     output: distance matrix (not distance object)
#  5. build.strat.obj(morphdat)
#     input: morphdat data frame
#     output: strat column object including occurrance data, taxon names, and ages
#  6. combine.samples(morphdat, combines)
#     input: morphdat data frame and list of vectors defining ages of samples to bin. Each list object is one bin
#     output: morphdat data frame with age values for the average age of each combined age bin. 
#  7. Plot_ConvexHull(xcoord, ycoord, border.line, lcolor, lwd, shade, scolor)
#     input: x and y coordinents, a color for the lines, and a color for shading
#     output: adds convex hulls to NMDS plots. 
#  8. rangechart3(range.sc)
#     Generates a range chart figure for a given stratigraphic column (sc) object
#     Includes a number of graphical parameters that can be tewaked
#     input: a full stratigraphic object OR a counts matrix and depths vector.
#  9. make.inp(x, filename, header=NULL) 
#     Makes an .inp file for use with RMark. the input 'x' is the table of counts 
#        from strat column object sc$counts. 
#     Requires filename (and [relative] path) input to run

##### 1. call.traits.csvs() #####
call.traits.csvs<-function() {
   #Find the relevant Trait.csv files (note, must be capitalized Trait not trait)
   files <- list.files(pattern="Trait*", recursive=T)  #Locations of trait connectivity matrices, called Trait*.csv
   
   # call in .csv distance matrix files, normalize, and build traits list
   traits<-list()                        #Make empty list for traits list
   
   for (i in 1:length(files)) {           #Start loop to call in .csv files
      foo<-read.csv(files[i], header=F)   #read in the file as data frame
      norm<-foo/max(foo)                  #normalize the denticle traits
      # assign(as.character(split), norm)   #Give name to normalized matrix
      traits[[i]]<-norm                    #Generate list of all trait matrices
   } 
   
   # Name the traits matrices within the traits list
   split<-as.character(strsplit(files, ".csv"))   #call name of text before '.csv'
   split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))  
   split<-sapply(split, split_path)[1,]
   names(split)<-split 
   names(traits)<-split  #Note that the list traits is not capitalized
   
   #output
   return(traits)
}

##### 2. toothdat.cleanup(toothdat, fix_dat=FALSE) #####
# fix_dat = TRUE means add 1 to 0/1 presence absence. If presence/absence traits are 1/2, fix_dat = FALSE
# sortby options include: 
#     'age-obj' (sort by object ID that includes sample ID (age) and object ID)
#     'age' - sort by ageID
#     'morph' - sort by morphotype ID by row
#     'original' - return to original spreadsheet order
toothdat.cleanup<-function(toothdat, fix_dat=FALSE, sortby='age-obj') {
   teeth_all<-subset(toothdat, A==1)  #Clean up data to be just the teeth (TraitA = 1)
   teeth<-teeth_all[complete.cases(teeth_all),]  #clean up data to include just teeth that have been described
   teeth<-subset(teeth, B!=4) #Get rid of poor quality teeth if anything is still there... 
   if(fix_dat==TRUE) {  #for spreadsheets I compiled which don't use the numerical values in the published figure
      teeth$K1 <- teeth$K1 + 1   #flange presence/absence trait
      teeth$O <- teeth$O + 1 #Root presence/absence
   }
   teeth.g106<-subset(teeth, WID>=100)  #Only ones with at least one dimension >100 um; note that here wid was by definition the smaller dimension
   teeth.nowid<-subset(teeth, WID==0 & AR!=0) # Samples where I've put in AR (including AR=99) but no wid/len measurements
   teeth.dat<-rbind(teeth.g106, teeth.nowid)
   
   #sort the samples... 
   if(sortby=='age-obj') {
      teeth.dat<-teeth.dat[order(teeth.dat[,4]), ]  #age/object order (!) 
   }
   else if(sortby == 'age') {
      teeth.dat<-teeth.dat[order(teeth.dat[,2]), ] } #age-only order, objects may be scrambled
   else if(sortby == 'morph') {
      teeth.dat<-teeth.dat[order(teeth.dat[,3]), ] } #sort by morphotype number
   else if(sortby == 'original') {
      teeth.dat <- teeth.dat[order(as.numeric(row.names(teeth.dat))),] }  #put them in original (csv) order again
   else teeth.dat<-teeth.dat  #no ordering
   
   teeth.dat$AR[teeth.dat$AR == 99]<-0  #reset dummy AR to 0 for diststances.clust analyses. 
   return(teeth.dat)
}

##### 3. distances.clust(morph, traits, weights) #####
distances.clust<-function(morph, traits, weights, coresFree=2) {
   #call traits distance matrices directly from the working directory
   if(missing(traits)) {
      traits<-call.traits.csvs()
   }
   
   #assign equal weight to each trait if no weights
   if(missing(weights)) {
      weights<-rep(1, length(traits)+3)
   }   
   
   # Create matrix of just traits, so that names can be maintained in the function properly
   morph.mat<-morph[,7:length(morph)]
   
   # Create sets of pairs: 
   #species<-teeth$Obj_num  # length should match number of teeth in analysis - uses object numbers from subset
   species<-c(1:length(morph.mat[,4])) #doesn't matter what column this calls, b/c just giving integers 1-n
   combs.rows<-t(combn(species,2)) #2xn Matrix of all possible pairwise combinations of species
   toothID<-as.numeric(rownames(morph))
   objID<-as.character(morph$ID)
   combs.toothID<-t(combn(toothID,2))
   
   #Set up cores and cluster to run loop
   cores<-detectCores()
   cl<-makeCluster(cores-coresFree)  #detect the number of cores and leave some open
   registerDoParallel(cl)
   
   #Set up parallel loop: 
   dat<-foreach(i=1:length(combs.rows[,1]), .combine='rbind') %dopar% { #for all pairwise comparisons, in parallel
      cc<-combs.rows[i,]  #Call the pairwise comparison to look at 
      Traitscomb<-as.vector(c())     #create empty vector to fill with comparissons of whichever 2 species are specified for the loop
      
      for(j in 1:length(traits)) {   #For however many traits there are...
         w<-weights[j]              #trait weight value
         tt<-data.frame(traits[j])    #calls relevant trait matrix and turns it into dataframe for easier handling
         foo<-as.numeric(tt[morph.mat[cc[1],j], morph.mat[cc[2],j]]) #calls the distance value of trait[i] by looking up first Dent A's trait[i] and then Dent B's Trait[i]
         foo.w<-foo*w
         Traitscomb<-c(Traitscomb, foo.w)   #Puts the trait[i] distance value in the comparisson vector
      }
      
      #add in continuous distances calculated based on **normalized** input vectors;
      #Must be last columns in trait matrix (cannot be mixed into the discrete traits)
      #Called: AR, LEN, WID
      
      if(morph$AR[cc[1]] == 0 ) { #| morph$AR[cc[2]] == 0 ) {
         Traitscomb<-Traitscomb } #else {
      else if(morph$AR[cc[2]] == 0) {
         Traitscomb<-Traitscomb }
      else {
         ar.dist<-dist(morph$AR) # distance calculation of aspect ratio distances for teeth considered in this analysis
         ar.dist<-ar.dist/max(ar.dist) #normalized
         ar<-ar.dist[i]*weights[length(traits)+1]
         Traitscomb<-c(Traitscomb, ar) }
      
      if(morph$LEN[cc[1]] == 0) {
         Traitscomb<-Traitscomb }
      else if(morph$LEN[cc[2]] == 0) {
         Traitscomb<-Traitscomb }
      else {
         len.dist<-dist(morph$LEN) # distance calculation of length differences for teeth considered in this analysis
         len.dist<-len.dist/max(len.dist) #normalized
         len<-len.dist[i]*weights[length(traits)+2]
         Traitscomb<-c(Traitscomb, len) }
      
      if(morph$WID[cc[1]] == 0) {
         Traitscomb<-Traitscomb }
      else if(morph$WID[cc[2]] == 0) {
         Traitscomb<-Traitscomb }
      else {
         wid.dist<-dist(morph$WID) # distance calculation of width differences for teeth considered in this analysis
         wid.dist<-wid.dist/max(wid.dist) #normalized
         wid<-wid.dist[i]*weights[length(traits)+3]
         Traitscomb<-c(Traitscomb, wid) }
      
      Traitscomb<-unlist(Traitscomb)
      Traitscomb<-as.numeric(Traitscomb)
      Dist<-sum(Traitscomb)  #Adds up distances calculated in inner loop
      Dist.avg<-mean(Traitscomb) #Mean of distance between 2 species
      Num.traits<-length(Traitscomb)
      to.dat<-data.frame(cc[1], cc[2], combs.toothID[i,1], combs.toothID[i,2], Dist, Dist.avg, Num.traits, objID[cc[1]], objID[cc[2]])
      to.dat #output of foreach loop
   }
   
   #stop cluster and remove traces of it (otehrwise the comptuer hangs)
   stopCluster(cl)
   rm(cl, cores)
   
   #Get final data table to put into distmat
   #col<-length(unlist(dat[1]))
   #dist.df<-data.frame(matrix(unlist(dat), ncol=col, byrow=T))
   dist.df<-data.frame(dat)
   names(dist.df)<-c("1", "2", "ToothA","ToothB","dist.sum","dist.avg", "traits.length", 'objID.A', 'objID.B')
   return(dist.df)
}

##### 4. distmat(distpairs) #####
# can be type='sum' or type='avg' for different kinds of distances. 
distmat.fn<-function(distpairs, type='avg') {
   # Make Distance Matrix
   objs<-unique(c(distpairs[,3], distpairs[,4]))  #make a list of all the unique objects considered
   dist.mat<-matrix(data=0, nrow=length(objs), ncol=length(objs))  #matrix of 0-values of dimensions [objs x objs]
   
   if(type=='avg') {
      for (i in 1:length(distpairs[,3])) {    
         dat<-as.numeric(distpairs[i,])        # call the first row of the data frame and coerce it into numeric values
         dist.mat[dat[1], dat[2]]<-dat[6] # 6th column is dist.avg; [1] and [2] are matrix locations 
         dist.mat[dat[2], dat[1]]<-dat[6]
      }
   }
   else if(type == 'sum') {
      for (i in 1:length(distpairs[,3])) {    
         dat<-as.numeric(distpairs[i,])        # call the first row of the data frame and coerce it into numeric values
         dist.mat[dat[1], dat[2]]<-dat[5] # 5th column is dist.sum; [1] and [2] are matrix locations
         dist.mat[dat[2], dat[1]]<-dat[5]
      }
   }
   
   rownames(dist.mat)<-objs
   colnames(dist.mat)<-objs
   return(dist.mat)
}

##### 5. build.strat.obj(morphdat) #####
build.strat.obj<-function(morphdat) {
   #define unique age bins for the samples
   AgeID.unique<-unique(morphdat$AgeID)
   #call each unique taxa name
   taxa<-levels(morphdat$Morphotype_name)[2:length(levels(morphdat$Morphotype_name))] #excludes the '' blank; should have 136 items
   #build the strat column object for use with package stratigraph, and as storage for ranges
   sub<-data.frame(morphdat$AgeID, morphdat$Morphotype_name)
   ranges<-table(sub, exclude='')  #make occurrance table#
   range.sc<-strat.column(ranges, absolute.ages=AgeID.unique, taxa=taxa)  #make "strat column" object
   #return the strat column object for manipulation
   return(range.sc)
}


##### 6. combine.samples(morphdat, combines) #####
# inputs: morphdat; combines - list of vectors defining samples in each bin 
# Note that the combines vectors are currently defined within the function for 
#   clarity of this code, but are unique to this specific dataset
combine.samples<-function(morphdat, combines) {
   
   if(missing(combines)) {
      ## this was done manually for this paper
      # Define samples to bin
      c09<-c(57.12, 57.34, 57.88, 58.14, 58.4, 58.66)  #104, 105, 106, 107, 108, 109     #52, 53, 54, 55, 56, 57  #note not all of the slides were described, due to computer errors
      c10<-c(58.90, 59.13, 59.36, 59.59, 59.82, 60.05) #110, 111, 112, 113, 114, 115     #58, 59, 60, 61, 62, 63
      c11<-c(60.28, 60.5, 60.7, 60.9, 61.1)       #116, 117, 118, 119, 120       #64, 65, 66, 67, 68
      c12<-c(61.53, 61.98, 62.20, 62.40, 62.60)   #122, 123, 124, 125, 126       #70, 71, 72, 73, 74
      c13<-c(62.81, 63.01, 63.22, 63.42, 63.61)   #127, 128, 129, 130, 131       #75, 76, 77, 78, 79
      c14<-c(63.8, 64, 64.4, 64.59, 64.77)        #132, 133, 134, 135, 136       #80, 81, 82, 83, 84
      c15<-c(64.94, 65.3, 65.5, 65.7, 65.9)       #137, 139, 140, 141, 142       #85, 87, 88, 89, 90
      c16<-c(66.31,66.51, 66.72, 66.92, 67.12)    #143, 144, 145, 146, 147       #91, 92, 93, 94, 95
      c17<-c(67.31, 67.51, 67.71, 67.90, 68.33)   #148, 149, 150, 151, 153       #96, 97, 98, 99, 101
      c18<-c(68.78, 69.00, 69.20, 69.41, 69.61)   #155, 156, 157, 158, 159       #103, 104, 105, 106, 107        
      c19<-c(70.31, 70.53, 70.75, 70.97, 71.19)   #161, 162, 163, 164, 165       #109, 110, 111, 112, 113   
      c20<-c(71.40, 71.62, 71.83, 72.05, 72.26)   #166, 167, 168, 169, 170       #114, 115, 116, 117, 118   
      c21<-c(72.47, 72.68, 72.89, 73.10, 73.30)   #171, 172, 173, 174, 175       #119, 120, 121, 122, 123   
      combines<-list(c09, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20, c21) 
   }
   
   #substitutions
   for(i in 1:length(combines)) {
      c.sub<-combines[[i]]
      morphdat$AgeID<-ifelse(morphdat$AgeID %in% c.sub, mean(c.sub), morphdat$AgeID)
   }
   
   return(morphdat)
}

##### 7. Plot_ConvexHull(xcoord, ycoord, lcolor, lwd) #####
# plots a convex hull around a field of points in a plot. Used for NMDS plots
Plot_ConvexHull<-function(xcoord, ycoord, border.line=TRUE, lcolor, lwd, shade=FALSE, scolor=NULL, ...){
   if(missing(ycoord)) {
      ycoord<-xcoord[,2]
      xcoord<-xcoord[,1]
   }
   hpts <- chull(x = xcoord, y = ycoord)
   hpts <- c(hpts, hpts[1])
   if(border.line==T) {
      lines(xcoord[hpts], ycoord[hpts], col = lcolor, lwd=lwd, ...)
   }
   if(shade==T) {
      if(missing(scolor)) {scolor = lcolor} 
      polygon(xcoord[hpts], ycoord[hpts], col=adjustcolor(scolor, alpha.f=0.5), border=NA)
   }
}  

##### 8. rangechart3(range.sc) #####
# riff off of rangechart() from stratigraph package to suit my particular needs
rangechart3<-function (x = NULL, counts = NULL, depths = NULL, sample.labels = NULL, 
   taxa = NULL, short.names = NULL, higher.grp = NULL, tax.cat = NULL, 
   reorder = NULL, plot.points = FALSE, plot.depths.increasing.down = TRUE, 
   llwd = 2, llcol = gray(0.5), llty = 1, cex.xaxis = 0.8, cex.yaxis = 1, cex.points = 1, pch.points = 1, 
   col.points = 'black', colors.vec = NULL, baselines = FALSE, 
   xaxis.labels = c('names', 'numbers', 'alphanum'), return.xaxis = FALSE, 
   legend = TRUE, legend.values = NULL, legend.loc = NULL, legend.horiz = FALSE, 
   legend.title = NULL, legend.bg = 'white', large.size = 1, count.group = FALSE, ...) 
{
   if (is.strat.column(x)) {
      counts <- x$counts
      if (is.null(depths)) 
         depths <- x$depths
      if (is.null(taxa)) 
         taxa <- x$taxa
      if (is.null(short.names)) 
         short.names <- x$short.names
      if (is.null(higher.grp)) 
         higher.grp <- x$higher.grp
      if (is.null(tax.cat)) 
         tax.cat <- x$tax.cat
   }
   if (is.data.frame(counts) || is.matrix(counts)) {
      counts <- apply(counts, 2, "as.numeric")
   }
   else if (is.character(counts) && length(counts) == 1) {
      counts <- read.csv(counts, header = TRUE, skip = 0, colClasses = "")
      counts <- apply(counts, 2, "as.numeric")
   }
   else {
      stop("argument to counts not understood")
   }
   if (is.null(depths)) {
      warning("no depths provided; plotting samples at regular intervals")
      depths <- 1:nrow(counts)
   }
   depths <- as.numeric(depths)
   if (length(depths) != nrow(counts)) {
      stop(paste(length(depths), " depths, and ", nrow(counts), 
         " rows in the count matrix.", sep = ""))
   }
   if (sum(is.na(counts)) > 0) {
      warning(paste(sum(is.na(counts)), "missing values in count matrix replaced with zeros"))
      counts[is.na(counts)] <- 0
   }
   emptycols <- !(colSums(counts, na.rm = TRUE) > 0)
   if (any(emptycols)) {
      counts <- counts[, !emptycols]
      tax.cat <- tax.cat[!emptycols]
      taxa <- taxa[!emptycols]
      warning(paste(sum(emptycols), "columns with zero counts at all levels removed"))
   }
   
   ######alphanumeric rewrite
   if(xaxis.labels == 'alphanum') {
      colnames(counts) <- 1:length(colnames(counts))
   }
   
   if (!is.null(reorder)) {
      if (length(reorder) == 1 && is.character(reorder)) {
         funny <- function(x) return((1:length(x))[x > 0])
         if (pmatch(reorder, "fad.by.category", nomatch = FALSE)) {
            fads <- depths[as.numeric(lapply(apply(counts, 
               2, funny), max))]
            reorder.vect <- sort(fads, decreasing = TRUE, 
               index.return = TRUE)$ix
            counts <- counts[, reorder.vect]
            reorder.vect <- sort(as.character(taxa), index.return = TRUE)$ix
            counts <- counts[, reorder.vect]
         }
         else if (pmatch(reorder, "lad.by.category", nomatch = FALSE)) {
            lads <- depths[as.numeric(lapply(apply(counts, 
               2, funny), min))]
            reorder.vect <- sort(lads, decreasing = TRUE, 
               index.return = TRUE)$ix
            counts <- counts[, reorder.vect]
            reorder.vect <- sort(as.character(taxa), index.return = TRUE)$ix
            counts <- counts[, reorder.vect]
         }
         else if (pmatch(reorder, "lad.by.fad", nomatch = FALSE)) {
            lads <- depths[as.numeric(lapply(apply(counts, 
               2, funny), min))]
            reorder.vect <- sort(lads, decreasing = TRUE, 
               index.return = TRUE)$ix
            counts <- counts[, reorder.vect]
            fads <- depths[as.numeric(lapply(apply(counts, 
               2, funny), max))]
            reorder.vect <- sort(fads, decreasing = TRUE, 
               index.return = TRUE)$ix
            counts <- counts[, reorder.vect]
            reorder.vect <- sort(as.character(taxa), index.return = TRUE)$ix
            counts <- counts[, reorder.vect]
         }
         else if (pmatch(reorder, "by.count", nomatch = FALSE)) {
            reorder.vect <- sort(colSums(counts), decreasing = TRUE, 
               index.return = TRUE)$ix
            counts <- counts[, reorder.vect]
         }
      }
      else if (length(reorder) == ncol(counts)) {
         reorder.vect <- as.numeric(reorder)
         counts <- counts[, reorder.vect]
      }
      else {
         stop("argument to reorder not understood")
      }
   }
   if (!is.null(tax.cat)) {
      if (length(tax.cat) != ncol(counts)) {
         warning("taxon category labels seem to be the wrong length")
         tax.cat <- NULL
      }
   }
   if (is.null(tax.cat)) {
      tax.cat <- rep("", ncol(counts))
   }
   
   if (xaxis.labels == 'names') { 
      xaxis.labels <- colnames(counts)
   }
   if (xaxis.labels == 'numbers') {
      colnum <- dim(counts)[2]
      xaxis.labels <- as.character(c(1:colnum))
   }
   
   if(xaxis.labels == 'alphanum') {
      xaxis.labels <- colnames(counts)
   }
   
   ad <- a.datums(strat.column(counts = counts, depths = depths))
   plot(1:ncol(counts), ylim = c(max(depths), min(depths)), 
      type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "Age (Ma)")
   segments(1:ncol(counts), ad[, 1], 1:ncol(counts), ad[, 2], 
      lwd = llwd, col = llcol, lty = llty)
   if (baselines == TRUE) { segments(1:ncol(counts), ad[, 2], 1:ncol(counts), rep(par()$usr[3], 
      ncol(counts)), col = 'lightblue', lty = 3, lwd = 0.5) } #col=grey(0.5)
   axis(1, at = 1:ncol(counts), labels = xaxis.labels, cex.axis = cex.xaxis, ##### labels = colnames(counts)
      las = 3)
   axis(2, las = 1, cex.axis = cex.yaxis)
   axis(2, at = depths, labels = FALSE, tck = -0.01)
   
   #added code: 
   if (is.null(colors.vec)) colors<-rainbow(max(counts), end=5/6) else colors<-colors.vec
   sizes<-seq(0,(max(counts)-1), 1)
   sizes<-((sizes/(max(sizes)/large.size)) + large.size)
   #groups
   group.fn<-function(x) {
      for (i in 1:length(x)) {
         if(x[i]>=2 & x[i]<=5) x[i]<-4
         if(x[i]>=6 & x[i]<=10) x[i]<-7
         if(x[i]>10) x[i] <- 15
      }
      return(x)
   }
   if (count.group == TRUE) counts.plt<-group.fn(counts) else counts.plt<-counts
   
   # make the plots! 
   for (i in 1:ncol(counts.plt)) {
      num.val<-c(counts.plt[,i][counts.plt[,i]>0])
      if (pch.points == 'by.count') points<-num.val else points<-pch.points
      if (col.points == 'by.count') cols<-c(colors[num.val]) else cols<-col.points
      if (cex.points == 'by.count') points.cex<-c(sizes[num.val]) else points.cex<-cex.points
      plocs <- depths[(counts.plt > 0)[, i]]
      points(rep(i, length(plocs)), plocs, cex=points.cex, pch=points, col=cols,
         ...)
   }
   
   # add a legend
   if(is.null(legend.values)) {legend.values <- seq(1:max(counts.plt))}
   if(is.null(legend.loc)) {legend.loc='topleft'} 
   if (legend == TRUE & count.group == FALSE) { 
      if (pch.points == 'by.count') leg.pch<-c(1:max(counts.plt)) else leg.pch<-pch.points
      if (col.points == 'by.count') leg.col<-colors else leg.col<-col.points
      if (cex.points == 'by.count') leg.cex<-sizes else leg.cex<-cex.points
      legend(legend.loc, legend=legend.values, pch=leg.pch, 
         col=leg.col, pt.cex=leg.cex, cex=large.size, title=legend.title,
         horiz = legend.horiz, bg = legend.bg) }
   if (legend == TRUE & count.group == TRUE) {
      if (pch.points == 'by.count') leg.pch<-sort(unique(as.vector(counts.plt)))[2:5] else leg.pch<-pch.points
      if (col.points == 'by.count') leg.col<-colors[sort(unique(as.vector(counts.plt)))[2:5]] else leg.col<-col.points
      if (cex.points == 'by.count') leg.cex<-sizes[sort(unique(as.vector(counts.plt)))[2:5]] else leg.cex<-cex.points
      legend(legend.loc, legend=c('1', '2-5', '6-10', '10+'), pch=leg.pch, 
         col=leg.col, pt.cex=leg.cex, cex=large.size, title=legend.title,
         horiz = legend.horiz, bg = legend.bg) }
   
   if(return.xaxis == TRUE) {
      return(colnames(counts))
   }
}

##### 9. make.inp(x, filename, header=NULL) #####
make.inp<-function(x, filename, header=NULL) {  #make input file for RMark, x is table sc$counts 
   mat<-as.vector(x)
   mat<-matrix(mat, nrow=length(x[1,]), ncol=length(x[,1]), byrow=T)
   rownames(mat)<-colnames(x)
   mat<-ifelse(mat>0, 1, 0)  #change to 1's and 0's
   if(is.null(header)) header=filename
   filefoo<-file(paste(filename, sep=''), 'w')
   writeLines(paste('/*', header, '*/', sep=''), filefoo)
   for(i in 1:length(mat[,1])) {
      mm<-mat[i,]
      writeLines(paste('/* ', rownames(mat)[i], ' */ ', paste(mm, sep='', collapse=''), ' 1;', sep=''), filefoo) }
   close(filefoo)
}

##### 10. new.nmds.plotting.function(time.splits) #####
### new plotting data frame function ###

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

# nmds.plot.list.new <- new.nmds.plotting.dfs(time.splits)
# 
# df <- nmds.plot.list.new$NMDS3.con
# 
# sum(apply(df[, 16:19],1, sum) > 1) #check eocene, seems to be good (should be = 0)
# sum(apply(df[, 12:15],1, sum) >1 ) #check P2, seems to be good (should be = 0)
# sum(apply(df[, 8:11],1, sum) > 1) #check P1
# sum(apply(df[, 4:7],1, sum) > 1) #check cretaceous
# sum(apply(df[, 4:19], 1, sum) == 0) #check if there are any which do not appear at all (0 values) - this should also be 0 at the end
# 
# # which(apply(df[, 4:19], 1, sum) ==0) #find the ones that have no spot in this. 
# # # the function misses one where fad and lad are time 5 (aka originated after last point in the Eocene. 
# #in a perfect world, I think we're actually set here..... 
