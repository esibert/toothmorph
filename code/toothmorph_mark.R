########################################################################
#
#     Toothmorph Evolutionary Rates: CMR Analyses
#
#
#
######################################################################


##### Load/save datasets #####
source('code/toothmorph/toothmorph_functions.R')  #

# load('RData/toothmorph_mark.RData')
# save.image('../RData/toothmorph_mark.RData')

# accum.596<-read.csv('csv/596_IAR_IMAGES.csv', header=T)
# accum.596<-accum.596[,-1]
# ma3 <- c(1,1,1)/3

library(RMark)
library(Hmisc)
# MarkPath='C:/Software/Apps/MARK'   #Tell it where MARK is located!
setwd('mark/') #need to put the mark files in a separate file, otherwise everything is a mess

## Useful functions for handling RMark input
#generates .inp file from counts table
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

#calls and extracts parameter estimates from a model-list
par.estimates<-function(model.list, rank) {         
   #model.list is the marklist to pull the model from
   #rank is model rank within model table
   #pars is vector of parameters, eg c('Phi', 'p', 'f')
   
   model<-model.list[[as.numeric(row.names(model.list$model.table)[rank])]]
   pars<-names(model$parameters)
   results<-model$results$real
   par.names<-rownames(results)
   par.names<-word(par.names,1,1) #extract first word of rownames
   
   for ( i in 1:length(pars)) {
      p <- pars[i]                                 #parameter [i] name
      par.index<-which(p == par.names)             #indices of parameter [i] 
      par.output<-results[par.index,]              #call the output for parameter [i]
      assign(paste(p, '.df', sep=''), par.output)  #create object of par[i].df
   }
   
   pars.list<-c(paste(pars, '.df', sep=''))     # names for parameter dataframes
   pars.sep<-mget(pars.list)                    # make list from the objects in pars.list
   names(pars.sep)<-pars
   return(pars.sep)
}

# indexing call: 
morphdat.orig[grep(tax, morphdat.orig$Morphotype_name),]


#### PICK A DATASET ####
sc<-sc.con



# step 1: set up MARK, and relevant variables to run
# a) define AgeID.rev (dt.rev from foote)
#AgeID.rev<-AgeID.unique[(length(AgeID.unique)-1):1] #for plotting, mostly

effort<-c(rep(0,13), rep(1,8)) #allows for different sampling effort in Eocene vs. Pal/Cret
samplesize<-rev(as.vector(apply(sc$counts, 1, sum))) #number of teeth considered in each time bin

# b) make *.inp file for the dataset we are working with then call it in

range.counts<-sc$counts #order by default is Eocene to Cretaceous
range.counts.rev<-range.counts[length(range.counts[,1]):1,] #reverse it so that order is Cretaceous to Eocene
make.inp(x=range.counts.rev, filename='mark/rangesrev_con.inp', header='ages reversed full dataset, 21 age points')
rangesrev<-convert.inp('mark/rangesrev_con')
rm(range.counts, range.counts.rev)

# Step 2: Run MARK on different models: Pradel-recruitment, Pradel-lambda, POPAN
# process data into "processed" (dp) and design (ddl); function and run MARK on all possible models. 


##### a) Pradel-recruitment #####
dp.pradrec<-process.data(rangesrev, model='Pradrec', time.intervals=dt.rev)  
ddl.pradrec<-make.design.data(dp.pradrec)
ddl.pradrec$p$effort<-effort  #add "effort" variable to test for variable p
ddl.pradrec$p$samplesize<-samplesize #add sample size variable to test for variable p
# ddl.pradrec$Phi$samplesize <- samplesize[2:21]
# ddl.pradrec$f$samplesize <- samplesize[1:20]

Pradel.recruit<-function() {
   setwd('C:/Elizabeth_Files/Personal/Research/R_Files/Manuscripts/tooth_morph/mark/') #need to put the mark files in here... 
   #Create formulas for Phi, p, and f
   Phi.time<-list(formula=~time)
   Phi.dot<-list(formula=~1)
   p.time<-list(formula=~time)
   p.dot<-list(formula=~1)
   p.effort<-list(formula=~effort)
   #p.samplesize <- list(formula =~ samplesize)
   #p.effortplussamplesize<-list(formula =~ effort + samplesize)
   f.time<-list(formula=~time)
   f.dot<-list(formula=~1)
   #create models to run by combining the above formulas, looking for objects with Phi. p. and f. at beginning
   cml<-create.model.list('Pradrec')
   #run all the models
   results<-mark.wrapper(cml,data=dp.pradrec,ddl=ddl.pradrec,output=FALSE,silent=TRUE)
   return(results)
}

#save output from LSS and sample size giant models... 

## currently officially saved in the right locations
# Pradrec.con.p.samplesize.effort<-Pradrec.con 


#Original dataset
# Pradrec.original<-Pradel.recruit()   #run models
# Pradrec.original  #displays table

## Run Mark properly: Conservatively trimmed dataset
Pradrec.con<-Pradel.recruit()    #run models
Pradrec.con  #displays table

# save output:
Pradrec.con.p.allcombos<-Pradrec.con
save.image('../RData/toothmorph_mark.RData')

#Calculations useful for plotting
Pradrec.con.avg<-model.average(Pradrec.con, vcv=TRUE) # vcv=TRUE includes 
Pradrec.con.avg<-Pradrec.con.avg$estimates
### Corrected SE values for extinction just for pr.con ##
pr.con.ext<-Pradrec.con.avg[1:20,]
pr.con.ext$estimate<-1-pr.con.ext$estimate
pr.con.ext$ucl<-1-pr.con.ext$ucl
pr.con.ext$lcl<-1-pr.con.ext$lcl
pr.con.orig<-Pradrec.con.avg[42:61,]

pr.con.orig[pr.con.orig == Inf] <- 0

# 
# # #Liberally trimmed dataset
# Pradrec.lib<-Pradel.recruit()   #run models
# Pradrec.lib  #display table 


##### b) Pradel-Lambda model #####
dp.pradlam<-process.data(rangesrev, model='Pradlambda', time.intervals=dt.rev)  
ddl.pradlam<-make.design.data(dp.pradlam)
ddl.pradlam$p$effort<-effort

Pradel.lambda<-function() {
   #Create formulas for Phi, p, and f
   Phi.time<-list(formula=~time)
   Phi.dot<-list(formula=~1)
   p.time<-list(formula=~time)
   p.dot<-list(formula=~1)
   p.effort<-list(formula=~effort)
   Lambda.time<-list(formula=~time)
   Lambda.dot<-list(formula=~1)
   #create models to run by combining the above formulas, looking for objects with Phi. p. and f. at beginning
   cml<-create.model.list('Pradlambda')
   #run all the models
   results<-mark.wrapper(cml,data=dp.pradlam,ddl=ddl.pradlam,output=FALSE,silent=TRUE)
   return(results)
}

# #Original dataset
# Pradlam.original<-Pradel.lambda()   #run models
# Pradlam.original  #displays table

# #Conservatively trimmed dataset 
# Pradlam.con<-Pradel.lambda()    #run models
# Pradlam.con  #displays table
# 
# #Liberally trimmed dataset
Pradlam.lib<-Pradel.lambda()   #run models
Pradlam.lib  #display table 


##### c) POPAN Model #####
dp.popan<-process.data(rangesrev, model='POPAN', time.intervals=dt.rev)  
ddl.popan<-make.design.data(dp.popan)
ddl.popan$p$effort<-effort

popan.fn<-function() {
   #Create formulas for Phi, p, and f
   Phi.time<-list(formula=~time)
   Phi.dot<-list(formula=~1)
   p.time<-list(formula=~time)
   p.dot<-list(formula=~1)
   p.effort<-list(formula=~effort)
   pent.time<-list(formula=~time)
   pent.dot<-list(formula=~1)
   N.dot<-list(formula=~1)
   #create models to run by combining the above formulas, looking for objects with Phi. p. and f. at beginning
   cml<-create.model.list('POPAN')
   #run all the models
   results<-mark.wrapper(cml,data=dp.popan,ddl=ddl.popan,output=FALSE,silent=TRUE)
   return(results)
}
# #Original dataset
# popan.original<-popan.fn()   #run models
# popan.original  #displays table

# #Conservatively trimmed dataset 
# popan.con<-popan.fn()    #run models
# popan.con  #displays table
# 
# #Liberally trimmed dataset
popan.lib<-popan.fn()   #run models
popan.lib  #display table 


##### Step 3: Calculate Model Averages #####
Pradrec.lib.avg<-model.average(Pradrec.lib)
Pradrec.con.avg<-model.average(Pradrec.con)
Pradrec.orig.avg<-model.average(Pradrec.original)

Pradlam.lib.avg<-model.average(Pradlam.lib)
Pradlam.con.avg<-model.average(Pradlam.con)
Pradlam.orig.avg<-model.average(Pradlam.original)

Popan.lib.avg<-model.average(popan.lib)
Popan.con.avg<-model.average(popan.con)
Popan.orig.avg<-model.average(popan.original)

##### Step 4: Extract parameters for each set of models: #####

#Pradel Recruitment: 61 parameters: 
# Phi: 1-20; 
# p: 21-41; (includes all AgeID's
# f: 42-61;
pr.orig.ext<-Pradrec.orig.avg[1:20,]
pr.orig.ext[,2] <- 1-pr.orig.exit[,2]  #1-survival is extinction
pr.orig.orig<-Pradrec.orig.avg[42:61,]

### Corrected SE values for extinction just for pr.con ##
pr.con.ext<-Pradrec.con.avg[1:20,]
pr.con.ext[,2]<-1-pr.con.ext[,2]
pr.con.orig<-Pradrec.con.avg[42:61,]

pr.lib.ext<-1-Pradrec.lib.avg[1:20,]
pr.lib.orig<-Pradrec.lib.avg[42:61,]

# Pradel Lambda: 61 parameters
# Phi: 1-20; 
# p: 21-41; (includes all AgeID's
# Lambda: 42-61;
pl.orig.ext<-1-Pradlam.orig.avg[1:20,]
pl.orig.orig<-Pradlam.orig.avg[42:61,]-Pradlam.orig.avg[1:20,] #lambda - phi
pl.con.ext<-1-Pradlam.con.avg[1:20,]
pl.con.orig<-Pradlam.con.avg[42:61,]-Pradlam.con.avg[1:20,]
pl.lib.ext<-1-Pradlam.lib.avg[1:20,]
pl.lib.orig<-Pradlam.lib.avg[42:61,]-Pradlam.lib.avg[1:20,]

# POPAN: 62 parameters
# Phi: 1-20; 
# p: 21-41; (includes all AgeID's
# pent: 42-61;
# N #metapopulation size
pop.orig.ext<-1-Popan.orig.avg[1:20,]
pop.orig.orig<-Popan.orig.avg[42:61,]
pop.con.ext<-1-Popan.con.avg[1:20,]
pop.con.orig<-Popan.con.avg[42:61,]
pop.lib.ext<-1-Popan.lib.avg[1:20,]
pop.lib.orig<-Popan.lib.avg[42:61,]


##### Plot this #####

#this function isn't working - ack! 
 error.bars <- function(x, y, upper, lower=upper, length.arr=0.1, ...){
   if(length(x) != length(y) | length(y) != length(lower) | length(lower) != length (upper))
      stop("vectors must be the same length")
   arrows(x, y+upper, x, y-lower, angle=90, code=3, length=length.arr, ...) 
}

AgeID.int <- c() 
for(i in 1:length(AgeID.unique)-1) {
   AgeID.int[i] <- mean(c(AgeID.unique[i], AgeID.unique[i+1]))
}
AgeID.int<-rev(AgeID.int) #go from older to younger, as in capture history calculations 

##### Cleaned Figure #####

# select appropriate model output
Pradrec.con<-Pradrec.con.p.effort
# calculate model averages for plotting
Pradrec.con.avg<-model.average(Pradrec.con, vcv=TRUE)
# extract modeled origination and extinction parameters
pr.con.ext<-Pradrec.con.avg[1:20,]
pr.con.ext[,2]<-1-pr.con.ext[,2]
pr.con.orig<-Pradrec.con.avg[42:61,]

### NOTE FIGURE SAVE LOCATION ASSUMES THE WD IS THE MARK WD. THIS IS NOT ALWAYS THE CASE ##
pdf('../figures/cmr_orig_ext_accum_errbar_peffort.pdf', width=8, height=5, useDingbats = FALSE)
# Set up plot parameters
par(mar=c(5,4,4,4))
xlims <- c(73, 43)
ext.color<-'firebrick'
orig.color<-'blue3'
accum.color<-'gray30'

# ext.color<-'red'
# orig.color<-'green'
# accum.color<-'gray30'


#make blank plot and annotations
plot(AgeID.int, pr.con.orig$estimate, type='n', xlim=xlims, ylim=c(-0.03,0.36), 
   main='CMR Evolutionary Rate Estimates', 
   xlab='', ylab='')
rect(64.3, -0.05, 63, 0.39, col='gray80', border=NA)
rect(61.8, -0.05, 55.8, 0.39, col='gray90', border=NA)
rect(59.2, -0.05, 58, 0.39, col='gray79', border=NA)
mtext(text='Per-morphotype rate estimate', side=2, line=2.5) #x-axis label
mtext(text='Age (Ma)', side=1, line=2.5)
text(63.5, 0.365, 'Pulse 1')
text(58.5, 0.365, 'Pulse 2')
abline(v=66.5, col='black', lwd=1)
abline(h=0, col='gray60', lty=2, lwd=0.8)
box()

## add extinction ##
errbar(AgeID.int, pr.con.ext$estimate, type='o', add=TRUE, 
   yplus=pr.con.ext$estimate+pr.con.ext$se, yminus=pr.con.ext$estimate-pr.con.ext$se, 
   errbar.col=ext.color, col=ext.color, pch=16, lwd=1)

#with ucl/lcl
errbar(AgeID.int, pr.con.ext$estimate, type='o', add=TRUE,
   yplus=pr.con.ext$lcl, yminus=pr.con.ext$ucl, 
   errbar.col=ext.color, col=ext.color, pch=17, lwd=1)

#with se/lcl
errbar(AgeID.int, pr.con.ext$estimate, type='o', add=TRUE,
   yplus=pr.con.ext$estimate+pr.con.ext$se, yminus=pr.con.ext$ucl, 
   errbar.col=ext.color, col=ext.color, pch=17, lwd=1)

## add origination ##
errbar(AgeID.int, pr.con.orig$estimate, type='o', add=TRUE,
   yplus=pr.con.orig$estimate+pr.con.orig$se, yminus=pr.con.orig$estimate-pr.con.orig$se, 
   errbar.col=orig.color, col=orig.color, pch=17, lwd=1)

#with ucl/lcl
errbar(AgeID.int, pr.con.orig$estimate, type='o', add=TRUE,
   yplus=pr.con.orig$ucl, yminus=pr.con.orig$lcl, 
   errbar.col=orig.color, col=orig.color, pch=17, lwd=1)

#with se/lcl
errbar(AgeID.int, pr.con.orig$estimate, type='o', add=TRUE,
   yplus=pr.con.orig$estimate+pr.con.orig$se, yminus=pr.con.orig$lcl, 
   errbar.col=orig.color, col=orig.color, pch=17, lwd=1)




#legend
legend('topright', legend=c('Origination', 'Extinction', 'Accumulation'), 
   col=c(orig.color, ext.color, accum.color), 
   pch=c(16, 16, -1), lty=c(1,1,2), lwd=c(1,1,1.2), cex=0.9)

# add accumulation curve
par(new=T)
# Accum plot
plot(accum.596$age, filter(accum.596$teeth.mar, ma3), xlim=xlims, ylim=c((-250/12), 250), #y-axis scaled to match 0-value from above
   lty=2, type='l', col=accum.color, lwd=1.5, 
   axes=F, xlab='', ylab='')
axis(4)
mtext(text=parse(text='Teeth %.% cm^-2 %.% myr^-1'), side=4, line=2.5)
dev.off()


##### Plot code used for manuscript figure #####
### P.effort only! ###

# select appropriate model output
Pradrec.con<-Pradrec.con.p.effort
# calculate model averages for plotting and pull parameter estimates
Pradrec.con.avg<-model.average(Pradrec.con, vcv=TRUE)
Pradrec.con.avg<-Pradrec.con.avg$estimates
### Corrected SE values for extinction just for pr.con ##
pr.con.ext<-Pradrec.con.avg[1:20,]
pr.con.ext$estimate<-1-pr.con.ext$estimate
pr.con.ext$ucl<-1-pr.con.ext$ucl
pr.con.ext$lcl<-1-pr.con.ext$lcl
pr.con.orig<-Pradrec.con.avg[42:61,]

pr.con.orig[pr.con.orig == Inf] <- 0

# 

### NOTE FIGURE SAVE LOCATION ASSUMES THE WD IS THE MARK WD. THIS IS NOT ALWAYS THE CASE ##
pdf('../figures/cmr_orig_ext_accum_errbar_peffort.pdf', width=8, height=5, useDingbats = FALSE)
# Set up plot parameters
par(mar=c(5,4,4,4))
xlims <- c(73, 43)
ext.color<-'firebrick'
orig.color<-'blue3'
accum.color<-'gray30'

# ext.color<-'red'
# orig.color<-'green'
# accum.color<-'gray30'


#make blank plot and annotations
plot(AgeID.int, pr.con.orig$estimate, type='n', xlim=xlims, ylim=c(-0.03,0.36), 
   main='CMR Evolutionary Rate Estimates', 
   xlab='', ylab='')
rect(64.3, -0.05, 63, 0.39, col='gray80', border=NA)
rect(60.8, -0.05, 55.8, 0.39, col='gray90', border=NA)
rect(59.2, -0.05, 58, 0.39, col='gray79', border=NA)
mtext(text='Per-morphotype rate estimate', side=2, line=2.5) #x-axis label
mtext(text='Age (Ma)', side=1, line=2.5)
text(63.5, 0.365, 'Pulse 1')
text(58.5, 0.365, 'Pulse 2')
abline(v=66.5, col='black', lwd=1)
abline(h=0, col='gray60', lty=2, lwd=0.8)
box()

# # add extinction
# errbar(AgeID.int, pr.con.ext$estimate, type='o', add=TRUE, 
#    yplus=pr.con.ext$estimate+pr.con.ext$se, yminus=pr.con.ext$estimate-pr.con.ext$se, 
#    errbar.col=ext.color, col=ext.color, pch=16, lwd=1)
# 
#with ucl/lcl
errbar(AgeID.int, pr.con.ext$estimate, type='o', add=TRUE,
   yplus=pr.con.ext$lcl, yminus=pr.con.ext$ucl,
   errbar.col=ext.color, col=ext.color, pch=17, lwd=1)

# #with se/lcl
# errbar(AgeID.int, pr.con.ext$estimate, type='o', add=TRUE,
#    yplus=pr.con.ext$estimate+pr.con.ext$se, yminus=pr.con.ext$ucl, 
#    errbar.col=ext.color, col=ext.color, pch=16, lwd=1)

# # add origination 
# errbar(AgeID.int, pr.con.orig$estimate, type='o', add=TRUE,
#    yplus=pr.con.orig$estimate+pr.con.orig$se, yminus=pr.con.orig$estimate-pr.con.orig$se, 
#    errbar.col=orig.color, col=orig.color, pch=17, lwd=1)
# 
#with ucl/lcl
errbar(AgeID.int, pr.con.orig$estimate, type='o', add=TRUE,
   yplus=pr.con.orig$ucl, yminus=pr.con.orig$lcl,
   errbar.col=orig.color, col=orig.color, pch=17, lwd=1)

# #with se/lcl
# errbar(AgeID.int, pr.con.orig$estimate, type='o', add=TRUE,
#    yplus=pr.con.orig$estimate+pr.con.orig$se, yminus=pr.con.orig$lcl, 
#    errbar.col=orig.color, col=orig.color, pch=17, lwd=1)


#legend
legend('topright', legend=c('Origination', 'Extinction', 'Accumulation'), 
   col=c(orig.color, ext.color, accum.color), 
   pch=c(17, 16, -1), lty=c(1,1,2), lwd=c(1,1,1.2), cex=0.9)

# add accumulation curve
par(new=T)
# Accum plot
plot(accum.596$age, filter(accum.596$teeth.mar, ma3), xlim=xlims, ylim=c((-250/12), 250), #y-axis scaled to match 0-value from above
   lty=2, type='l', col=accum.color, lwd=1.5, 
   axes=F, xlab='', ylab='')
axis(4)
mtext(text=parse(text='Teeth %.% cm^-2 %.% myr^-1'), side=4, line=2.5)
dev.off()
