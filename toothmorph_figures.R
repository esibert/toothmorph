####################################
#
#  Toothmorph Figures
#  1. Range Chart(s)
#  2. Range Chart (normalized and grouped by % teeth; colors and legend good) 
#  3. Final NMDS Figure
#  4. CMR origination/Extinction figure
#
################################

library(viridis)
library(stratigraph)
source('code/toothmorph/toothmorph_functions.R')
#setwd("")


##### 1. Range Charts - by count (not used in final manuscript) ##### 


# save.image('RData/toothmorph_setup.RData')
# load('RData/toothmorph_setup.RData')

## note that this figure was  modified in illustrator for color scheme, and to include the accumulation rate curve from Sibert et al 2016.

pdf(file='figures/range_chart_conservative.pdf', width=55, height=30, useDingbats = F)
par(mar=c(38, 6, 4, 4))
rangechart3(sc.con, reorder='lad.by.fad', depths=sc.con$absolute.ages, 
   cex.xaxis=1.8, cex.yaxis=3, cex.points='by.count', llwd=3, col.points='by.count',
   pch.points=19, baselines=TRUE, large.size=4, count.group=FALSE) 
abline(h=66, col='red')
abline(h=56, col='blue')
dev.off()


pdf(file='figures/range_chart_liberal.pdf', width=55, height=30, useDingbats = F)
par(mar=c(38, 6, 4, 4))
rangechart3(sc.lib, reorder='lad.by.fad', depths=sc.lib$absolute.ages, 
   cex.xaxis=1.8, cex.yaxis=3, cex.points='by.count', llwd=3, col.points='by.count',
   pch.points=19, baselines=TRUE, large.size=4, count.group=FALSE) 
abline(h=66, col='red')
abline(h=56, col='blue')
dev.off()

pdf(file='figures/range_chart_original.pdf', width=55, height=30, useDingbats = F)
par(mar=c(38, 6, 4, 4))
rangechart3(sc.orig, reorder='lad.by.fad', depths=sc.orig$absolute.ages, 
   cex.xaxis=1.8, cex.yaxis=3, cex.points='by.count', llwd=3, col.points='by.count',
   pch.points=19, baselines=TRUE, large.size=4, count.group=FALSE) 
abline(h=66, col='red')
abline(h=56, col='blue')
dev.off()


##### 2. Range Chart - normalized #####

# save.image('RData/toothmorph_setup.RData')
# load('RData/toothmorph_setup.RData')

## functions for normalizing and then grouping by % ranges
normalize.sc<-function(sc) {  #returns a normalized "counts" matrix. All rows (ages) should add up to 1. 
   norm.row<-function(row) {
      row/sum(row)
   }
   norm.sc<-apply(sc$counts, 1, norm.row)  #normalize the matrix by rows
   norm.sc<-t(norm.sc)                    #for some reason I have to transpose the output back to the normal "counts" form
   norm.sc
}
# note that breaks.fn() uses normalize.sc() in it, so by running breaks.fn(), 
# it will automatically normalize the dataset. If it is already normalized, this 
# is fine - as everything will just be divided by 1, and therefore unchanged in value
breaks.fn <- function(sc, splits) {
   norm.sc<-100*normalize.sc(sc)
   splits<-c(0, splits, 100) #make breaks have full length
   for(i in 1:length(splits)-1) {
      norm.sc[norm.sc > splits[i] & norm.sc <= splits[i+1] ] = i
   }
   sc$counts <- norm.sc 
   return(sc)
}

##for plotting simplicity, take the output of normalize.sc, multiply by 100, and truncate values, so that new "counts" go from 0-18

#choose the one to plot

# create new counts matrix and range chart
sc.plot<-sc.con #select which dataset to use
splits<-c(3,6,9,12)
sc.plot <- breaks.fn (sc.plot, splits)

### Relative Abundance Plot ###

# define legend input
splits.legend<-paste(c(0,splits), '-', c(splits,100), '%', sep='')
splits.legend[1]<-paste('<', splits[1], '%', sep='')
splits.legend[length(splits.legend)] <- paste('>', splits[length(splits)], '%', sep='')

# # cols.vec <- c('gray60', 'gray30', 'goldenrod2', 'seagreen', 'darkblue')
# cols.vec <- viridis(5)
# 
# #names version
# pdf(file='figures/range_chart_conservative_relative.pdf', width=55, height=30, useDingbats = F)
# par(mar=c(38, 6, 4, 4))
# rangechart3(sc.plot, reorder='lad.by.fad', depths=sc.con$absolute.ages, 
#    cex.xaxis=1.8, cex.yaxis=3, cex.points="by.count", llwd=1, 
#    col.points="by.count", colors.vec = cols.vec, xaxis.labels = 'names', #xaxis.labels = c('numbers', 'names')
#    pch.points=16, baselines=TRUE, large.size=4, count.group=FALSE, legend.loc = 'bottomright', 
#    legend.values = splits.legend, legend.horiz = TRUE, return.xaxis = TRUE) 
# abline(h=66, col='red')
# abline(h=56, col='blue')
# dev.off()
# 
# #numbers version
# pdf(file='figures/range_chart_conservative_relative_xaxisNum.pdf', width=55, height=30, useDingbats = F)
# par(mar=c(6, 6, 4, 4))
# rangechart3(sc.plot, reorder='lad.by.fad', depths=sc.con$absolute.ages, 
#    cex.xaxis=1.8, cex.yaxis=3, cex.points="by.count", llwd=1, 
#    col.points="by.count", colors.vec = cols.vec, xaxis.labels = 'numbers', #xaxis.labels = c('numbers', 'names')
#    pch.points=16, baselines=FALSE, large.size=4, count.group=FALSE, legend.loc = 'bottomright', 
#    legend.values = splits.legend, legend.horiz = TRUE, return.xaxis = FALSE) 
# abline(h=66, col='red')
# abline(h=56, col='blue')
# dev.off()


## view in R test version

# cols.vec <- rev(viridis(5))
# cols.vec <- rev(magma(6))[2:6]
# cols.vec <- rev(gray.colors(5, start = 0.2, end = 0.7))
# cols.vec <-(viridis(5))

##### USED THIS VERSION # the one that Matt liked
cols.vec <- rev(viridis(5))
cols.vec[1] <- 'gray70'

pdf(file='figures/range_chart_conservative_relative_xaxisNum.pdf', width=55, height=26, useDingbats = F)
sc.plot<-sc.con #select which dataset to use
splits<-c(3,6,9,12)
sc.plot <- breaks.fn (sc.plot, splits)

# pdf(file='figures/range_chart_liberal_relative_xaxisNum.pdf', width=55, height=26, useDingbats = F)
# sc.plot<-sc.lib #select which dataset to use
# splits<-c(3,6,9,12)
# sc.plot <- breaks.fn (sc.plot, splits)

# pdf(file='figures/range_chart_original_relative_xaxisNum.pdf', width=55, height=26, useDingbats = F)
# sc.plot<-sc.orig #select which dataset to use
# splits<-c(3,6,9,12)
# sc.plot <- breaks.fn (sc.plot, splits)

par(mar=c(6, 6, 4, 4))
rangechart3(sc.plot, reorder='lad.by.fad', depths=sc.con$absolute.ages, 
   cex.xaxis=2.5, cex.yaxis=3, cex.points="by.count", llwd=1, llcol = 'lightgray', llty = 3,
   col.points="by.count", colors.vec = cols.vec, xaxis.labels = 'alphanum', #xaxis.labels = c('numbers', 'names')
   pch.points=16, baselines=FALSE, large.size=4, count.group=FALSE, legend.loc = 'bottomright', 
   legend.values = splits.legend, legend.horiz = TRUE, return.xaxis = FALSE) 
abline(h=66, col='red')
abline(h=56, col='blue')
dev.off()

# to get the rownames...# 

####


##### 3. Final NMDS figure (origination only) #####
load('RData/toothmorph_nmds.RData')
#save.image('RData/toothmorph_nmds.RData')


##### Step 1: define data to plot (nmds3.con)
nmds.plot.list.new <- new.nmds.plotting.dfs(time.splits) #make sure you have the right data
df <- nmds.plot.list.new$NMDS3.con  #pull the appropriate NMDS plot

# colors
col.up <- 'black'  #pch = 24
bg.up <- 'blue'
col.down <- 'gray40' #pch = 25
bg.down <- 'gray80'
col.diamond <- 'gray40'  #pch = 23
bg.diamond <- 'gray80'
col.square <- 'gray40'   #pch = 22
bg.square <- 'gray80' 
prior.color <- 'gray60'
#cex value for points
zoom <- 1
hull.border <- TRUE  #border line? 
line.br <- 3  #line type around convex hull
lwd.br <- 1  #convex hull border width

#### Step 2: Define output location
# jpeg(filename='figures/NMDS_figure_triangles_combined.jpg', 
#    width=16, height=8, units='in', res=150)

pdf(file='figures/NMDS_figure_origination_only.pdf', 
   width=16, height=4, useDingbats = FALSE)   #16x4 makes each plot 4x4


##### Step 3: Make the plots ###

par(mfrow = c(1,4))

#par(mfrow=c(2,2))

## Cretaceous plot ##
plot(df$MDS1, df$MDS2, type='n', main = 'Cretaceous (72-66 Ma)', xlab = 'MDS1', ylab = 'MDS2') #create same blank NMDS plot for every plot

# plot Cretaceous teeth
# points(subset(df, c.diamond == 1, select = c(MDS1, MDS2)), pch=23, col = col.diamond, bg = bg.diamond)
points(subset(df, c.diamond == 1, select = c(MDS1, MDS2)), pch=21, col = col.diamond, bg = bg.diamond, cex = zoom)
# points(subset(df, c.down == 1, select = c(MDS1, MDS2)), pch=25, col = col.down, bg = bg.down, cex = zoom)
points(subset(df, c.down == 1, select = c(MDS1, MDS2)), pch=25, col = 'black', bg = 'red', cex = zoom)
points(subset(df, c.up == 1, select = c(MDS1, MDS2)), pch=21, col = col.square, bg = bg.square, cex = zoom)  #no ups for Cretaceous
points(subset(df, c.square == 1, select = c(MDS1, MDS2)), pch=21, col = col.square, bg = bg.square, cex = zoom)

## P1 plot ##
plot(df$MDS1, df$MDS2, type='n', main = 'Pulse 1 (66-60 Ma)', xlab = 'MDS1', ylab = 'MDS2') #create same blank NMDS plot for every plot

# plot Cretaceous convex hull 
Plot_ConvexHull(xcoord=df[which(df$c.up==1 | df$c.down == 1 | df$c.diamond == 1 | df$c.square == 1),]$MDS1,
   ycoord = df[which(df$c.up==1 | df$c.down == 1 | df$c.diamond == 1 | df$c.square == 1),]$MDS2,
   border.line = hull.border, lcolor=prior.color, lwd=lwd.br, lty = line.br, 
   shade=T, scolor=adjustcolor(prior.color,alpha.f = 0.25))

# Plot P1 teeth
points(subset(df, p1.diamond == 1, select = c(MDS1, MDS2)), pch=21, col = col.diamond, bg = bg.diamond)
points(subset(df, p1.down == 1, select = c(MDS1, MDS2)), pch=21, col = col.down, bg = bg.down, cex = zoom)
points(subset(df, p1.up == 1, select = c(MDS1, MDS2)), pch=24, col = col.up, bg = bg.up, cex = zoom)
points(subset(df, p1.square == 1, select = c(MDS1, MDS2)), pch=21, col = col.square, bg = bg.square, cex = zoom)

## P2 plot ##
plot(df$MDS1, df$MDS2, type='n', main = 'Pulse 2 (60-56 Ma)', xlab = 'MDS1', ylab = 'MDS2') #create same blank NMDS plot for every plot

# Plot P1 convex hull
Plot_ConvexHull(xcoord=df[which(df$p1.up==1 | df$p1.down == 1 | df$p1.diamond == 1 | df$p1.square == 1),]$MDS1,
   ycoord = df[which(df$p1.up==1 | df$p1.down == 1 | df$p1.diamond == 1 | df$p1.square == 1),]$MDS2,
   border.line = hull.border, lcolor=prior.color, lwd=lwd.br, lty = line.br,
   shade=T, scolor=adjustcolor(prior.color,alpha.f = 0.25))

# plot P2 teeth
points(subset(df, p2.diamond == 1, select = c(MDS1, MDS2)), pch=21, col = col.diamond, bg = bg.diamond)
points(subset(df, p2.down == 1, select = c(MDS1, MDS2)), pch=21, col = col.down, bg = bg.down, cex = zoom)
points(subset(df, p2.up == 1, select = c(MDS1, MDS2)), pch=24, col = col.up, bg = bg.up, cex = zoom)
points(subset(df, p2.square == 1, select = c(MDS1, MDS2)), pch=21, col = col.square, bg = bg.square, cex = zoom)

## Eocene plot ##
plot(df$MDS1, df$MDS2, type='n', main = 'Eocene (55-43 Ma)', xlab = 'MDS1', ylab = 'MDS2') #create same blank NMDS plot for every plot

# Plot P2 convex hull
Plot_ConvexHull(xcoord=df[which(df$p2.up==1 | df$p2.down == 1 | df$p2.diamond == 1 | df$p2.square == 1),]$MDS1,
   ycoord = df[which(df$p2.up==1 | df$p2.down == 1 | df$p2.diamond == 1 | df$p2.square == 1),]$MDS2,
   border.line = hull.border, lcolor=prior.color, lwd=lwd.br, lty = line.br,
   shade=T, scolor=adjustcolor(prior.color,alpha.f = 0.25))

# plot eocene teeth
# points(subset(df, e.diamond == 1, select = c(MDS1, MDS2)), pch=23, col = col.diamond, bg = bg.diamond)
points(subset(df, e.diamond == 1, select = c(MDS1, MDS2)), pch=21, col = col.diamond, bg = bg.diamond)
# points(subset(df, e.up == 1, select = c(MDS1, MDS2)), pch=24, col = col.up, bg = bg.up, cex = zoom)
points(subset(df, e.down == 1, select = c(MDS1, MDS2)), pch=21, col = col.down, bg = bg.down, cex = zoom)
points(subset(df, e.up == 1, select = c(MDS1, MDS2)), pch=24, col = col.up, bg = bg.up, cex = zoom)
points(subset(df, e.square == 1, select = c(MDS1, MDS2)), pch=24, col = col.up, bg = bg.up, cex = zoom)

dev.off()

##### 4. Origination/Extinction: CMR #####
pdf('figures/cmr_orig_ext_accum_errbar.pdf', width=8, height=5, useDingbats = FALSE)
# Set up plot parameters
par(mar=c(5,4,4,4))
xlims <- c(73, 43)
ext.color<-'firebrick'
orig.color<-'blue3'
accum.color<-'gray30'

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

# add origination 
errbar(AgeID.int, pr.con.orig$estimate, type='o', add=TRUE,
   yplus=pr.con.orig$estimate+pr.con.orig$se, yminus=pr.con.orig$estimate-pr.con.orig$se, 
   errbar.col=orig.color, col=orig.color, lwd=1)

# add extinction
errbar(AgeID.int, pr.con.ext$estimate, add=TRUE, type='o',
   yplus=pr.con.ext$estimate+pr.con.ext$se, yminus=pr.con.ext$estimate-pr.con.ext$se, 
   errbar.col=ext.color, col=ext.color, lwd=1)

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
