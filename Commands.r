### Load Pareto/NBD model code
folder<- "C:/Users/Brian Do/Desktop/Eventbrite/"
source(paste(folder,"paretonbd.r",sep=""))

### Load Data
eventbrite <- read.csv(paste(folder,"eventbrite_Jan09_10weeks.csv", sep=""),header=TRUE)
data.eb<-data.frame(p1x=eventbrite$p1x,tx<-eventbrite$t_x,T=eventbrite$T,p2x=eventbrite$p2x)
actual <- read.csv(paste(folder,"eventbrite_Jan09_10weeksCumSls.csv", sep=""), header=TRUE)
actual <-actual$Cumsls
wks<-max(data.eb$T)

### Remove outliers
data.eb<-data.eb[which(data.eb[,1]<9),]


### Commands ###

### Calibrate model and get parameter estimates
results.eb<-PNBD.est(data.eb,c(1,1,1,1))

### Compute P(alive)
compute_pactive(data.eb,results.eb$params)

### Compute Conditional Expectation
compute_ce(59,data.eb,results.eb$params)

### Create data for k-means clustering
x<-compute_ce(59,data.eb,results.eb$params)$ce

### Create clusters
pc<-plot_clusters2(data.eb,results.eb$params,x,7)

### 80/20 computation
tresh=1
sum(pc$trans[which(pc$ce>=tresh)])/sum(pc$trans)
sum(pc$predtrans[which(pc$ce>=tresh)])/sum(pc$predtrans)
sum(pc$size[which(pc$ce>=tresh)])/sum(pc$size)

### Create tracking plots
create_tracking_plots(data.eb, results.eb$params)
