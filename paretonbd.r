### h2f1 Gaussian Hypergeometric Function
h2f1 <- function(a,b,c,z){
	lenz <- length(z)
	j = 0
	uj <- 1:lenz
	uj <- uj/uj
	y <- uj
	lteps <- 0

	while (lteps<lenz){
   		lasty <- y
   		j <- j+1
   		uj <- uj*(a+j-1)*(b+j-1)/(c+j-1)*z/j
   		y <- y + uj
   		lteps <- sum(y==lasty)
	}
	y
}


### Pareto NBD Function
pareto_nbd_ll <- function(paramslog,data){
	param <- exp(paramslog)
	r = param[1]
	alpha = param[2]
	s = param[3]
	beta = param[4]

	p1x<-data$p1x
	T<-data$T
	tx<-data$tx

	maxab <- max(alpha, beta)
	absab <- abs(alpha - beta)
	param2 <- s+1
	if (alpha<beta) {
		param2 <- r + p1x
	}
	Lpart1 <- lgamma(r+p1x)+r*log(alpha)+s*log(beta)-lgamma(r)
	part2 <- 1/((alpha+T)^(r+p1x)*(beta+T)^s)
	if (absab==0){
		F1 <- 1/((maxab+tx)^(r+s+p1x))
		F2 <- 1/((maxab+T)^(r+s+p1x))
	} else {
		F1=h2f1(r+s+p1x,param2,r+s+p1x+1,absab/(maxab+tx))/((maxab+tx)^(r+s+p1x))
   		F2=h2f1(r+s+p1x,param2,r+s+p1x+1,absab/(maxab+T))/((maxab+T)^(r+s+p1x))
	}
	f <- -sum(Lpart1+log(part2+(s/(r+s+p1x))*(F1-F2)))
	f
}


### Computer P(active)
compute_pactive<-function(data,params){
	r<-params[1]
	alpha<-params[2]
	s<-params[3]
	beta<-params[4]

	p1x<-data$p1x
	p2x<-data$p2x
	T<-data$T
	tx<-data$tx

	maxab = max(alpha,beta)
	absab = abs(alpha-beta)
	param2 = s+1
	if (alpha < beta){
    		param2 = r+p1x
	}

	F0 = (alpha+T)^(r+p1x)*(beta+T)^s;
	F1=h2f1(r+s+p1x,param2,r+s+p1x+1,absab/(maxab+tx))/((maxab+tx)^(r+s+p1x))
	F2=h2f1(r+s+p1x,param2,r+s+p1x+1,absab/(maxab+T))/((maxab+T)^(r+s+p1x))
	pactive = 1/(1+(s/(r+s+p1x))*F0 *(F1-F2))

	pa_actual <- 1:(max(p1x)+1)
	pa_actual <- pa_actual-pa_actual
	pa_est <- 1:(max(p1x)+1)
	pa_est <- pa_est-pa_est
	np1x <- 1:(max(p1x)+1)
	np1x <- np1x-np1x

	for (y in unique(p1x)){
		isx <- which(p1x==y)
		np1x[y+1] <- length(isx)
		pa_actual[y+1]<-sum(p2x[isx]>0)/np1x[y+1]
		pa_est[y+1]<-sum(pactive[isx])/np1x[y+1]
	}

	censor <- 7
	denom <- sum(np1x[(censor+1):length(np1x)])

	pa_act_cen <- pa_actual[1:censor]
	pa_act_cen[censor+1] <- sum((np1x[(censor+1):length(np1x)]*pa_actual[(censor+1):length(np1x)])/denom)

	pa_est_cen <- pa_est[1:censor]
	pa_est_cen[censor+1] <- sum((np1x[(censor+1):length(np1x)]*pa_est[(censor+1):length(np1x)])/denom)

	par(mfrow=c(1,1))
	plot(0:censor,pa_act_cen, type="l", axes=FALSE,xlab=paste("Activity in weeks 1 - ",wks, sep=""),ylab="P(Alive)")
	lines(0:censor, pa_est_cen, type="o", lty=2, pch=22)
	axis(1, at=0:7,lab=c("0","1","2","3","4","5","6","7+"))
	axis(2, las=1, at=.1*0:10)
	title(main="Probability Alive Given Prior Activity")
	legend(0, max(pa_est_cen,pa_act_cen), c("Actual","Pareto/NBD"), cex=.8,pch=-1:22,lty=1:2)

	list(actual=pa_actual,est=pa_est,pactive=pactive)
}


### Computer ce
compute_ce<-function(t, data,params){
	p1x<-data$p1x
	p2x<-data$p2x
	T<-data$T
	tx<-data$tx
	pactive <- compute_pactive(data,params)$pactive
	r<-params[1]
	alpha<-params[2]
	s<-params[3]
	beta<-params[4]

	tmp1 <- (r+p1x)*(beta+T)/((alpha+T)*(s-1))
	tmp2 <- ((beta+T)/(beta+T+t))^(s-1)
	ce <- tmp1*(1-tmp2)*pactive

	ce_act <- rep(0,max(p1x)+1)
	ce_est <- rep(0,max(p1x)+1)
	np1x <- rep(0,max(p1x)+1)
	
	for (y in unique(p1x)){
		isx <- which(p1x==y)
		np1x[y+1] <- length(isx)
		ce_act[y+1] <- sum(p2x[isx])/np1x[y+1]
		ce_est[y+1] <- sum(ce[isx])/np1x[y+1]
	}
	
	censor <- 7
	denom <- sum(np1x[(censor+1):length(np1x)])

	ce_act_cen <- ce_act[1:censor]
	ce_act_cen[censor+1] <- sum((np1x[(censor+1):length(np1x)]*ce_act[(censor+1):length(np1x)])/denom)

	ce_est_cen <- ce_est[1:censor]
	ce_est_cen[censor+1] <- sum((np1x[(censor+1):length(np1x)]*ce_est[(censor+1):length(np1x)])/denom)
	par(mfrow=c(1,1))
	plot(0:censor,ce_act_cen, type="l", axes=FALSE, xlab=paste("Activity in weeks 1 - ", wks, sep=""), ylab=paste("Activity in weeks", wks+1, "-", length(actual), sep=" "))
	lines(0:censor,ce_est_cen, lty=2, type="o", pch=22)
	axis(1, at=0:7,lab=c("0","1","2","3","4","5","6","7"))
	axis(2, las=1, at=5*0:max(ce_est_cen,ce_act_cen))
	title(main="Expected Activity Given Prior Activity")
	legend(0, max(ce_est_cen,ce_act_cen), c("Actual","Pareto/NBD"), cex=.8,pch=-1:22,lty=1:2)
	list(ce=ce, actual=ce_act, est=ce_est, ce_act_cen=ce_act_cen, ce_est_cen=ce_est_cen)
}

### Create Tracking Plots
create_tracking_plots <- function(data,params){
	r<-params[1]
	alpha<-params[2]
	s<-params[3]
	beta<-params[4]
	p1x<-data$p1x
	p2x<-data$p2x
	T<-data$T
	tx<-data$tx
	ns<-rep(0,31)
	num <- round(max(T)*7+1)
	for (i in 1:31){
		ns[i]<-sum(round(T,2)==round((num-i)/7,2))
	}
	
	endwk<-length(actual)
	#endwk<-117
	endday<-endwk*7

	tmp1<-r*beta/(alpha*(s-1))
	tmpcumsls1<-rep(0,endday)
	for (i in 1:endday){
		tmp2 <- (beta/(beta+i/7))^(s-1)
		tmpcumsls1[i] <- tmp1*(1-tmp2)
	}

	tmpcumsls2 <- matrix(0, 31, endday)
	for (i in 1:31){
		tmpcumsls2[i,]<-c(rep(0,i),tmpcumsls1[1:(endday-i)])
	}
	
	cumrptsls <- rep(0,endwk)
	dailysls1 <- ns*tmpcumsls2
	dailysls<-rep(0,endday)
	for (i in 1:endday){
		dailysls[i]<-sum(dailysls1[,i])
	}
	for (i in 1:endwk){
		cumrptsls[i] <- dailysls[i*7]
	}
	
	incrptsls <- c(cumrptsls[1], diff(cumrptsls))
	incactual <- c(actual[1], diff(actual))
	
	par(mfrow=c(2,1))
	
	plot(1:endwk,actual[1:endwk], type="l", xlab="Week", ylab="Cum repeat activity")
	lines(cumrptsls, lty=2)
	abline(v=wks, lty=2)
	title(main="Cumulative Repeat Activity")
	legend(endwk*.75, max(actual[1:endwk],cumrptsls)*.4, c("Actual","Pareto/NBD"), cex=.8,lty=1:2)	

	plot(1:endwk,incactual[1:endwk], type="l", xlab="Week", ylab="Weekly repeat activity",ylim=c(0,max(incactual,incrptsls)))
	lines(incrptsls, lty=2)
	abline(v=wks, lty=2)
	title(main="Weekly Repeat Activity")
	legend(endwk*.75, max(incactual,incrptsls)*.95, c("Actual","Pareto/NBD"), cex=.8,lty=1:2)			

	list(actual=actual, cumsls=cumrptsls)
}

### Estimate Parameters

PNBD.est<-function(data,params.start){
	p1x<-data$p1x
	T<-data$T
	tx<-data$tx

	inits<-log(params.start)
	results<-optim(inits,pareto_nbd_ll, data=data)
	params <- exp(results$par)
	ll<-results$value
	cat("\n Estimates",
		"\n r =", round(params[1],3),
		"\n alpha =",round(params[2],3),
		"\n s =",round(params[3],3),
		"\n beta =",round(params[4],3),
		"\n LL =",round(ll,3),"\n\n"
		)
	list(params=params,ll=ll)
}





wss_plot<-function(x){	
	par(mfrow=c(1,1))
	wss<-rep(0,10)
	n<-length(x)
	#n<-nrow(x)
	#wss[1]<-(n-1)*var(x)
	wss[1]<-(n-1)*sum(apply(x,2,var))
	for (i in 2:10){
		wss[i]<-sum(kmeans(x,i,9999,100)$withinss)
	}
	plot(1:10,wss, type="b")
}

plot_clusters2<-function(data,params,x,k){
	cl<-kmeans(x,k, 99999, 100)
	p2x<-data$p2x
	pa<-compute_pactive(data,params)$pactive
	p2x_act<-rep(0,k)
	pa_act<-rep(0,k)
	pa_est<-rep(0,k)
	mean<-matrix(0,k,4)
	for (i in 1:k){
		p2x_act[i]<-mean(p2x[which(cl$cluster==i)])
		pa_act[i]<-sum(p2x[which(cl$cluster==i)]>0)/length(p2x[which(cl$cluster==i)])
		pa_est[i]<-mean(pa[which(cl$cluster==i)])
	}
	b<-matrix(0,k,5)
	b<-cbind(cl$size,pa_est,cl$centers,p2x_act,pa_act,1:k)
	b<-b[order(b[,1]),]
	predtrans=round(b[,1]*b[,3])
	trans=b[,1]*b[,4]
	ptpre=matrix(0,k,k)
	tpre=matrix(0,k,k)
	cumclu=matrix(0,k,k)
	for (i in 1:k){
		ptpre[1:i,i]=predtrans[1:i]
		tpre[1:i,i]=trans[1:i]
		cumclu[1:i,i]=b[1:i,1]
		mean[i,]<-apply(data[which(cl$cluster==b[i,6]),],2,mean)
	}
	cumpredtrans=apply(ptpre,2,sum)
	cumtrans=apply(tpre,2,sum)
	cumclu=apply(cumclu,2,sum)
	par(mfrow=c(2,2))
	plot(cbind(jitter(pa,1,0),jitter(x,1,0)), col=cl$cluster, xlab="P(Alive)", ylab="Conditional Expectation", main="Clusters")
	points(cbind(pa_est,cl$centers),col=sort(unique(cl$cluster)),cex=2.5, pch=13)
	legend(x="topleft",legend=b[,1],bty="n", fill=b[,6])
	plot(b[,4], type="l", xlab="Cluster (smallest to largest)", ylab="Conditional Expectation", ylim=c(0,max(b[,3:4])), main="Conditional Expectation By Cluster")
	lines(b[,3], lty=2, type="o", pch=22)
	legend(x="topright",legend=c("Actual","Pareto/NBD"),lty=1:2,pch=-1:22)
	#plot(b[,5], type="l",ylim=0:1, xlab="Cluster (smallest to largest)", ylab="P(Alive)", main="P(Alive) Prediction By Cluster")
	#lines(b[,2], lty=2, type="o", pch=22)
	#legend(x="bottomleft",legend=c("Actual","Pareto/NBD"),lty=1:2,pch=-1:22)
	plot(trans, type="l", xlab="Cluster (smallest to largest)", ylab="Future Activity", ylim=c(0,1.5*max(trans,predtrans)), main="Future Activity By Cluster")
	lines(predtrans, lty=2, type="o", pch=22)
	legend(x="topleft",legend=c("Actual","Pareto/NBD"),lty=1:2,pch=-1:22)
	plot(cumtrans, type="l", xlab="Cluster (smallest to largest)", ylab="Cumulative Future Activity", ylim=c(0,max(cumtrans,cumpredtrans,cumclu)), main="Cumulative Future Activity By Cluster")
	lines(cumpredtrans, lty=2, type="o", pch=22)
	lines(cbind(1:k,cumclu), type="b")
	legend(x="topleft",legend=c("Actual","Pareto/NBD","Cum Cluster Size"),lty=c(1,2,1),pch=-1:22)
	list(size=b[,1],ce=b[,3],p2x=b[,4],predtrans=predtrans,trans=trans,pa_stat=sum((b[,5]-b[,2])^2/b[,2]), ce_stat=sum((b[,4]-b[,3])^2/b[,3]), mean=mean)
}
