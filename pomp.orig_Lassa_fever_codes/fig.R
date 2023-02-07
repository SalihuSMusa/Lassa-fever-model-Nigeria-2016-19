len <- 30
tempco <- rep(1e-12,7*len)
tempfp <- tempco
x<-read.csv('precip_nigeria.csv')
xout<-(1:1156)/365.25+2016
precip<-approx(x[,1],y=(x[,2]),xout)$y
temp1<-precip[which(xout>2016.34-1/12 & xout<2017.12)]
temp2<-precip[which(xout>2017.34-1/12 & xout<2018.12)]
temp3<-precip[which(xout>2018.34-1/12 & xout<2019.12)]
temp<-cbind(temp1,temp2,c(temp3,NA))
monthlab <- as.Date(c('2016-11-1','2017-1-1','2017-3-1','2017-5-1'))
month <- as.numeric(monthlab-as.Date('2015-12-31'))/366+2016
#print(monthlab)
#print(month)
  require('pomp.orig')
  require('pracma')
  d1 <- as.numeric(as.Date(c('2016-10-31','2017-6-7'))-as.Date("2016-10-31"))/366+2016.834
  dir <- getwd()
  dir <- substr(dir,nchar(dir)-4,nchar(dir))
  z<-c('','','')
type<-'pdf'
  file <- paste('~/fitting_',dir,'.',type,sep="")

flag <- 0
setEPS()
if(type=='eps')
  postscript(file,width=16,height=ifelse(flag==0,5.5,10))
if(type=='pdf')
  pdf(file,width=12,height=ifelse(flag==0,5,14))
  par(fig=c(0,1,0,1),las=1,xaxs='i',yaxs='i')
  plot.new()

##################################################################################################

##################################################################################################


load('mle_lassa002.rda')
infectious<-as.null()
for(i in 1:3)
infectious[[i]]<-apply(mle[[1]]@pred.mean[(i-1)*8+c(3,5,6),],2,sum)


load('allregion.rda')
print(length(all))
x <- read.csv('best_lassa.csv',row=1)
y<-read.csv('model_parameters.csv',r=1)
x$np<-sum(y$type=='par') -10 +x$nm
np<-x$np
k<-read.csv('nigeria_LHF.csv',head=1)
nd<-sum(is.na(k[,c(1,3,5)])==FALSE)
#  bic<- -2*x$loglik + np*log(nd)
bic<- -2*x$loglik+2*np*nd/(nd-np-1)
x$bic<-bic
y <- x
  j <-  which.min(bic[1:10])
print(j)

	for(i in 1){
	y[i,] <- x[j[i],]
        all[[i]] <- all[[j[i]]]
	}
  dyn.load('lassamodel.so')
  source('lassamodel.R')
  mle <- all[[1]]
  params <- as.numeric(y[1,])
  names(params) <- colnames(y)
  mle@coef <- par.trans(params[match(names(coef(mle)),names(params))])
  n <- 1000
  s <- simulate(mle,nsim=n,seed=100)
  print("succ")
  xa <- c(0,0.34,0.64,0.96)
  if(flag==0)
  ya <- c(0,0.98)[2:1]
  if(flag==1)
  ya <- c(0,0.55)[2:1]


for(i in c(1:3))
{
j<-i
  par(fig=c(xa[j],xa[j+1],ya[2],ya[1]),
  mar=c(3,ifelse(j==1,4,1),1.5,ifelse(j==3,3,2)),new=T)
  times <- mle@data[1,]
#print(times)
  cases <- mle@data[1+i,]
  simul <- apply(sapply(s,function(x)x@data[1+i,]),1,median)
  s1 <- apply(sapply(s,function(x)x@data[1+i,]),1,function(x)quantile(x,0.025))
  s2 <- apply(sapply(s,function(x)x@data[1+i,]),1,function(x)quantile(x,0.975))

  print(i)

  matplot(times,(cbind(cases,simul,s1,s2)),type='l',xlim=d1,ylim=c(0,90),
  ann=F,axe=F,frame=T,xlab='',ylab='')
  polygon(c(times,times[length(times):1]),(c(s1,s2[length(s2):1])),col='lightgrey',border=NA)
time2<-times
  points(times,(cases),type='b',col='black',lwd=2)
  lines(times,(simul),col='red',lwd=2)
  axis(1,at=month,lab=substr(months(monthlab),1,3),cex.axis=1,tck=-0.01,padj=-1.2)
 # if(j==3)axis(1,at=2014:2016,lab=2014:2016,padj=0.8)
 # if(j==1)axis(1,at=2014:2016,lab=2014:2016,padj=0.8)
mtext(side=1,line=2,adj=0.5,paste(2015+i,'-',2016+i))
  axis(2,hadj=0.6,tck=-0.02)

if(j==1) mtext(side=2,line=2.2,text=expression(Weekly~reported~cases),cex=1.2,las=0,col='black')
  time <- seq(0,1,len=y$nm[1])
 # slope <- a[length(a)]
 # a[length(a)]<-a[1]
  times <- seq(min(times),max(times),by=1/365.25)
  a <- as.numeric(y[1,match('log.muv',names(y))+c(1:y$nm[1])-1])
  my <-cubicspline(seq(0,1,len=y$nm[1]),a,(times-min(times))/(max(times)-min(times)),endp2nd=TRUE)
par(new=T)
  ch<-c('BSh.0','CSh.0','DSh.0')
  para<-y[1,]
	        m<-exp(my)/365

# ms <- myspline(cbind(seq(0,1,len=length(infectious[[i]])),infectious[[i]]/7),a=0,b=0) 
# my <- mypredict(ms,(times-min(times))/(max(times)-min(times)))
  my <-cubicspline(seq(0,1,len=length(infectious[[i]])),infectious[[i]]/7,(times-min(times))/(max(times)-min(times)),endp2nd=TRUE)
 print(mean(m))
  plot(times,m, type='l',lty=2,lwd=2,cex=0.85,xlim=d1,ylim=c(-20,35),col='blue',ann=F,axe=F,frame=F,xlab='',ylab='')
#  points(seq(min(times),max(times),len=y$nm[1]),exp(a),col='red',lwd=2)
lines(times,my*x$cmuv[which.min(x$bic)]/365,lwd=2,col='brown',lty=4)
print(mean(my*x$cmuv[which.min(x$bic)]/365))
days=floor(365.25*(times-2016.324-0.33));
print(days)
lines(times,temp[days,j]/100)
axis(4,tck=-0.02,hadj=0.2,at=c(0,10,20,30),lab=ifelse(j==3,T,F),col='blue',col.ticks='blue')
 # lines(times,R0,type='b',col='blue')
  if(j==3) mtext(side=4,line=2,adj=0.42,text=expression(lambda[r](t)),cex=1.2,las=0,col='blue')
  if(j==3) mtext(side=4,line=2,adj=0.5,text='&',cex=1.2,las=0)
  if(j==3) mtext(side=4,line=2,adj=0.59,text=expression(lambda[h](t)),cex=1.2,las=0,col='brown')
#  R0 <- rbind(R0,c(mean(exp(my)/y[i,'gamma'])*(1+y[i,'phi']),y[i,'S.0']))
  mtext(side=3,line=0,adj=0,text=paste('(',letters[j+ifelse(flag==1,1,0)],')',z[i]),font=1,cex=1.2)

if(j==1){      col<-c('black','red','blue','brown')
        legend("topright",bty='n',lty=c(1,1,2,4),lwd=c(1,1,2,2),pch=c(21,NA,NA,NA),pt.bg=c('white',NA,NA,NA),cex=1,pt.cex=1,col=col,text.col=col,
        legend=c('reported cases','simulation','',''))
text(mean(d1)+0.12,27,expression(lambda[r](t)),col='blue')
text(mean(d1)+0.12,24,expression(lambda[h](t)),col='brown')
}

if(j==2) {
	m <- c(1:10)
	m1 <- par("fig")
	par(new=T,fig=c(m1[1]+(m1[2]-m1[1])*0.46,m1[2]*1.006,m1[3]+(m1[4]-m1[3])*0.65,m1[4]),
		mar=par("mar")+c(0,0,0.2,0.2))

	aic<-sapply(split(x$bic,x$nm),min)
	plot(as.numeric(names(aic)),aic,log='',type='b',col='black',pch=1,ann=F,axe=F,frame=T,xlab='',ylab='')
	axis(1,padj=-0.8,tck=-0.05)
	at <- axTicks(side=2)
	if(length(at)>4) at <- at[seq(1,length(at),by=2)]
	axis(2,at=at,hadj=0.8,tck=-0.05)
	m2 <- which.min(bic[m])
	points(x$rho[m[m2]],bic[m[m2]],col='red',pch=24)
	mtext(side=1,line=2,adj=0.5,text=expression(n[z]))
	mtext(side=2,line=2.7,text='AICc',las=0)
	}
}
print(y[1,])

print(min(bic))

  dev.off()
#  write.csv(R0,file='R0.csv')
  #system(paste('scp',file,'madhhe@158.132.174.118:'))
