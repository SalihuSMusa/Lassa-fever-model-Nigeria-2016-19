Npar<-2000
require(pomp.orig)

Args <- commandArgs(trailingOnly=TRUE)
print(Args)
ci<-eval(parse(text = Args))
nreps <- 1
#imagefile <- sprintf("serljob%3.3d.rda", ci) 
id <- 'fit1f1_lassa'
mpfile <- 'mp_lassa.csv'
bestfile <- sprintf("best_lassa%3.3d.csv", ci)  
mlefile <- sprintf("mle_lassa%3.3d.rda", ci) 
filtnames <- c(
'BSH','BEH','BIH','BEQ','BIQ','BH','BR','BC',
'CSH','CEH','CIH','CEQ','CIQ','CH','CR','CC',
'DSH','DEH','DIH','DEQ','DIQ','DH','DR','DC'
)
#pfilt<-TRUE
pfilt<-FALSE


x<-read.csv('nigeria_LHF.csv',head=1)[,c(1,3,5,2)]
lassacases<-as.data.frame(x[,-4])
names(lassacases)<-c('B','C','D')
lassacases$time<-x[,4]
print(lassacases)
x<-read.csv('precip_nigeria.csv')
xout<-(1:1156)/365.25+2016
precip<-approx(x[,1],y=smooth(x[,7]),xout)$y
temp1<-precip[which(xout>2016.34-1/12 & xout<2016.931)]
temp2<-precip[which(xout>2017.34-1/12 & xout<2017.931)]
temp3<-precip[which(xout>2018.34-1/12 & xout<2018.931)]


#print(lassacases)

parameters <- read.csv(file=mpfile,row.names=1,as.is=TRUE)
districts <- names(parameters)[-c(1,2)]
districts <- districts[ifelse(ci %% 10 ==0, 10, ci%%10)] 

print(districts)
#print(cbind(rownames(parameters),parameters[,districts]))
ifelse(pfilt==TRUE, test.sigma <- parameters$sd/1e20, test.sigma <- parameters$sd)
names(test.sigma) <- rownames(parameters)

source('lassamodel.R')


parameters[grep('log.muv',rownames(parameters)),'type']<-'fixed'
n<-as.numeric(parameters['nm',districts])
parameters[grep('log.muv',rownames(parameters)),'type'][1:n]<-'par'
parameters[grep('log.muv',rownames(parameters)),districts][-c(1:n)]<-1


print(parameters[,c('type',districts)])

ndists <- length(districts)
njobs <- nreps*ndists
#seeds <- lavaseed(n=njobs)
seeds<-ceiling(10^runif(n=njobs,min=7,max=9))
tic<-Sys.time()
joblist <- vector(mode='list',length=njobs)
j <- 0
for (d in districts) {
  mle <- district.pomp(lassacases,temp1,temp2,temp3,rownames(parameters),parameters$type,parameters['nm',d])
  test.params <- parameters[[d]]
  names(test.params) <- rownames(parameters)
  theta.x <- as.matrix(particles(mle,Np=nreps,center=par.trans(test.params),sd=test.sigma)[rownames(parameters),])
  for (r in 1:nreps) {
    j <- j+1
    x <- theta.x[,r]
    names(x) <- rownames(theta.x)
	print(names(x))
    joblist[[j]] <- list(
                         mle=mle,
                         district=d,
                         seed=seeds[j],
                         theta=x,
                         sigma=test.sigma,
                         ivpnames= rownames(parameters)[parameters$type=='ivp'],
                         estnames=rownames(parameters)[parameters$type=='par'],
                         filtnames=filtnames
                         )
  }
}
dyn.load('lassamodel.so')
done <- 0
for(i in 1:nreps){
joblist[[i]] <-
with(joblist[[i]],
                    {
                      save.seed <- .Random.seed
                      set.seed(seed)
                      x <- mif(
                               mle,
                               Nmif=1,
                               ivps=ivpnames,
                               pars=estnames,
                               stvs=filtnames,
                               alg.pars=list(Np=Npar,CC=3,T0=nrow(lassacases),cooling.factor=0.95),
                               start=theta,
                               rw.sd=sigma,
                               max.fail=1000,
                               weighted=T,
                               warn=F
                               )
                      result <- list(
                                     mle=x,
                                     district=district,
                                     seed=seed,
                                     rngstate=.Random.seed
                                     )
                      .Random.seed <<- save.seed
                      result
                    }
                    )
}

joblist <- joblist[!sapply(joblist,is.null)]
#save.image(file=imagefile)
done <- done+1

while (done < ifelse(pfilt==TRUE,1,10)) {
for(i in 1:nreps){
joblist[[i]] <-
with(joblist[[i]],
                    {
                        save.seed <- .Random.seed
                        .Random.seed <<- rngstate
                        if (done < 2) {
                          x <- continue(mle,Nmif=11,weighted=T,max.fail=1000)
                        } else {
                          x <- continue(mle,Nmif=11,weighted=T,max.fail=1000)
                        }
                        result <- list(
                                       mle=x,
                                       district=district,
                                       seed=seed,
                                       rngstate=.Random.seed
                                       )
                        .Random.seed <<- save.seed
                        result
                      }
                      )
}
  joblist <- joblist[!sapply(joblist,is.null)]
 # save.image(file=imagefile)
  done <- done+1
}

for(i in 1:nreps){
joblist[[i]] <-
with(joblist[[i]],
                    {
                      save.seed <- .Random.seed
                      .Random.seed <<- rngstate
                      ff <- vector(mode='list')
                      for (n in 1:10)
                        ff[[n]] <- pfilter(mle,max.fail=1000)
                      result <- list(
                                     mle=mle,
                                     district=district,
                                     seed=seed,
                                     rngstate=rngstate,
                                     nfail=sapply(ff,function(x)x$nfail),
                                     loglik=sapply(ff,function(x)x$loglik)
                                     )
                      .Random.seed <<- save.seed
                      result
                    },
                    joblist=joblist
                    )
}
joblist <- joblist[!sapply(joblist,is.null)]
#save.image(file=imagefile)


x <- do.call(
             rbind,
             lapply(
                    joblist,
                    function (x) {
                      cbind(
                            data.frame(
                                       district=x$district,
                                       loglik=mean(x$loglik),
                                       loglik.sd=sd(x$loglik),
                                       nfail.max=max(x$nfail),
                                       nfail.min=min(x$nfail)
                                       ),
                            as.data.frame(t(par.untrans(coef(x$mle))))
                            )
                    }
                    )
             )

best.ind <- tapply(1:nrow(x),x$district,function(k)k[which.max(x$loglik[k])])
best <- x[best.ind,]
rownames(best) <- as.character(best$district)
best$district <- NULL
best <- best[order(rownames(best)),]
write.csv(best,file=bestfile)

mle <- lapply(
              joblist[best.ind],
              function(x)x$mle
              )
names(mle) <- sapply(joblist[best.ind],function(x)x$district)
mle <- mle[order(names(mle))]
save(list='mle',file=mlefile)
toc <- Sys.time()
print(toc-tic)


dyn.unload('lassamodel.so')



quit()

