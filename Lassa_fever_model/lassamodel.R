logit <- function(p){log(p/(1-p))}
expit <- function(r){1/(1+exp(-r))}
expp<-function(gamma,n1,dt){365*dt*n1/(1.0-exp(-n1*gamma*dt))}
logg<-function(gamma,n1,dt){-log(1.0-(365*dt*n1)/gamma)/(n1*dt)}


par.trans <- function (params) {
  x <- params
  r <- length(dim(x))
  if (r > 1) {
    x <- apply(x,2:r,par.trans)
  } else {
    x['beta1'] <- log(params['beta1'])
    x['beta2'] <- log(params['beta2'])
    x['beta3'] <- log(params['beta3'])
    x['muh'] <- log(params['muh'])
    x['sigma1'] <- log(params['sigma1'])
    x['gamma1'] <- log(params['gamma1'])
    x['tau1'] <- log(params['tau1'])
    x['alpha1'] <- log(params['alpha1'])
    x['delta1'] <- log(params['delta1'])
    x['sigma2'] <- log(params['sigma2'])
    x['tau2'] <- log(params['tau2'])
    x['alpha2'] <- log(params['alpha2'])
    x['delta2'] <- log(params['delta2'])
    x['tau3'] <- log(params['tau3'])
    x['delta3'] <- log(params['delta3'])
    x['xi'] <- log(params['xi'])
    x['rho'] <- logit(params['rho'])
    x['t0'] <- logit(params['t0'])
    x['tau'] <- log(params['tau'])
    x['bmuv'] <- log(params['bmuv'])
    x['cmuv'] <- log(params['cmuv'])
  }
  x
}

par.untrans <- function (params) {
  x <- params
  r <- length(dim(x))
  if (r > 1) {
    x <- apply(x,2:r,par.untrans)
  } else {
    x['beta1'] <- exp(params['beta1'])
    x['beta2'] <- exp(params['beta2'])
    x['beta3'] <- exp(params['beta3'])
    x['muh'] <- exp(params['muh'])
    x['sigma1'] <- exp(params['sigma1'])
    x['gamma1'] <- exp(params['gamma1'])
    x['tau1'] <- exp(params['tau1'])
    x['alpha1'] <- exp(params['alpha1'])
    x['delta1'] <- exp(params['delta1'])
    x['sigma2'] <- exp(params['sigma2'])
    x['tau2'] <- exp(params['tau2'])
    x['alpha2'] <- exp(params['alpha2'])
    x['delta2'] <- exp(params['delta2'])
    x['tau3'] <- exp(params['tau3'])
    x['delta3'] <- exp(params['delta3'])
    x['xi'] <- exp(params['xi'])
    x['rho'] <- expit(params['rho'])
    x['t0'] <- expit(params['t0'])
    x['tau'] <- exp(params['tau'])
    x['bmuv'] <- exp(params['bmuv'])
    x['cmuv'] <- exp(params['cmuv'])
  }
  x
}

district.pomp <- function (lassacases, temp1, temp2, temp3, paramnames, types, nb) {

  pomp(
       nstep=7,
       data=rbind(
         time=lassacases$time,
         cases.b=lassacases[,'B'],
         cases.c=lassacases[,'C'],
         cases.d=lassacases[,'D']
         ),
       t0=as.double(2*lassacases$time[1]-lassacases$time[2]),
       rprocess = function (X, t1, t2, nstep, tbp, bp, tbasis, basis, ...) {
         nvar <- nrow(X)
         np <- ncol(X)
         stateindex <- match(c(
	'BSH','BEH','BIH','BEQ','BIQ','BH','BR','BC',
	'CSH','CEH','CIH','CEQ','CIQ','CH','CR','CC',
	'DSH','DEH','DIH','DEQ','DIQ','DH','DR','DC'
	),rownames(X))-1
         parindex <- match(c('beta1','beta2','beta3','muh','sigma1','gamma1','tau1','alpha1','delta1','sigma2',
'tau2','alpha2','delta2','tau3','delta3','xi','rho','nm',
               'log.muv','bmuv','cmuv','pop','t0'),rownames(X))-1
         result <- .C("basic_sirs_pois",
                      X = as.double(X),
                      t1 = as.double(t1),
                      t2 = as.double(t2),
                      nstep = as.integer(nstep),
		      t_start=as.double(2*lassacases$time[1]-lassacases$time[2]),
		      t_end=as.double(lassacases$time[length(lassacases$time)]),
                      nvar = as.integer(nvar),
                      np = as.integer(np),
                      stateindex = as.integer(stateindex),
                      parindex = as.integer(parindex),
		      temp1 =  as.double(temp1),
		      temp2 =  as.double(temp2),
		      temp3 =  as.double(temp3),
		      seas_dim = as.integer(nb),
                      NAOK = TRUE,
                      DUP = FALSE,
                      PACKAGE = "lassamodel"
                      )$X
         array(result,dim=dim(X),dimnames=dimnames(X))
       },
       dmeasure = function (X, Y, ...) {
         n <- dim(X)
         np <- n[2]
         lassaindex <- match(c('BC','CC','DC','tau'),rownames(X))-1
         .C("negbin_dmeasure",
            n = as.integer(n),
            index = as.integer(lassaindex),
            X = as.double(X),
            y1 = as.double(Y[2]),
            y2 = as.double(Y[3]),
            y3 = as.double(Y[4]),
            f = double(np),
            DUP = FALSE,
            NAOK = TRUE,
            PACKAGE = "lassamodel"
            )$f
       },
       rmeasure = function (X, time, ...) {
         n <- dim(X)
         nv <- n[1]
         np <- n[2]
         lassaindex <- match(c('BC','CC','DC','tau'),rownames(X))-1
         matrix(
                .C("negbin_rmeasure",
                   n = as.integer(n),
                   index = as.integer(lassaindex),
                   X = as.double(X),
                   cases = double(3*np),   ##################################   number of columns of data
                   DUP = FALSE,
                   NAOK = TRUE,
                   PACKAGE = "lassamodel"
                   )$cases,
                3,    #################### number of columns of data
                np
                )
       },
       particles = function (Np, center, sd, parnames, fixnames, ivpnames, statenames, bp, ...) {
         pop <- rep(center['pop'], 3)
         infe.ini <- center['infe.ini'] 
         susc.ini <- center['susc.ini'] 
	 neqn <- 7
	 ncity <- 3
	 names.cases <- c('BC','CC','DC')
         X <- matrix(
                     data=0,
                     nrow=length(names.cases)+length(statenames)+length(parnames)+length(fixnames)+length(ivpnames),
			#### add number of columns of data
                     ncol=Np,
                     dimnames=list(
                       c(statenames,parnames,fixnames,ivpnames,names.cases),
                       NULL
                       )
                     )
#	a<-20
#        center[grep('log.muv',names(center))]<-ifelse(center[grep('log.muv',names(center))]>log(a),log(a),center[grep('log.muv',names(center))])
#	center['rho']<-ifelse(center['rho']>logit(0.8),logit(0.8),center['rho'])

         X[parnames,] <- rnorm(
                               n=Np*length(parnames),
                               mean=center[parnames],
                               sd=sd[parnames]
                               )
         X[fixnames,] <- center[fixnames]
         
         X[ivpnames,] <- apply(
                               matrix(
                                      center[ivpnames]*rlnorm(
                                                              n=Np*length(ivpnames),
                                                              meanlog=-0.5*sd[ivpnames]^2,
                                                              sdlog=sd[ivpnames]
                                                              ),
                                      nrow=length(ivpnames),
                                      ncol=Np
                                      ),
                               2,
                               function(x)
                                 {
				for(i in 1:ncity) {
                                x[1+(i-1)*neqn]<-ifelse(x[1+(i-1)*neqn]<0.01,0.01,x[1+(i-1)*neqn])
                                x[1+(i-1)*neqn]<-ifelse(x[1+(i-1)*neqn]>0.9,0.9,x[1+(i-1)*neqn])
			#	if(i==1) x[1+(i-1)*neqn]<-0.95-susc.ini
			#	else x[1+(i-1)*neqn]<-0.95
                                x[2+(i-1)*neqn]<-ifelse(x[2+(i-1)*neqn]<infe.ini/pop[i],infe.ini/pop[i],x[2+(i-1)*neqn])
                                x[2+(i-1)*neqn]<-ifelse(x[2+(i-1)*neqn]>10*infe.ini/pop[i],10*infe.ini/pop[i],x[2+(i-1)*neqn])
                                x[3+(i-1)*neqn]<-ifelse(x[3+(i-1)*neqn]<infe.ini/pop[i],infe.ini/pop[i],x[3+(i-1)*neqn])
                                x[3+(i-1)*neqn]<-ifelse(x[3+(i-1)*neqn]>10*infe.ini/pop[i],10*infe.ini/pop[i],x[3+(i-1)*neqn])
				x[c(4:6)+(i-1)*neqn]<-0
				x[neqn+(i-1)*neqn]<-1-sum(x[1:(neqn-1)+(i-1)*neqn])
				}
                                x
                                }
                               )
  	for(i in 1:ncity)
         X[statenames[c(1:neqn)+(i-1)*neqn],] <- round(pop[i]*X[ivpnames[c(1:neqn)+(i-1)*neqn],]) # must round to nearest integer!
         X
       },
       parnames=paramnames[types=='par'],
       ivpnames=paramnames[types=='ivp'],
       fixnames=paramnames[types=='fixed'],
       statenames=sub('.0','',paramnames[types=='ivp'])
       )
}
