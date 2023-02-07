require('pomp.orig')
all<-list(len=10)
j<-1
for(i in 1:10)
{
if(file.exists(sprintf('mle_lassa%3.3d.rda',i))){
load(sprintf('mle_lassa%3.3d.rda',i))
all[[j]]<-mle[[1]]
j<-j+1
}
}
save(all,file='allregion.rda')
