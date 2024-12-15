#### Fonctions de base #### 
library(FreeSortR)
MdsDiss= function (MatDissimil, ndim = 2, metric = TRUE, ties = "primary", 
          itmax = 5000, eps = 1e-06) 
{
  res <- smacof::smacofSym(MatDissimil, ndim = ndim, type = "interval", 
                   verbose = FALSE, itmax = itmax, eps = eps)
  Config <- scale(res$conf, center = TRUE, scale = FALSE)
  W <- Config %*% t(Config)
  bid <- svd(W)
  if(ndim==1){
    Config <- as.matrix(bid$u[, 1:ndim]) %*% sqrt(diag(as.matrix(bid$d[1:ndim])))}else{
      Config <- bid$u[, 1:ndim] %*% sqrt(diag(bid$d[1:ndim])) 
    }
  if (metric == FALSE) {
    res <- smacof::smacofSym(MatDissimil, ndim = ndim, type = "ordinal", 
                     ties = ties, init = Config, verbose = FALSE, itmax = itmax, 
                     eps = eps)
    Config <- scale(res$conf, center = TRUE, scale = FALSE)
    W <- Config %*% t(Config)
    bid <- svd(W)
    bid$d[1:ndim]/sum(bid$d)
    if(ndim==1){
      Config <- as.matrix(bid$u[, 1:ndim]) %*% sqrt(diag(as.matrix(bid$d[1:ndim])))}else{
        Config <- bid$u[, 1:ndim] %*% sqrt(diag(bid$d[1:ndim])) 
      }
  }
  Percent <- bid$d/sum(bid$d)
  Stress <- sqrt(sum((res$dhat - res$confdist)^2)/sum(res$confdist^2))
  res <- list(Config = Config, Percent = Percent, Stress = Stress)
  return(res)
}

MdsSort=function (Part, ndim = 2, nboot = 0, metric = FALSE, ties = "primary", 
          itmax = 5000, eps = 1e-06) 
{
  if (!class(Part) == "SortingPartition") {
    return("This is not an object of class SortingPartition")
  }
  else {
    nstimuli <- Part@nstimuli
    nsubjects <- Part@nsubjects
    ListDissim <- Dissimil(Part)
    MatDissim <- apply(simplify2array(ListDissim), c(1, 2), 
                       "sum")
    res <- MdsDiss(MatDissim, ndim = ndim, metric = metric, 
                   ties = ties, itmax, eps)
    Config <- res$Config
    rownames(Config) <- Part@LabStim
    colnames(Config) <- paste(rep("Dim"), 1:ndim)
    Stress <- res$Stress
    Percent <- res$Percent[1:ndim]
    if (nboot > 0) {
      ResBoot <- vector("list", nstimuli)
      for (b in 1:nboot) {
        panel <- sample(1:nsubjects, replace = TRUE)
        MatDistTotPanel = matrix(0, nstimuli, nstimuli)
        for (k in 1:nsubjects) {
          MatDistTotPanel <- MatDistTotPanel + ListDissim[[panel[k]]]
        }
        if (metric == TRUE) {
          res <- smacof::smacofSym(MatDistTotPanel, ndim = ndim, 
                           type = "ratio", ties = "primary", init = Config, 
                           verbose = F, itmax = itmax, eps = eps)
        }
        else {
          res <- smacof::smacofSym(MatDistTotPanel, ndim = ndim, 
                           type = "ordinal", ties = "primary", init = Config, 
                           verbose = F, itmax = itmax, eps = eps)
        }
        ConfigP <- res$conf
        res <- procrustes(Config, ConfigP)
        for (i in 1:nstimuli) {
          ResBoot[[i]] <- rbind(ResBoot[[i]], res$Yrot[i, 
          ])
        }
      }
      return(new("SortingMds", nstimuli = nstimuli, nsubjects = nsubjects, 
                 LabStim = Part@LabStim, LabSubj = Part@LabSubj, 
                 ndim = ndim, Config = Config, Percent = Percent, 
                 Stress = Stress, ResBoot = ResBoot))
    }
    else {
      ResBoot <- NULL
      return(new("SortingMds", nstimuli = nstimuli, nsubjects = nsubjects, 
                 LabStim = Part@LabStim, LabSubj = Part@LabSubj, 
                 ndim = ndim, Config = Config, Percent = Percent, 
                 Stress = Stress))
    }
  }
}


plot_MDS=function(resMds,j,k){
  N=resMds@nstimuli
  couleur=rainbow(N)
  df=as.data.frame(resMds@ResBoot[[1]])[,c(j,k)]
  colnames(df)=c("x","y")
  p<-ggplot(df,aes(x=x,y=y))+stat_ellipse(color=couleur[1])
  
  for(i in 2:N){
    df=as.data.frame(resMds@ResBoot[[i]])[,c(j,k)]
    colnames(df)=c("x","y")
    p<-p+stat_ellipse(data=df,aes(x=x,y=y),color=couleur[i])
  }
  config2<-data.frame(matrix(NA,nrow=nrow(resMds@Config),ncol=2))
  rownames(config2)=rownames(resMds@Config)
  for(i in 1:N){
    config2[i,]=t(apply(resMds@ResBoot[[i]][,c(j,k)],2,mean))
  }
  colnames(config2)=c("Dim1","Dim2")
  Pourcentage=round(resMds@Percent[c(j,k)]*100,1)
  p<-p+geom_text(data=config2,aes(x=Dim1,y=Dim2,
                                  label=rownames(config2)),color=couleur)+
    theme_minimal()+
    xlab(paste("Dim",j,Pourcentage[1],"%",sep=" "))+
    ylab(paste("Dim",k,Pourcentage[2],"%",sep=" "))+
    ggtitle(paste("MDS, Stress=",round(resMds@Stress,3),sep=" "))+
    theme(legend.position = "n")
  return(p)
}

#### MDS ####

data("AromaSort")

## La base est notée df
df=AromaSort
View(df)

## Définition de l'objet de type tri-libre

tl<-SortingPartition(df)

# contenu de tl :

tl@nsubjects
tl@Partition

DescriptionPartition(tl,subject=1)

ListDiss<-Dissimil(tl)
T1=as.data.frame(ListDiss[[1]])
rownames(T1)=tl@LabStim
colnames(T1)=tl@LabStim

D<-DissTot(tl)
rownames(D)=tl@LabStim

stress<-rep(NA,6)
# getStress donne le stress de Kruskal (ie normalisé)

for(i in 1:6){
  resMds<-MdsSort(tl,ndim=i,metric = T)
  stress[i]<-getStress(resMds)  
}

plot(1:6,stress,type="b",pch=20,
     xlab="Nombre de dimensions",
     col="purple")


ndim=3
resMds<-MdsSort(tl,ndim=ndim,metric = T)
config<-getConfig(resMds)
config=as.data.frame(config)
colnames(config)=paste("Dim",1:ndim,sep="")
p1<-ggplot(config,aes(x=Dim1,y=Dim2,label=rownames(config)))+geom_text()+
  ggtitle(paste("MDS, Stress=",round(resMds@Stress,3)))+
  xlab(paste("Dim1",round(100*resMds@Percent[1],1),"%",sep=" "))+
  ylab(paste("Dim2",round(100*resMds@Percent[2],1),"%",sep=" "))+
  theme_minimal()

p2<-ggplot(config,aes(x=Dim2,y=Dim3,label=rownames(config)))+geom_text()+
  ggtitle(paste("MDS, Stress=",round(resMds@Stress,3)))+
  xlab(paste("Dim2",round(100*resMds@Percent[2],1),"%",sep=" "))+
  ylab(paste("Dim3",round(100*resMds@Percent[3],1),"%",sep=" "))+
  theme_minimal()

library(gridExtra)
grid.arrange(p1,p2,ncol=2)

# Ellipses

resMds<-MdsSort(tl,ndim=ndim,metric = T,nboot=500)
p1=plot_MDS(resMds,1,2)
p2=plot_MDS(resMds,1,3)
grid.arrange(p1,p2,ncol=2)

### Informations additionnelles ###

data("AromaTerms")
ia=as.data.frame(AromaTerms)

beta<-data.frame(matrix(NA,nrow=dim(ia)[2],ncol=dim(config)[2]))

for(i in 1:dim(beta)[1]){
  beta[i,]<-coef(lm(AromaTerms[,i]~as.matrix(config)))[-1]
}
rownames(beta)<-colnames(ia)
colnames(beta)[1:ndim]<-paste("Dim",1:ndim,sep="")

ggplot(beta, aes(x=1.05*Dim1, y=1.05*Dim2,label=rownames(beta))) +
  geom_text() +
  geom_segment(data=beta,aes(x=0, y=0, xend=Dim1, yend=Dim2), arrow = arrow(length=unit(.1, 'cm')))+
  theme_minimal()+
  #theme(element_text(size=10))+
  xlab(paste("Dim1=",round(resMds@Percent[1]*100,3),"%"))+
  ylab(paste("Dim2=",round(resMds@Percent[2]*100,3),"%"))+
  ggtitle(paste("MDS, stress=",round(resMds@Stress,3)))


### Partition consensuelle ###

cons=ConsensusPartition(tl,ngroups=0,type="fusion")

# Contenu de cons :
cons$Crit

cons$Consensus

config=getConfig(resMds)
colnames(config)=paste("Dim",1:ndim,sep="")
Consensus=as.factor(cons$Consensus)
ggplot(config,aes(x=Dim1,y=Dim2,label=rownames(config),color=Consensus))+
  geom_text()+
  xlab(paste("Dim1=",round(resMds@Percent[1]*100,1),"%"))+
  ylab(paste("Dim2=",round(resMds@Percent[2]*100,1),"%"))+
  ggtitle(paste("MDS, stress=",round(resMds@Stress,3)))+
  theme_minimal()+theme(legend.position = "none")


d_ARI<-function(df,i,j){
  return(sqrt(1-RandIndex(df[,i],df[,j])$AdjustedRand))
}
N<-dim(df)[2]
M<-matrix(NA,nrow=N,ncol=N)
for(i in 1:N){
  for(j in 1:N){
    M[i,j]<-d_ARI(df,i,j)
  }
}

dd <- as.dist(M)

res<-hclust(dd,method="ward.D2")
plot(res,hang=-1,ylim=c(0,1.5))

k.gr=2
gr<-cutree(res,k.gr)


df.g<-df[,gr==1]
tl.g<-SortingPartition(df.g)
resMds.g<-MdsSort(tl.g,ndim=3)
config.g<-getConfig(resMds.g)
cons.g<-ConsensusPartition(tl.g,ngroups=0)
config.g=getConfig(resMds.g)
colnames(config.g)=paste("Dim",1:ndim,sep="")
Consensus.g=as.factor(cons.g$Consensus)

p1=ggplot(config.g,aes(x=Dim1,y=Dim2,label=rownames(config.g),color=Consensus.g))+
  geom_text()+
  xlab(paste("Dim1=",round(resMds.g@Percent[1]*100,3),"%"))+
  ylab(paste("Dim2=",round(resMds.g@Percent[2]*100,3),"%"))+
  labs(
    title=paste("MDS, stress=",round(resMds.g@Stress,3),", N=",sum(gr==2),sep=""),
    subtitle = paste("Consensus=",round(cons.g$Crit,3),", k gr=",max(cons.g$Consensus)),sep="")+
  theme_minimal()+theme(legend.position = "none")
