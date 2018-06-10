prot<-read.table("top_scoring_disordered_PXIXIT_IUme0")
#prot<-read.table("top_scoring_disordered_LXVP_IUme0")
dr<-density(prot$V4,bw="SJ")
xr<-hist(prot$V4,breaks=100)    # Activate for Px
#xr<-hist(prot$V4,breaks=50)      # Activate for Lx
xr$counts=log10(xr$counts)
print("Dealing with IU >= 0")
#print(summary(prot$V4))

bioid<-read.table("bioid_top_scoring_disordered_PXIXIT_IUme0")
#bioid<-read.table("bioid_top_scoring_disordered_LXVP_IUme0")
d5<-density(bioid$V4,bw="SJ")
x5<-hist(bioid$V4,breaks=100)
x5$counts=log10(x5$counts)
print("Dealing with IU >= 0")
#print(summary(bioid$V3))

pdf("bioID_vs_entire_proteome_topPx_Mann_whitney_IUme0_shift_hist.pdf",width=6,height=6,useDingbats=F)
#pdf("bioID_vs_entire_proteome_topLx_Mann_whitney_IUme0_shift_hist.pdf",width=6,height=6,useDingbats=F)
plot(xr,ylim=c(0,4),xlim=c(-50,70),xlab="PxIxIT score",cex.axis=1.5,cex.lab=1.5,ylab="Log10 (Frequency)",cex=1.2)
#plot(xr,ylim=c(0,4),xlim=c(-50,70),xlab="LxVP score",cex.axis=1.5,cex.lab=1.5,ylab="Log10 (Frequency)",cex=1.2)
par(new=T)
plot(x5,ylim=c(0,4),xlim=c(-50,70),xlab="",cex.axis=1.5,cex.lab=1.5,ylab="",col=rgb(1,0,0),cex=1.2)
abline(v=mean(prot$V4),lty=2,lwd=1.5,col="black")
abline(v=mean(bioid$V4),col="red",lty=2,lwd=1.5)
box(which="plot",lty=1)
dev.off()

print(wilcox.test(prot$V4,bioid$V4,alternative="less",correct=TRUE)) # alt H1 => x < y
print(wilcox.test(prot$V5,bioid$V5,correct=TRUE))

for(i in c(0.1, 0.2, 0.3))
{
 prot<-read.table(paste("top_scoring_disordered_PXIXIT_IUme",i,sep=""),sep="\t")
 #prot<-read.table(paste("top_scoring_disordered_LXVP_IUme",i,sep=""),sep="\t")
 dr<-density(prot$V4,bw="SJ")
 xr<-hist(prot$V4,breaks=1500) # Activate for Px
 #xr<-hist(prot$V4,breaks=50)    # Activate for Lx
 xr$counts=log10(xr$counts)
 print(paste("Dealing with IU >= ",i))
 #print(summary(prot$V4))

 bioid<-read.table(paste("bioid_top_scoring_disordered_PXIXIT_IUme",i,sep=""),sep="\t")
 #bioid<-read.table(paste("bioid_top_scoring_disordered_LXVP_IUme",i,sep=""),sep="\t")
 d5<-density(bioid$V4,bw="SJ")
 x5<-hist(bioid$V4,breaks=100)
 x5$counts=log10(x5$counts)
 print(paste("Dealing with IU >= ",i))
 #print(summary(bioid$V3))

 pdf(paste("bioID_vs_entire_proteome_topPx_Mann_whitney_IUme",i,"_shift_hist.pdf",sep=""),width=6,height=6,useDingbats=F)
 #pdf(paste("bioID_vs_entire_proteome_topLx_Mann_whitney_IUme",i,"_shift_hist.pdf",sep=""),width=6,height=6,useDingbats=F)
 plot(xr,ylim=c(0,4),xlim=c(-50,70),xlab="PxIxIT score",cex.axis=1.5,cex.lab=1.5,ylab="Log10 (Frequency)",cex=1.2)
 #plot(xr,ylim=c(0,4),xlim=c(-50,70),xlab="LxVP score",cex.axis=1.5,cex.lab=1.5,ylab="Log10 (Frequency)",cex=1.2)
 par(new=T)
 plot(x5,ylim=c(0,4),xlim=c(-50,70),xlab="",cex.axis=1.5,cex.lab=1.5,ylab="",col=rgb(1,0,0),cex=1.2)
 abline(v=mean(prot$V4),lty=2,lwd=1.5,col="black")
 abline(v=mean(bioid$V4),col="red",lty=2,lwd=1.5)
 box(which="plot",lty=1)
 dev.off()

 #pdf(paste("bioID_vs_entire_proteome_topPx_Mann_whitney_IUme",i,"_shift_density.pdf",sep=""),width=6,height=6,useDingbats=F)
 #pdf(paste("bioID_vs_entire_proteome_topLx_Mann_whitney_IUme",i,"_shift_density.pdf",sep=""),width=6,height=6,useDingbats=F)
 #plot(dr$x,dr$y,col="black",xlim=c(-25,50),type="l",xlab="PxIxIT score",ylab="Density",cex.axis=1.5,cex.lab=1.5,cex=1.2)
 #plot(dr$x,dr$y,col="black",xlim=c(-25,50),type="l",xlab="LxVP score",ylab="Density",cex.axis=1.5,cex.lab=1.5,cex=1.2)
 #lines(d5$x,d5$y,col=rgb(1,0,0))
 #box(which="plot",lty=1)
 #dev.off()

 print(paste("Dealing with IU >= ",i))
 # Null H0 => diff in means of the two distros = 0
 print(wilcox.test(prot$V4,bioid$V4,alternative="less",correct=TRUE)) # alt H1 => x < y
 print(wilcox.test(prot$V5,bioid$V5,correct=TRUE))
}
