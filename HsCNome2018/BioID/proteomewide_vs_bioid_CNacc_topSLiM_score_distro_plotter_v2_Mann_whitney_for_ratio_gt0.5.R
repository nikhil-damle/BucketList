#prot<-read.table("CNacc_proteome_CNacc_PXIXIT_shortIU_me0_topScorer")
##prot<-read.table("CNacc_proteome_CNacc_LXVP_shortIU_me0_topScorer")
#dr<-density(prot$V6,bw="SJ")
#xr<-hist(prot$V6,breaks=1000)    # Activate for Px
##xr<-hist(prot$V6,breaks=50)      # Activate for Lx
#xr$counts=log10(xr$counts)
#print("Dealing with IU >= 0")
##print(summary(prot$V6))
#
#bioid<-read.table("SWATH_n_DDA_WT_BFDR_le0.05_CNacc_PXIXIT_shortIU_me0")
##bioid<-read.table("SWATH_n_DDA_WT_BFDR_le0.05_CNacc_LXVP_shortIU_me0")
#d5<-density(bioid$V6,bw="SJ")
#x5<-hist(bioid$V6,breaks=100)
#x5$counts=log10(x5$counts)
#print("Dealing with IU >= 0")
##print(summary(bioid$V3))
#
#pdf("bioID_vs_entire_proteome_top_CNaccPx_Mann_whitney_IUme0_shift_hist.pdf",width=6,height=6,useDingbats=F)
##pdf("bioID_vs_entire_proteome_top_CNaccLx_Mann_whitney_IUme0_shift_hist.pdf",width=6,height=6,useDingbats=F)
#plot(xr,ylim=c(0,4),xlim=c(-50,70),xlab="PxIxIT score",cex.axis=1.5,cex.lab=1.5,ylab="Log10 (Frequency)",cex=1.2)
##plot(xr,ylim=c(0,4),xlim=c(-50,70),xlab="LxVP score",cex.axis=1.5,cex.lab=1.5,ylab="Log10 (Frequency)",cex=1.2)
#par(new=T)
#plot(x5,ylim=c(0,4),xlim=c(-50,70),xlab="",cex.axis=1.5,cex.lab=1.5,ylab="",col=rgb(1,0,0),cex=1.2)
#abline(v=mean(prot$V6),lty=2,lwd=1.5,col="black")
#abline(v=mean(bioid$V6),col="red",lty=2,lwd=1.5)
#box(which="plot",lty=1)
#dev.off()
#
#print(wilcox.test(prot$V6,bioid$V6,alternative="less",correct=TRUE)) # alt H1 => x < y
#print(wilcox.test(prot$V5,bioid$V5,correct=TRUE))

for(i in c(0, 0.1, 0.2, 0.3))
{
 #prot<-read.table(paste("CNacc_proteome_CNacc_PXIXIT_shortIU_me",i,"_topScorer",sep=""),sep="\t")
 prot<-read.table(paste("CNacc_proteome_CNacc_LXVP_shortIU_me",i,"_topScorer",sep=""),sep="\t")
 dr<-density(prot$V6,bw="SJ")
 #xr<-hist(prot$V6,breaks=1500) # Activate for Px
 xr<-hist(prot$V6,breaks=50)    # Activate for Lx
 xr$counts=log10(xr$counts)
 print(paste("Dealing with IU >= ",i))
 print(summary(prot$V6))

 #bioid<-read.table(paste("SWATH_n_DDA_WT_BFDR_le0.05_ratio_gt0.5_CNacc_PXIXIT_shortIU_me",i,sep=""),sep="\t")
 bioid<-read.table(paste("SWATH_n_DDA_WT_BFDR_le0.05_ratio_gt0.5_CNacc_LXVP_shortIU_me",i,sep=""),sep="\t")
 d5<-density(bioid$V6,bw="SJ")
 x5<-hist(bioid$V6,breaks=100)
 x5$counts=log10(x5$counts)
 print(paste("Dealing with IU >= ",i))
 print(summary(bioid$V6))

 #pdf(paste("bioID_vs_entire_proteome_top_CNaccPx_ratio_gt0.5_Mann_whitney_IUme",i,"_shift_hist.pdf",sep=""),width=6,height=6,useDingbats=F)
 pdf(paste("bioID_vs_entire_proteome_top_CNaccLx_ratio_gt0.5_Mann_whitney_IUme",i,"_shift_hist.pdf",sep=""),width=6,height=6,useDingbats=F)
 #plot(xr,ylim=c(0,4),xlim=c(-50,70),xlab="PxIxIT score",cex.axis=1.5,cex.lab=1.5,ylab="Log10 (Frequency)",cex=1.2)
 plot(xr,ylim=c(0,4),xlim=c(-50,70),xlab="LxVP score",cex.axis=1.5,cex.lab=1.5,ylab="Log10 (Frequency)",cex=1.2)
 par(new=T)
 plot(x5,ylim=c(0,4),xlim=c(-50,70),xlab="",cex.axis=1.5,cex.lab=1.5,ylab="",col=rgb(1,0,0),cex=1.2)
 abline(v=mean(prot$V6),lty=2,lwd=1.5,col="black")
 abline(v=mean(bioid$V6),col="red",lty=2,lwd=1.5)
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
 print(wilcox.test(prot$V6,bioid$V6,alternative="less",correct=TRUE)) # alt H1 => x < y
 #print(median.test(prot$V6,bioid$V6,correct=TRUE))
 #print(wilcox.test(prot$V5,bioid$V5,correct=TRUE))
}
