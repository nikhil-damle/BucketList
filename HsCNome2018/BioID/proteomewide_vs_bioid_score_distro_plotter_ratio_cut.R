for(i in c(0, 0.1, 0.2, 0.3))
{
 #ransam<-read.table(paste("mean10000X_PxPSSM_topScorer_distro_frm_REST_o_HsProteome_IUme",i,"_for_ratio_gt0.5",sep=""))
 ransam<-read.table(paste("mean10000X_LxPSSM_topScorer_distro_frm_REST_o_HsProteome_IUme",i,"_for_ratio_gt0.5",sep=""))
 dr<-density(ransam$V2,bw="SJ")
 xr<-hist(ransam$V2,breaks=15)
 xr$counts=log10(xr$counts)
 print(paste("Dealing with IU >= ",i))
 summary(ransam$V2)

 #bioid<-read.table(paste("SnD_highest_scoring_uniq12mer_IUme",i,"_ratio_gt0.5",sep=""))
 bioid<-read.table(paste("SnD_highest_scoring_uniq7mer_IUme",i,"_ratio_gt0.5",sep=""))
 d5<-density(bioid$V3,bw="SJ")
 x5<-hist(bioid$V3,breaks=50)
 x5$counts=log10(x5$counts)
 print(paste("Dealing with IU >= ",i))
 summary(bioid$V3)

 #pdf(paste("bioID_vs_rest_o_proteome_PxPSSM_IUme",i,"_shift_hist_ratio_gt0.5.pdf",sep=""),width=6,height=6,useDingbats=F)
 pdf(paste("bioID_vs_rest_o_proteome_LxPSSM_IUme",i,"_shift_hist_ratio_gt0.5.pdf",sep=""),width=6,height=6,useDingbats=F)
 #plot(xr,ylim=c(0,4),xlim=c(-10,50),xlab="PxIxIT score",cex.axis=1.5,cex.lab=1.5,ylab="Log10 (Frequency)",cex=1.2)
 plot(xr,ylim=c(0,4),xlim=c(-10,50),xlab="LxVP score",cex.axis=1.5,cex.lab=1.5,ylab="Log10 (Frequency)",cex=1.2)
 par(new=T)
 plot(x5,ylim=c(0,4),xlim=c(-10,50),xlab="",cex.axis=1.5,cex.lab=1.5,ylab="",col=rgb(1,0,0,0.3),cex=1.2)
 box(which="plot",lty=1)
 dev.off()

 #pdf(paste("bioID_vs_rest_o_proteome_PxPSSM_IUme",i,"_shift_density_ratio_gt0.5.pdf",sep=""),width=6,height=6,useDingbats=F)
 pdf(paste("bioID_vs_rest_o_proteome_LxPSSM_IUme",i,"_shift_density_ratio_gt0.5.pdf",sep=""),width=6,height=6,useDingbats=F)
 #plot(dr$x,dr$y,col="black",xlim=c(-10,50),type="l",xlab="PxIxIT score",ylab="Density",cex.axis=1.5,cex.lab=1.5,cex=1.2)
 plot(dr$x,dr$y,col="black",xlim=c(-10,50),type="l",xlab="LxVP score",ylab="Density",cex.axis=1.5,cex.lab=1.5,cex=1.2)
 lines(d5$x,d5$y,col=rgb(1,0,0))
 box(which="plot",lty=1)
 dev.off()

 print(paste("Dealing with IU >= ",i))
 # Null H0 => diff in means of the two distros = 0
 #(t.test(ransam$V2,bioid$V3,alternative="two.sided", var.equal=TRUE)) # alt H0 => diff in means ne 0
 print(t.test(bioid$V3,ransam$V2,alternative="greater", var.equal=TRUE)) # alt H0 => x > y
}

#data<-read.table("LxPSSM_topScorer_frm_REST_o_HsProteome_IUme0")
#data<-read.table("LxPSSM_topScorer_frm_REST_o_HsProteome_IUme0.1")
#data<-read.table("LxPSSM_topScorer_frm_REST_o_HsProteome_IUme0.2")
#data<-read.table("LxPSSM_topScorer_frm_REST_o_HsProteome_IUme0.3")
#h<-data$V3
#d<-density(h,bw="SJ")
#x<-hist(h,breaks=100)
#x$counts=log10(x$counts)
#summary(h)
#
##bioid<-read.table("WT_WF_ratio_highest_scoring_7mer_IUme0")
##bioid<-read.table("WT_WF_ratio_highest_scoring_7mer_IUme0.1")
##bioid<-read.table("WT_WF_ratio_highest_scoring_7mer_IUme0.2")
#bioid<-read.table("WT_WF_ratio_highest_scoring_7mer_IUme0.3")
#h5<-bioid$V12
#d5<-density(h5,bw="SJ")
#x5<-hist(h5,breaks=100)
#x5$counts=log10(x5$counts)
#summary(h5)
#
##bioid2<-read.table("WT_WF_ratio_highest_scoring_7mer_IUme0_ratio_mt1")
##bioid2<-read.table("WT_WF_ratio_highest_scoring_7mer_IUme0.1_ratio_mt1")
##bioid2<-read.table("WT_WF_ratio_highest_scoring_7mer_IUme0.2_ratio_mt1")
#bioid2<-read.table("WT_WF_ratio_highest_scoring_7mer_IUme0.3_ratio_mt1")
#hselected<-bioid2$V12
#dselected<-density(hselected,bw="SJ")
#xselected<-hist(hselected,breaks=100)
#xselected$counts=log10(xselected$counts)
#summary(hselected)
#
##bioid3<-read.table("WT_WF_ratio_highest_scoring_7mer_IUme0_w_log2_WF_log2_mt0.5")
##bioid3<-read.table("WT_WF_ratio_highest_scoring_7mer_IUme0.1_w_log2_WF_log2_mt0.5")
##bioid3<-read.table("WT_WF_ratio_highest_scoring_7mer_IUme0.2_w_log2_WF_log2_mt0.5")
#bioid3<-read.table("WT_WF_ratio_highest_scoring_7mer_IUme0.3_w_log2_WF_log2_mt0.5")
#hselected2<-bioid3$V12
#dselected2<-density(hselected2,bw="SJ")
#xselected2<-hist(hselected2,breaks=100)
#xselected2$counts=log10(xselected2$counts)
#summary(hselected2)
#
##png("bioID_vs_rest_o_proteome_LxPSSM_IUme0_shift_hist.png",width=6,height=6,units="in",res=300)
##png("bioID_vs_rest_o_proteome_LxPSSM_IUme0.1_shift_hist.png",width=6,height=6,units="in",res=300)
##png("bioID_vs_rest_o_proteome_LxPSSM_IUme0.2_shift_hist.png",width=6,height=6,units="in",res=300)
#png("bioID_vs_rest_o_proteome_LxPSSM_IUme0.3_shift_hist.png",width=6,height=6,units="in",res=300)
#plot(x,ylim=c(0,3.5),xlim=c(-30,50),xlab="LxVP score",cex.axis=1.2,cex.lab=1.2,ylab="Log10 (Frequency)")
#plot(x5,ylim=c(0,3.5),xlim=c(-30,50),xlab="",cex.axis=1.2,cex.lab=1.2,ylab="",col=rgb(1,0,0),add=T)
#plot(xselected,ylim=c(0,35),xlim=c(-30,50),xlab="",cex.axis=1.5,cex.lab=1.5,ylab="",col=rgb(1,0,1),add=T,cex=1.2)
#plot(xselected2,ylim=c(0,35),xlim=c(-30,50),xlab="",cex.axis=1.5,cex.lab=1.5,ylab="",col=rgb(0,1,0),add=T,cex=1.2)
#dev.off()
#
##png("bioID_vs_rest_o_proteome_LxPSSM_IUme0_shift_density.png",width=6,height=6,units="in",res=300)
##png("bioID_vs_rest_o_proteome_LxPSSM_IUme0.1_shift_density.png",width=6,height=6,units="in",res=300)
##png("bioID_vs_rest_o_proteome_LxPSSM_IUme0.2_shift_density.png",width=6,height=6,units="in",res=300)
#png("bioID_vs_rest_o_proteome_LxPSSM_IUme0.3_shift_density.png",width=6,height=6,units="in",res=300)
#plot(d$x,d$y,col="black",xlim=c(-30,50),ylim=c(0,0.1),type="l",xlab="LxVP score",ylab="Density",cex.axis=1.2,cex.lab=1.2)
#lines(d5$x,d5$y,col=rgb(1,0,0))
#lines(dselected$x,dselected$y,col=rgb(1,0,1))
#lines(dselected2$x,dselected2$y,col=rgb(0,0,1))
#dev.off()
#
#t.test(data$V3,bioid$V12,alternative="two.sided", var.equal=TRUE)
#t.test(data$V3,bioid2$V12,alternative="two.sided", var.equal=TRUE)
#t.test(data$V3,bioid3$V12,alternative="two.sided", var.equal=TRUE)
