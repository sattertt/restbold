
#prestats directory is first input
prestatsdir<-(commandArgs(T)[1])
setwd(prestatsdir)
#getwd()
#load component timecourses-- note that this is all hard coded
wm.name<-Sys.glob("confound_timecourses/*wm_timecourse.txt")
csf.name<-Sys.glob("confound_timecourses/*csf_timecourse.txt")
mask.name<-Sys.glob("confound_timecourses/*mask_timecourse.txt")
par<-as.matrix(read.table("mc/prefiltered_func_data_mcf.par",header=FALSE))
wm<-as.matrix(read.table(wm.name,header=FALSE))
csf<-as.matrix(read.table(csf.name,header=FALSE))
global<-as.matrix(read.table(mask.name,header=FALSE))
#combine to 9 regressors
conf<-cbind(global,csf,wm)

#get derivative and padd with 0
par.deriv<- diff(par,lag=1,differences=1)
par.pad<-matrix(0,1,6)
par.deriv.padded<-rbind(par.pad,par.deriv)

conf.deriv<- diff(conf,lag=1,differences=1)
conf.pad<-matrix(0,1,3)
conf.deriv.padded<-rbind(conf.pad,conf.deriv)


#demean both raw and derivs and square
par.dm<-sweep(par,2,colMeans(par),"-")
par.deriv.padded.dm<-sweep(par.deriv.padded,2,colMeans(par.deriv.padded),"-")
par.dm.sq<-par.dm^2
par.deriv.padded.dm.sq<-par.deriv.padded.dm^2

conf.dm<-sweep(conf,2,colMeans(conf),"-")
conf.deriv.padded.dm<-sweep(conf.deriv.padded,2,colMeans(conf.deriv.padded),"-")
conf.dm.sq<-conf.dm^2
conf.deriv.padded.dm.sq<-conf.deriv.padded.dm^2


#save final output matricies
outmat24<-cbind(par,par.deriv.padded,par.dm.sq,par.deriv.padded.dm.sq)
outmat36<-cbind(conf,conf.deriv.padded,conf.dm.sq,conf.deriv.padded.dm.sq,par,par.deriv.padded,par.dm.sq,par.deriv.padded.dm.sq)

dir.create("confound_regress_36EV")
dir.create("confound_regress_24EV")

write.table(outmat24,"confound_regress_24EV/confmat24EV.txt",row.names=FALSE,col.names=FALSE,sep="\t")
write.table(outmat36,"confound_regress_36EV/confmat36EV.txt",row.names=FALSE,col.names=FALSE,sep="\t")


