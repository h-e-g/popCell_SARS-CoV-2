options(stringsAsFactors=FALSE, max.print=9999, width=300, datatable.fread.input.cmd.message=FALSE)
.libPaths("/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/R_libs/4.1.0")

library(data.table)
library(ggplot2)
library(tictoc)

EIP='/pasteur/zeus/projets/p02/evo_immuno_pop'
PARAMS=fread(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulation_parameters.tsv",EIP))
PARAMS[,RUN:=1:.N]

result_dir=sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/",EIP)
result_files=dir(result_dir,pattern='results')

source(sprintf("%s/single_cell/resources/template_scripts/processing_pipeline/00_set_colors.R",EIP))

dir.create(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots",EIP))

##########@@ plot population size and selection
MAX_GENERATION=4962
PARAMS_SIM=data.table(generation_ago=0:MAX_GENERATION)
for(POP in c('CEU','CHS','YRI')){
  Demography=fread(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/relate/%s/spiedel_popsize_%s.coal",EIP,POP,POP))
  Ngen=round(as.numeric(names(Demography)))
  Neff=round(0.5/as.numeric(unlist(Demography)))
  DEMO=data.table(start=Ngen[2:31],end=Ngen[1:30],Neff=Neff[3:32])
  PARAMS_SIM[,POP:=sapply(generation_ago,function(x){DEMO[end<=x & x<start,Neff]})]
  setnames(PARAMS_SIM,'POP',POP)
}
PARAMS_SIM=PARAMS_SIM[1:MAX_GENERATION,]
PARAMS_SIM=melt(PARAMS_SIM,variable.name='POP',value.name='Neff',id.vars='generation_ago')

PARAMS_SEL=list()
for(ONSET in 100*(1:4)){
  for(COEF in c(0.005,0.01,0.025,0.05)){
    select_group=sprintf('s%s_+_t%s-%s',COEF,round(ONSET),round(ONSET)+20)
    PARAMS_SEL[[select_group]]=data.table(generation_ago=0:MAX_GENERATION,
                                      select = ifelse(MAX_GENERATION:0>=ONSET*10 & MAX_GENERATION:0<=(10*ONSET+200),COEF,0),
                                      onset=ONSET,
                                      select_coeff=COEF)
    }
  }
PARAMS_SEL=rbindlist(PARAMS_SEL,idcol='select_group')

color_selection=c('s0_-_t0-0'=rgb(0,0,0,0.4),
                  's0.005_+_t100-120'=rgb(1,0,0,0.4),
                  's0.01_+_t100-120'=rgb(1,0,0,0.6),
                  's0.025_+_t100-120'=rgb(1,0,0,0.8),
                  's0.05_+_t100-120'=rgb(1,0,0,1),
                  's0.005_+_t200-220'=rgb(0,1,0,0.4),
                  's0.01_+_t200-220'=rgb(0,1,0,0.6),
                  's0.025_+_t200-220'=rgb(0,1,0,0.8),
                  's0.05_+_t200-220'=rgb(0,1,0,1),
                  's0.005_+_t300-320'=rgb(0,0,1,0.4),
                  's0.01_+_t300-320'=rgb(0,0,1,0.6),
                  's0.025_+_t300-320'=rgb(0,0,1,0.8),
                  's0.05_+_t300-320'=rgb(0,0,1,1),
                  's0.005_+_t400-420'=rgb(0.7,0.7,0,0.4),
                  's0.01_+_t400-420'=rgb(0.7,0.7,0,0.6),
                  's0.025_+_t400-420'=rgb(0.7,0.7,0,0.8),
                  's0.05_+_t400-420'=rgb(0.7,0.7,0,1),
                  's0.001_+_t1-496'=rgb(0,0.7,0.7,0.4),
                  's0.0025_+_t1-496'=rgb(0,0.7,0.7,0.7),
                  's0.005_+_t1-496'=rgb(0,0.7,0.7,1),
                  's0.005_-_t100-120'=rgb(0.5,0,0,0.4),
                  's0.01_-_t100-120'=rgb(0.5,0,0,0.6),
                  's0.025_-_t100-120'=rgb(0.5,0,0,0.8),
                  's0.05_-_t100-120'=rgb(0.5,0,0,1),
                  's0.005_-_t200-220'=rgb(0,0.5,0,0.4),
                  's0.01_-_t200-220'=rgb(0,0.5,0,0.6),
                  's0.025_-_t200-220'=rgb(0,0.5,0,0.8),
                  's0.05_-_t200-220'=rgb(0,0.5,0,1),
                  's0.005_-_t300-320'=rgb(0,0,0.5,0.4),
                  's0.01_-_t300-320'=rgb(0,0,0.5,0.6),
                  's0.025_-_t300-320'=rgb(0,0,0.5,0.8),
                  's0.05_-_t300-320'=rgb(0,0,0.5,1),
                  's0.005_-_t400-420'=rgb(0.4,0.4,0,0.4),
                  's0.01_-_t400-420'=rgb(0.4,0.4,0,0.6),
                  's0.025_-_t400-420'=rgb(0.4,0.4,0,0.8),
                  's0.05_-_t400-420'=rgb(0.4,0.4,0,1),
                  's0.001_-_t1-496'=rgb(0,0.4,0.4,0.4),
                  's0.0025_-_t1-496'=rgb(0,0.4,0.4,0.7),
                  's0.005_-_t1-496'=rgb(0,0.4,0.4,1))

color_populations_1kg=setNames(color_populations,c('YRI','CEU','CHS'))

all_results_SLiM=list()
# for(file in result_files){
for (task in PARAMS[,which(COEF_SELECTION>=0 & ONSET_SELECTION!=1 & COEF_SELECTION!=0.005)]){
  file=sprintf('results_simTask%s.tsv.gz',task)
  cat(file,'\n')
  all_results_SLiM[[file]]=try(fread(sprintf('%s/%s',result_dir,file)))
  if(class(all_results_SLiM[[file]])=='try-error'){
    all_results_SLiM[[file]]=NULL
  }
}

all_results_SLiM=rbindlist(all_results_SLiM)
all_results_SLiM[,length(unique(paste(replicate,RUN)))]
notfixed=all_results_SLiM[generation_ago==0 & Freq>0 & Freq<1,paste(replicate,RUN)]
all_results_SLiM=all_results_SLiM[paste(replicate,RUN)%chin%notfixed,]
all_results_SLiM=all_results_SLiM[order(RUN,replicate,generation),]

all_results_SLiM[,DeltaFreq:=c(NA,diff(Freq)),keyby=.(RUN,replicate)]
all_results_SLiM[,DeltaFreqNorm:=DeltaFreq/sqrt(Freq*(1-Freq)),keyby=.(RUN,replicate)]
NormFactors=all_results_SLiM[COEF_SELECTION==0 & paste(replicate,RUN)%chin%notfixed,.(SD=sd(DeltaFreqNorm,na.rm=TRUE)),keyby=.(POP,generation,generation_ago)]
fwrite(NormFactors,file=sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/NormFactors.tsv.gz",EIP),sep='\t')
all_results_SLiM=merge(all_results_SLiM,NormFactors,by=c('POP','generation','generation_ago'))
all_results_SLiM[,Zscore:=DeltaFreqNorm/SD]

all_results_SLiM[,FreqFit:=pmin(1,pmax(0,predict(loess(Freq~generation,span=0.1),generation))),keyby=.(RUN,replicate)]
all_results_SLiM[,DeltaFreqFit:=c(NA,diff(FreqFit)),keyby=.(RUN,replicate)]
all_results_SLiM[,DeltaFreqFitNorm:=DeltaFreqFit/sqrt(Freq*(1-Freq)),keyby=.(RUN,replicate)]
NormFactorsFit=all_results_SLiM[COEF_SELECTION==0 & paste(replicate,RUN)%chin%notfixed,.(SDFit=sd(DeltaFreqFitNorm,na.rm=TRUE)),keyby=.(POP,generation,generation_ago)]
fwrite(NormFactorsFit,file=sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/NormFactorsFit.tsv.gz",EIP),sep='\t')
all_results_SLiM=merge(all_results_SLiM,NormFactorsFit,by=c('POP','generation','generation_ago'))
all_results_SLiM[,ZscoreFit:=DeltaFreqFitNorm/SDFit]

all_results_SLiM[,SELECTED:=generation>=ONSET_SELECTION & generation<=END_SELECTION]
all_results_SLiM[,simID:=paste(RUN,replicate)]

library(zoo)
all_results_SLiM[,FreqFit2:=c(rep(NA,4),rollmean(Freq,10),rep(NA,5)),keyby=.(RUN,replicate)]
all_results_SLiM[,DeltaFreqFit2:=c(NA,diff(FreqFit2)),keyby=.(RUN,replicate)]
all_results_SLiM[,DeltaFreqFit2Norm:=DeltaFreqFit2/sqrt(FreqFit2*(1-FreqFit2)),keyby=.(RUN,replicate)]
NormFactorsFit2=all_results_SLiM[COEF_SELECTION==0 & paste(replicate,RUN)%chin%notfixed,.(SDFit2=sd(DeltaFreqFit2Norm,na.rm=TRUE)),keyby=.(POP,generation,generation_ago)]
fwrite(NormFactorsFit2,file=sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/NormFactorsFit2.tsv.gz",EIP),sep='\t')
all_results_SLiM=merge(all_results_SLiM,NormFactorsFit2,by=c('POP','generation','generation_ago'))
all_results_SLiM[,ZscoreFit2:=DeltaFreqFit2Norm/SDFit2]

all_results_SLiM[,FreqFit3:=c(rep(NA,19),rollmean(Freq,40),rep(NA,20)),keyby=.(RUN,replicate)]
all_results_SLiM[,DeltaFreqFit3:=c(NA,diff(FreqFit3)),keyby=.(RUN,replicate)]
all_results_SLiM[,DeltaFreqFit3Norm:=DeltaFreqFit3/sqrt(FreqFit3*(1-FreqFit3)),keyby=.(RUN,replicate)]
NormFactorsFit3=all_results_SLiM[COEF_SELECTION==0 & paste(replicate,RUN)%chin%notfixed,.(SDFit3=sd(DeltaFreqFit3Norm,na.rm=TRUE)),keyby=.(POP,generation,generation_ago)]
fwrite(NormFactorsFit3,file=sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/NormFactorsFit3.tsv.gz",EIP),sep='\t')
all_results_SLiM=merge(all_results_SLiM,NormFactorsFit3,by=c('POP','generation','generation_ago'))
all_results_SLiM[,ZscoreFit3:=DeltaFreqFit3Norm/SDFit3]

all_results_SLiM[,FreqFit4:=c(rep(NA,39),rollmean(Freq,80),rep(NA,40)),keyby=.(RUN,replicate)]
all_results_SLiM[,DeltaFreqFit4:=c(NA,diff(FreqFit4)),keyby=.(RUN,replicate)]
all_results_SLiM[,DeltaFreqFit4Norm:=DeltaFreqFit4/sqrt(FreqFit4*(1-FreqFit4)),keyby=.(RUN,replicate)]
NormFactorsFit4=all_results_SLiM[COEF_SELECTION==0 & paste(replicate,RUN)%chin%notfixed,.(SDFit4=sd(DeltaFreqFit4Norm,na.rm=TRUE)),keyby=.(POP,generation,generation_ago)]
fwrite(NormFactorsFit4,file=sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/NormFactorsFit4.tsv.gz",EIP),sep='\t')
all_results_SLiM=merge(all_results_SLiM,NormFactorsFit4,by=c('POP','generation','generation_ago'))
all_results_SLiM[,ZscoreFit4:=DeltaFreqFit4Norm/SDFit4]

fwrite(all_results_SLiM,file=sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/all_results_SLiM.tsv.gz",EIP),sep='\t')

Power_th3=all_results_SLiM[,.(Zcore_power=mean(abs(Zscore)>3), ZcoreFit_power=mean(abs(ZscoreFit)>3),ZcoreFit2_power=mean(abs(ZscoreFit2)>3), ZcoreFit3_power=mean(abs(ZscoreFit3)>3), ZcoreFit4_power=mean(abs(ZscoreFit4)>3)), by=.(generation, generation_ago, COEF_SELECTION, ONSET_SELECTION, END_SELECTION, SELECTED, POP)]
Power_th3[,select_group:=paste0('s',abs(COEF_SELECTION),'_',ifelse(COEF_SELECTION>0,'+','-'),'_t',ONSET_SELECTION,'-',END_SELECTION)]
Power_th3[,TH:=3]

Power_th3.5=all_results_SLiM[,.(Zcore_power=mean(abs(Zscore)>3.5), ZcoreFit_power=mean(abs(ZscoreFit)>3.5),ZcoreFit2_power=mean(abs(ZscoreFit2)>3.5), ZcoreFit3_power=mean(abs(ZscoreFit3)>3.5), ZcoreFit4_power=mean(abs(ZscoreFit4)>3.5)), by=.(generation, generation_ago, COEF_SELECTION, ONSET_SELECTION, END_SELECTION, SELECTED, POP)]
Power_th3.5[,select_group:=paste0('s',abs(COEF_SELECTION),'_',ifelse(COEF_SELECTION>0,'+','-'),'_t',ONSET_SELECTION,'-',END_SELECTION)]
Power_th3.5[,TH:=3.5]

Power_th4=all_results_SLiM[,.(Zcore_power=mean(abs(Zscore)>4), ZcoreFit_power=mean(abs(ZscoreFit)>4),ZcoreFit2_power=mean(abs(ZscoreFit2)>4), ZcoreFit3_power=mean(abs(ZscoreFit3)>4), ZcoreFit4_power=mean(abs(ZscoreFit4)>4)), by=.(generation, generation_ago, COEF_SELECTION, ONSET_SELECTION, END_SELECTION, SELECTED, POP)]
Power_th4[,select_group:=paste0('s',abs(COEF_SELECTION),'_',ifelse(COEF_SELECTION>0,'+','-'),'_t',ONSET_SELECTION,'-',END_SELECTION)]
Power_th4[,TH:=4]

Power_allTH=rbind(Power_th3,Power_th3.5,Power_th4)
fwrite(Power_allTH,file=sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/Power_allTH.tsv.gz",EIP),sep='\t')
# Power_allTH=fread(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/Power_allTH.tsv.gz",EIP),sep='\t')

Zth=unique(c(seq(0,2,by=.1),seq(2,4,by=.01),seq(4,6,by=.1)))
typeIerror_SNPlevel=sapply(Zth,function(TH){all_results_SLiM[ONSET_SELECTION==0 & abs(Zscore)>TH, length(unique(paste(replicate,RUN)))/100000,by=POP]})
typeIerror_SNPlevel=rbindlist(apply(typeIerror_SNPlevel,2,as.data.table),idcol='TH')
typeIerror_SNPlevel[,TH:=Zth[TH]]
setnames(typeIerror_SNPlevel,'V1','typeIerror_SNPlevel')
typeIerror_SNPlevel[,MAX_GENERATION:=4962]
typeIerror_SNPlevel[,SCORE:='absZ']

typeIerror_SNPlevel2000=sapply(Zth,function(TH){all_results_SLiM[ONSET_SELECTION==0 & abs(Zscore)>TH & generation_ago<2000, length(unique(paste(replicate,RUN)))/100000,by=POP]})
typeIerror_SNPlevel2000=rbindlist(apply(typeIerror_SNPlevel2000,2,as.data.table),idcol='TH')
typeIerror_SNPlevel2000[,TH:=Zth[TH]]
setnames(typeIerror_SNPlevel2000,'V1','typeIerror_SNPlevel')
typeIerror_SNPlevel2000[,MAX_GENERATION:=2000]
typeIerror_SNPlevel2000[,SCORE:='absZ']


typeIerror_SNPlevel_Fit=sapply(Zth,function(TH){all_results_SLiM[ONSET_SELECTION==0 & abs(ZscoreFit)>TH, length(unique(paste(replicate,RUN)))/100000,by=POP]})
typeIerror_SNPlevel_Fit=rbindlist(apply(typeIerror_SNPlevel_Fit,2,as.data.table),idcol='TH')
typeIerror_SNPlevel_Fit[,TH:=Zth[TH]]
setnames(typeIerror_SNPlevel_Fit,'V1','typeIerror_SNPlevel')
typeIerror_SNPlevel_Fit[,MAX_GENERATION:=4962]
typeIerror_SNPlevel_Fit[,SCORE:='absZ_fit']

typeIerror_SNPlevel2000_Fit=sapply(Zth,function(TH){all_results_SLiM[ONSET_SELECTION==0 & abs(ZscoreFit)>TH & generation_ago<2000, length(unique(paste(replicate,RUN)))/100000,by=POP]})
typeIerror_SNPlevel2000_Fit=rbindlist(apply(typeIerror_SNPlevel2000_Fit,2,as.data.table),idcol='TH')
typeIerror_SNPlevel2000_Fit[,TH:=Zth[TH]]
setnames(typeIerror_SNPlevel2000_Fit,'V1','typeIerror_SNPlevel')
typeIerror_SNPlevel2000_Fit[,MAX_GENERATION:=2000]
typeIerror_SNPlevel2000_Fit[,SCORE:='absZ_fit']

typeIerror_SNPlevel_Fit2=sapply(Zth,function(TH){all_results_SLiM[ONSET_SELECTION==0 & abs(ZscoreFit2)>TH, length(unique(paste(replicate,RUN)))/100000,by=POP]})
typeIerror_SNPlevel_Fit2=rbindlist(apply(typeIerror_SNPlevel_Fit2,2,as.data.table),idcol='TH')
typeIerror_SNPlevel_Fit2[,TH:=Zth[TH]]
setnames(typeIerror_SNPlevel_Fit2,'V1','typeIerror_SNPlevel')
typeIerror_SNPlevel_Fit2[,MAX_GENERATION:=4962]
typeIerror_SNPlevel_Fit2[,SCORE:='absZ_fitMA']

typeIerror_SNPlevel2000_Fit2=sapply(Zth,function(TH){all_results_SLiM[ONSET_SELECTION==0 & abs(ZscoreFit2)>TH & generation_ago<2000, length(unique(paste(replicate,RUN)))/100000,by=POP]})
typeIerror_SNPlevel2000_Fit2=rbindlist(apply(typeIerror_SNPlevel2000_Fit2,2,as.data.table),idcol='TH')
typeIerror_SNPlevel2000_Fit2[,TH:=Zth[TH]]
setnames(typeIerror_SNPlevel2000_Fit2,'V1','typeIerror_SNPlevel')
typeIerror_SNPlevel2000_Fit2[,MAX_GENERATION:=2000]
typeIerror_SNPlevel2000_Fit2[,SCORE:='absZ_fitMA']

typeIerror_SNPlevel_all=rbind(typeIerror_SNPlevel,
                            typeIerror_SNPlevel2000,
                            typeIerror_SNPlevel_Fit,
                            typeIerror_SNPlevel2000_Fit,
                            typeIerror_SNPlevel_Fit2,
                            typeIerror_SNPlevel2000_Fit2)
fwrite(typeIerror_SNPlevel_all,file=sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/typeIerror_SNPlevel.tsv.gz",EIP),sep='\t')


DIFF_onset=list()
for (TH in c(3,3.5,4)){
  for (qTH in c(0.8,.9,.95,.99,1)){
    DIFF_onset[[paste(TH,qTH,'all')]]=all_results_SLiM[abs(ZscoreFit)>TH,.(DIFF=quantile(generation,1-qTH)-ONSET_SELECTION),by=.(replicate,POP,COEF_SELECTION,ONSET_SELECTION)][,.(DIFFavg=mean(abs(DIFF),na.rm=T),DIFFunder5=mean(abs(DIFF)<5,na.rm=T),DIFFunder10=mean(abs(DIFF)<10,na.rm=T)),keyby=.(ONSET_SELECTION,abs(COEF_SELECTION))]
    DIFF_onset[[paste(TH,qTH,'all')]][,Threshold:=TH]
    DIFF_onset[[paste(TH,qTH,'all')]][,PctIgnored:=(1-qTH)*100]
    DIFF_onset[[paste(TH,qTH,'2000gen')]]=all_results_SLiM[ONSET_SELECTION==100 & abs(ZscoreFit)>TH & generation_ago<2000,.(DIFF=quantile(generation,1-qTH)-ONSET_SELECTION),by=.(replicate,POP,COEF_SELECTION,ONSET_SELECTION)][,.(DIFFavg=mean(abs(DIFF),na.rm=T),DIFFunder5=mean(abs(DIFF)<5,na.rm=T),DIFFunder10=mean(abs(DIFF)<10,na.rm=T)),keyby=.(ONSET_SELECTION,abs(COEF_SELECTION))]
    DIFF_onset[[paste(TH,qTH,'2000gen')]][,Threshold:=TH]
    DIFF_onset[[paste(TH,qTH,'2000gen')]][,PctIgnored:=(1-qTH)*100]
  }
}
DIFF_onset=rbindlist(DIFF_onset)
fwrite(DIFF_onset,file=sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/DIFF_onset_selection_TRUE_vs_estimated.tsv.gz",EIP),sep='\t')


TH_dt=data.table(expand.grid(simID=all_results_SLiM[generation==110,simID], TH=Zth))
Power_thChoice_110gen=merge(all_results_SLiM[generation==110,],TH_dt,by='simID')
Power_thChoice_110gen=Power_thChoice_110gen[,.(Zcore_power=mean(abs(Zscore)>TH), ZcoreFit_power=mean(abs(ZscoreFit)>TH),ZcoreFit2_power=mean(abs(ZscoreFit2)>TH), ZcoreFit3_power=mean(abs(ZscoreFit3)>TH),ZcoreFit4_power=mean(abs(ZscoreFit4)>TH)), by=.(TH, generation, generation_ago, COEF_SELECTION, ONSET_SELECTION, END_SELECTION, SELECTED, POP)]
Power_thChoice_110gen[,select_group:=paste0('s',abs(COEF_SELECTION),'_',ifelse(COEF_SELECTION>0,'+','-'),'_t',ONSET_SELECTION,'-',END_SELECTION)]

TH_dt=data.table(expand.grid(simID=all_results_SLiM[generation==210,simID], TH=Zth))
Power_thChoice_210gen=merge(all_results_SLiM[generation==210,],TH_dt,by='simID')
Power_thChoice_210gen=Power_thChoice_210gen[,.(Zcore_power=mean(abs(Zscore)>TH), ZcoreFit_power=mean(abs(ZscoreFit)>TH),ZcoreFit2_power=mean(abs(ZscoreFit2)>TH), ZcoreFit3_power=mean(abs(ZscoreFit3)>TH),ZcoreFit4_power=mean(abs(ZscoreFit4)>TH)), by=.(TH, generation, generation_ago, COEF_SELECTION, ONSET_SELECTION, END_SELECTION, SELECTED, POP)]
Power_thChoice_210gen[,select_group:=paste0('s',abs(COEF_SELECTION),'_',ifelse(COEF_SELECTION>0,'+','-'),'_t',ONSET_SELECTION,'-',END_SELECTION)]

TH_dt=data.table(expand.grid(simID=all_results_SLiM[generation==310,simID], TH=Zth))
Power_thChoice_310gen=merge(all_results_SLiM[generation==310,],TH_dt,by='simID')
Power_thChoice_310gen=Power_thChoice_310gen[,.(Zcore_power=mean(abs(Zscore)>TH), ZcoreFit_power=mean(abs(ZscoreFit)>TH),ZcoreFit2_power=mean(abs(ZscoreFit2)>TH), ZcoreFit3_power=mean(abs(ZscoreFit3)>TH),ZcoreFit4_power=mean(abs(ZscoreFit4)>TH)), by=.(TH, generation, generation_ago, COEF_SELECTION, ONSET_SELECTION, END_SELECTION, SELECTED, POP)]
Power_thChoice_310gen[,select_group:=paste0('s',abs(COEF_SELECTION),'_',ifelse(COEF_SELECTION>0,'+','-'),'_t',ONSET_SELECTION,'-',END_SELECTION)]

TH_dt=data.table(expand.grid(simID=all_results_SLiM[generation==410,simID], TH=Zth))
Power_thChoice_410gen=merge(all_results_SLiM[generation==410,],TH_dt,by='simID')
Power_thChoice_410gen=Power_thChoice_410gen[,.(Zcore_power=mean(abs(Zscore)>TH), ZcoreFit_power=mean(abs(ZscoreFit)>TH),ZcoreFit2_power=mean(abs(ZscoreFit2)>TH), ZcoreFit3_power=mean(abs(ZscoreFit3)>TH),ZcoreFit4_power=mean(abs(ZscoreFit4)>TH)), by=.(TH, generation, generation_ago, COEF_SELECTION, ONSET_SELECTION, END_SELECTION, SELECTED, POP)]
Power_thChoice_410gen[,select_group:=paste0('s',abs(COEF_SELECTION),'_',ifelse(COEF_SELECTION>0,'+','-'),'_t',ONSET_SELECTION,'-',END_SELECTION)]

Power_thChoice=rbind(Power_thChoice_110gen[ONSET_SELECTION%in%c(0,100)& COEF_SELECTION>=0,],
                     Power_thChoice_210gen[ONSET_SELECTION%in%c(0,200)& COEF_SELECTION>=0,],
                     Power_thChoice_310gen[ONSET_SELECTION%in%c(0,300)& COEF_SELECTION>=0,],
                     Power_thChoice_410gen[ONSET_SELECTION%in%c(0,400)& COEF_SELECTION>=0,])
fwrite(Power_thChoice,file=sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/Power_thChoice.tsv.gz",EIP),sep='\t')
Power_thChoice=fread(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/Power_thChoice.tsv.gz",EIP))

Power_thChoice2=merge(Power_thChoice[ONSET_SELECTION==0,.(TH,generation,generation_ago,POP,ZcoreFit_power,ZcoreFit2_power,ZcoreFit3_power)],Power_thChoice[ONSET_SELECTION>0,],by=c('generation','POP','TH','generation_ago'),suffix=c('.neutral',''))



#check difference in onset :
#all_results_SLiM[ONSET_SELECTION>2,.(DIFF=min(generation[ZscoreFit>3])-ONSET_SELECTION),by=.(replicate,POP,COEF_SELECTION,ONSET_SELECTION)][,mean(abs(DIFF)<50,na.rm=T),keyby=.(ONSET_SELECTION,abs(COEF_SELECTION))]

# all_results_SLiM[,COEF_SELECTION]

p <- ggplot(all_results_SLiM[generation==1,][,.SD[sample(1:.N,min(10000,.N)),],by=.(POP,COEF_SELECTION)],aes(x=pmin(1,pmax(0,Freq+runif(length(Freq),0,.025)-0.0125))))+geom_histogram(breaks=seq(-.001,1.001,l=35))+facet_grid(COEF_SELECTION~POP)+theme_yann()
pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/histogram_freq_at_gen0.pdf",EIP))
print(p)
dev.off()

p <- ggplot(all_results_SLiM[generation_ago==0,][,.SD[sample(1:.N,min(10000,.N)),],by=.(POP,COEF_SELECTION)],aes(x=pmin(1,pmax(0,Freq+runif(length(Freq),0,.025)-0.0125))))+geom_histogram(breaks=seq(-.001,1.001,l=40))+facet_grid(COEF_SELECTION~POP)+theme_yann()
pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/histogram_freq_Now.pdf",EIP))
print(p)
dev.off()

# p <- ggplot(NormFactors,aes(x=log10(generation_ago*10),y=log10(10/SD^2),col=POP))+geom_line()+geom_point()+theme_yann()+scale_color_manual(values=color_populations_1kg)
# p <- p+ylab('log10(Ne)') +xlab('log10(generation ago)')
# pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/NormFactors_overTime.pdf",EIP))
# print(p)
# dev.off()

p <- ggplot(NormFactorsFit,aes(x=log10(generation_ago*10),y=log10(10/SDFit^2),col=POP))+geom_line()+geom_point()+theme_yann()+scale_color_manual(values=color_populations_1kg)
p <- p+ylab('log10(Ne)') +xlab('log10(generation ago)')
pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/NormFactorsFit_overTime.pdf",EIP))
print(p)
dev.off()
#
# p <- ggplot(NormFactorsFit2,aes(x=log10(generation_ago*10),y=log10(10/SDFit2^2),col=POP))+geom_line()+geom_point()+theme_yann()+scale_color_manual(values=color_populations_1kg)
# p <- p+ylab('log10(Ne)') +xlab('log10(generation ago)')
# pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/NormFactorsFit2_overTime.pdf",EIP))
# print(p)
# dev.off()
#
# p <- ggplot(NormFactorsFit3,aes(x=log10(generation_ago*10),y=log10(10/SDFit3^2),col=POP))+geom_line()+geom_point()+theme_yann()+scale_color_manual(values=color_populations_1kg)
# p <- p+ylab('log10(Ne)') +xlab('log10(generation ago)')
# pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/NormFactorsFit3_overTime.pdf",EIP))
# print(p)
# dev.off()


p <- ggplot()+geom_line(data=PARAMS_SIM,mapping=aes(x=generation_ago,y=log10(Neff),col=POP))+theme_yann()+scale_color_manual(values=color_populations_1kg)+
geom_rect(data=data.table(xmin=4960-1200,xmax=4960-1000,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(1,0,0),alpha=0.2)+
geom_rect(data=data.table(xmin=4960-2200,xmax=4960-2000,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,1,0),alpha=0.2)+
geom_rect(data=data.table(xmin=4960-3200,xmax=4960-3000,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,0,1),alpha=0.2)+
geom_rect(data=data.table(xmin=4960-4200,xmax=4960-4000,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0.7,0.7,0),alpha=0.2)
p <- p+ylab('log10(Ne)') +xlab('generation ago')+coord_cartesian(ylim=c(0,8))
pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Neff_overTime.pdf",EIP),width=4,height=4)
print(p)
dev.off()

p <- ggplot(data=PARAMS_SIM,mapping=aes(x=log10(generation_ago),y=log10(Neff),col=POP))+geom_line()+theme_yann()+scale_color_manual(values=color_populations_1kg)
p <- p+ylab('log10(Ne)') +xlab('log10(generation ago)')
pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Neff_overTime_log.pdf",EIP))
print(p)
dev.off()

p <- ggplot()+geom_line(data=PARAMS_SEL,mapping=aes(x=generation_ago,y=select,col=select_group))+theme_yann()+scale_color_manual(values=color_selection[names(color_selection)%in%PARAMS_SEL$select_group])
p <- p+ylab('selection coefficient') +xlab('generation ago')+coord_cartesian(ylim=c(0,.1))+
geom_rect(data=data.table(xmin=4960-1200,xmax=4960-1000,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(1,0,0),alpha=0.2)+
geom_rect(data=data.table(xmin=4960-2200,xmax=4960-2000,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,1,0),alpha=0.2)+
geom_rect(data=data.table(xmin=4960-3200,xmax=4960-3000,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,0,1),alpha=0.2)+
geom_rect(data=data.table(xmin=4960-4200,xmax=4960-4000,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0.7,0.7,0),alpha=0.2)
pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Selection_overTime.pdf",EIP),height=4,width=4)
print(p)
dev.off()

# Power_vs_threshold_ZscoreFit
p <- ggplot()+
  geom_line(data=Power_thChoice, mapping=aes(TH, ZcoreFit_power, group=select_group, color=select_group))+
  scale_color_manual(values=color_selection[names(color_selection)%in%Power_thChoice$select_group])+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
  xlab("Threshold")+ylab("Power to detect selection (|Z-score|>Threshold)")+
  coord_cartesian(ylim=c(0,1))+
  facet_grid(rows=vars(POP),col=vars(generation))+
  theme(panel.spacing=unit(0,'pt'))+theme_yann()

pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Power_vs_threshold_ZscoreFit_all4timePoints.pdf",EIP))
print(p)
dev.off()

# Power_vs_threshold_ZscoreFit CHS 410 gen
p <- ggplot()+
  geom_line(data=Power_thChoice_410gen[POP=='CHS',], mapping=aes(TH, ZcoreFit_power, group=select_group, color=select_group,linetype=as.character(COEF_SELECTION)))+scale_linetype_manual(values=c('0'=1,'0.05'=1,'0.025'=2,'0.01'=3))+
  scale_color_manual(values=color_selection[names(color_selection)%in%Power_thChoice$select_group])+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
  xlab("Threshold")+ylab("Power to detect selection (|Z-score|>Threshold)")+
  coord_cartesian(ylim=c(0,1))+
  facet_grid(rows=vars(POP),col=vars(generation))+
  theme(panel.spacing=unit(0,'pt'))+theme_yann()

pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Power_vs_threshold_ZscoreFit_410gen_CHS_v2.pdf",EIP),height=5,width=4)
print(p)
dev.off()


# Power_vs_threshold_ZscoreFit
# p <- ggplot()+
#   geom_line(data=Power_thChoice2, mapping=aes(x=ZcoreFit_power.neutral, ZcoreFit_power, group=select_group, color=select_group))+
#   scale_color_manual(values=color_selection[names(color_selection)%in%c(Power_thChoice2$select_group,'s0_-_t0-0')])+
#   scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
#   xlab("specificity")+ylab("Power to detect selection (|Z-score|>Threshold)")+
#   coord_cartesian(ylim=c(0,1))+
#   facet_grid(rows=vars(POP),col=vars(generation))+
#   theme(panel.spacing=unit(0,'pt'))+theme_yann()
#
# pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/ROC_Power_vs_alpha_ZscoreFit_all4timePoints.pdf",EIP))
# print(p)
# dev.off()


# Power_vs_threshold_Zscore_110_gen
# p <- ggplot()+
#   geom_line(data=Power_thChoice_110gen[ONSET_SELECTION%in%c(0,100),],mapping=aes(TH,Zcore_power,group=select_group,color=select_group))+
#   scale_color_manual(values=color_selection)+
#   scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
#   xlab("Threshold")+ylab("Power to detect selection (|Z-score|>TH) (110 generation)")+
# coord_cartesian(ylim=c(0,1))+
#   facet_grid(rows=vars(POP),col=vars(abs(COEF_SELECTION)))+
#   theme(panel.spacing=unit(0,'pt'))+theme_yann()
#
#   pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Power_vs_threshold_Zscore_110_gen.pdf",EIP))
#   print(p)
#   dev.off()
#
#   p <- ggplot()+
#     geom_line(data=Power_thChoice_110gen[ONSET_SELECTION%in%c(0,100),],mapping=aes(TH,ZcoreFit_power,group=select_group,color=select_group))+
#     scale_color_manual(values=color_selection)+
#     scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
#     xlab("Threshold")+ylab("Power to detect selection (|Z-score|>TH) (110 generation)")+
#   coord_cartesian(ylim=c(0,1))+
#     facet_grid(rows=vars(POP),col=vars(abs(COEF_SELECTION)))+
#     theme(panel.spacing=unit(0,'pt'))+theme_yann()
#
#     pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Power_vs_threshold_ZscoreFit_110_gen.pdf",EIP))
#     print(p)
#     dev.off()
#
#
# averaged_signal=all_results_SLiM[,.(Zcore_avg=mean(Zscore),ZcoreFit_avg=mean(ZscoreFit),ZcoreFit2_avg=mean(ZscoreFit2)),by=.(reps_group=(replicate-1)%/%100,generation,generation_ago,RUN,FREQ0,COEF_SELECTION,ONSET_SELECTION,END_SELECTION,SELECTED,POP)]
# averaged_signal[,select_group:=paste0('s',abs(COEF_SELECTION),'_',ifelse(COEF_SELECTION>0,'+','-'),'_t',ONSET_SELECTION,'-',END_SELECTION)]
#


#### Power Zcore th3
# p <- ggplot()+
# geom_rect(data=data.table(xmin=496-120,xmax=496-100,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(1,0,0),alpha=0.2)+
# geom_rect(data=data.table(xmin=496-220,xmax=496-200,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,1,0),alpha=0.2)+
# geom_rect(data=data.table(xmin=496-320,xmax=496-300,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,0,1),alpha=0.2)+
# geom_rect(data=data.table(xmin=496-420,xmax=496-400,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0.7,0.7,0),alpha=0.2)+
#   geom_line(data=Power_th3[ONSET_SELECTION!=1,],mapping=aes(generation_ago,Zcore_power,group=select_group,color=select_group))+
#   scale_color_manual(values=color_selection)+
#   scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
#   xlab("Generations from present")+ylab("Power to detect selection (|Z-score|>3)")+
# coord_cartesian(ylim=c(0,1))+
#   facet_grid(rows=vars(POP),cols=vars(COEF_SELECTION))+
#   theme(panel.spacing=unit(0,'pt'))+theme_yann()
#
#   pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Power_overTime_Zscore_th3.pdf",EIP))
#   print(p)
#   dev.off()
#
#   #### Power Zcore th3
#   p <- ggplot()+
#     geom_line(data=Power_th3[ONSET_SELECTION<2,],mapping=aes(generation_ago,Zcore_power,group=select_group,color=select_group))+
#     scale_color_manual(values=color_selection)+
#     scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
#     xlab("Generations from present")+ylab("Power to detect selection (|Z-score|>3)")+
#   coord_cartesian(ylim=c(0,1))+
#     facet_grid(rows=vars(POP),cols=vars(COEF_SELECTION))+
#     theme(panel.spacing=unit(0,'pt'))+theme_yann()
#
#     pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Power_overTime_Zscore_th3_continuous.pdf",EIP))
#     print(p)
#     dev.off()

#### Power Zcore Fit th3
p <- ggplot()+
geom_rect(data=data.table(xmin=496-120,xmax=496-100,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(1,0,0),alpha=0.2)+
geom_rect(data=data.table(xmin=496-220,xmax=496-200,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,1,0),alpha=0.2)+
geom_rect(data=data.table(xmin=496-320,xmax=496-300,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,0,1),alpha=0.2)+
geom_rect(data=data.table(xmin=496-420,xmax=496-400,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0.7,0.7,0),alpha=0.2)+
  geom_line(data=Power_th3[ONSET_SELECTION!=1 & (COEF_SELECTION==0 | COEF_SELECTION>=.01),],mapping=aes(generation_ago,ZcoreFit_power,group=select_group,color=select_group))+
  scale_color_manual(values=color_selection)+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
  xlab("Generations from present")+ylab("Power to detect selection (|Z-score|>3)")+
coord_cartesian(ylim=c(0,1))+
  facet_grid(rows=vars(POP),cols=vars(COEF_SELECTION))+
  theme(panel.spacing=unit(0,'pt'))+theme_yann()

  pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Power_overTime_ZscoreFit_th3_v2.pdf",EIP))
  print(p)
  dev.off()



  #### Power Zcore Fit th3
  p <- ggplot()+
  geom_rect(data=data.table(xmin=496-120,xmax=496-100,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(1,0,0),alpha=0.2)+
  geom_rect(data=data.table(xmin=496-220,xmax=496-200,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,1,0),alpha=0.2)+
  geom_rect(data=data.table(xmin=496-320,xmax=496-300,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,0,1),alpha=0.2)+
  geom_rect(data=data.table(xmin=496-420,xmax=496-400,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0.7,0.7,0),alpha=0.2)+
    geom_line(data=Power_th3[ONSET_SELECTION!=1 & (COEF_SELECTION==0 | COEF_SELECTION>=.01),],mapping=aes(generation_ago,ZcoreFit_power,group=select_group,color=select_group))+
    scale_color_manual(values=color_selection)+
    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
    xlab("Generations from present")+ylab("Power to detect selection (|Z-score|>3)")+
  coord_cartesian(ylim=c(0,1))+
    facet_grid(rows=vars(POP),cols=vars(COEF_SELECTION))+
    theme(panel.spacing=unit(0,'pt'))+theme_yann()

    pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Power_overTime_ZscoreFit_th3.pdf",EIP))
    print(p)
    dev.off()


#### Power Zcore Fit th3 continuous
p <- ggplot()+
  geom_line(data=Power_th3[ONSET_SELECTION<2,],mapping=aes(generation_ago,ZcoreFit_power,group=select_group,color=select_group))+
  scale_color_manual(values=color_selection)+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
  xlab("Generations from present")+ylab("Power to detect selection (|Z-score|>3)")+
coord_cartesian(ylim=c(0,1))+
  facet_grid(rows=vars(POP),cols=vars(COEF_SELECTION))+
  theme(panel.spacing=unit(0,'pt'))+theme_yann()

  pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Power_overTime_ZscoreFit_th3_continuous.pdf",EIP))
  print(p)
  dev.off()


#### Power Zcore Fit2 th3
p <- ggplot()+
geom_rect(data=data.table(xmin=496-120,xmax=496-100,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(1,0,0),alpha=0.2)+
geom_rect(data=data.table(xmin=496-220,xmax=496-200,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,1,0),alpha=0.2)+
geom_rect(data=data.table(xmin=496-320,xmax=496-300,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,0,1),alpha=0.2)+
geom_rect(data=data.table(xmin=496-420,xmax=496-400,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0.7,0.7,0),alpha=0.2)+
  geom_line(data=Power_th3[ONSET_SELECTION!=1,],mapping=aes(generation_ago,ZcoreFit2_power,group=select_group,color=select_group))+
  scale_color_manual(values=color_selection)+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
  xlab("Generations from present")+ylab("Power to detect selection (|Z-score|>3)")+
coord_cartesian(ylim=c(0,1))+
  facet_grid(rows=vars(POP),cols=vars(COEF_SELECTION))+
  theme(panel.spacing=unit(0,'pt'))+theme_yann()

  pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Power_overTime_ZscoreFit2_th3.pdf",EIP))
  print(p)
  dev.off()


#### Power Zcore Fit2 th3 continous
p <- ggplot()+
  geom_line(data=Power_th3[ONSET_SELECTION<2,],mapping=aes(generation_ago,ZcoreFit2_power,group=select_group,color=select_group))+
  scale_color_manual(values=color_selection)+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
  xlab("Generations from present")+ylab("Power to detect selection (|Z-score|>3)")+
coord_cartesian(ylim=c(0,1))+
  facet_grid(rows=vars(POP),cols=vars(COEF_SELECTION))+
  theme(panel.spacing=unit(0,'pt'))+theme_yann()

  pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Power_overTime_ZscoreFit2_th3_continuous.pdf",EIP))
  print(p)
  dev.off()


#### Power Zcore Fit3 th3
p <- ggplot()+
geom_rect(data=data.table(xmin=496-120,xmax=496-100,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(1,0,0),alpha=0.2)+
geom_rect(data=data.table(xmin=496-220,xmax=496-200,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,1,0),alpha=0.2)+
geom_rect(data=data.table(xmin=496-320,xmax=496-300,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,0,1),alpha=0.2)+
geom_rect(data=data.table(xmin=496-420,xmax=496-400,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0.7,0.7,0),alpha=0.2)+
  geom_line(data=Power_th3[ONSET_SELECTION!=1,],mapping=aes(generation_ago,ZcoreFit2_power,group=select_group,color=select_group))+
  scale_color_manual(values=color_selection)+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
  xlab("Generations from present")+ylab("Power to detect selection (|Z-score|>3)")+
coord_cartesian(ylim=c(0,1))+
  facet_grid(rows=vars(POP),cols=vars(COEF_SELECTION))+
  theme(panel.spacing=unit(0,'pt'))+theme_yann()

  pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Power_overTime_ZscoreFit3_th3.pdf",EIP))
  print(p)
  dev.off()


  #### Power Zcore Fit3 th3 CONTINOUS
  p <- ggplot()+
    geom_line(data=Power_th3[ONSET_SELECTION<2,],mapping=aes(generation_ago,ZcoreFit3_power,group=select_group,color=select_group))+
    scale_color_manual(values=color_selection)+
    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
    xlab("Generations from present")+ylab("Power to detect selection (|Z-score|>3)")+
  coord_cartesian(ylim=c(0,1))+
    facet_grid(rows=vars(POP),cols=vars(COEF_SELECTION))+
    theme(panel.spacing=unit(0,'pt'))+theme_yann()

    pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Power_overTime_ZscoreFit3_th3_continuous.pdf",EIP))
    print(p)
    dev.off()


#### Power Zcore Fit4 th3
p <- ggplot()+
  geom_line(data=Power_th3[ONSET_SELECTION<2,],mapping=aes(generation_ago,ZcoreFit4_power,group=select_group,color=select_group))+
  scale_color_manual(values=color_selection)+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
  xlab("Generations from present")+ylab("Power to detect selection (|Z-score|>3)")+
coord_cartesian(ylim=c(0,1))+
  facet_grid(rows=vars(POP),cols=vars(COEF_SELECTION))+
  theme(panel.spacing=unit(0,'pt'))+theme_yann()

  pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Power_overTime_ZscoreFit4_th3_continuous.pdf",EIP))
  print(p)
  dev.off()



  #### Power Zcore Fit th3.5
  p <- ggplot()+
  geom_rect(data=data.table(xmin=496-120,xmax=496-100,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(1,0,0),alpha=0.2)+
  geom_rect(data=data.table(xmin=496-220,xmax=496-200,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,1,0),alpha=0.2)+
  geom_rect(data=data.table(xmin=496-320,xmax=496-300,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,0,1),alpha=0.2)+
  geom_rect(data=data.table(xmin=496-420,xmax=496-400,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0.7,0.7,0),alpha=0.2)+
    geom_line(data=Power_th3.5[ONSET_SELECTION!=1,],mapping=aes(generation_ago,ZcoreFit_power,group=select_group,color=select_group))+
    scale_color_manual(values=color_selection)+
    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
    xlab("Generations from present")+ylab("Power to detect selection (|Z-score|>3)")+
  coord_cartesian(ylim=c(0,1))+
    facet_grid(rows=vars(POP),cols=vars(COEF_SELECTION))+
    theme(panel.spacing=unit(0,'pt'))+theme_yann()

    pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Power_overTime_ZscoreFit_th3.5.pdf",EIP))
    print(p)
    dev.off()

    #### Power Zcore Fit th4
    p <- ggplot()+
    geom_rect(data=data.table(xmin=496-120,xmax=496-100,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(1,0,0),alpha=0.2)+
    geom_rect(data=data.table(xmin=496-220,xmax=496-200,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,1,0),alpha=0.2)+
    geom_rect(data=data.table(xmin=496-320,xmax=496-300,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,0,1),alpha=0.2)+
    geom_rect(data=data.table(xmin=496-420,xmax=496-400,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0.7,0.7,0),alpha=0.2)+
      geom_line(data=Power_th4[ONSET_SELECTION!=1,],mapping=aes(generation_ago,ZcoreFit_power,group=select_group,color=select_group))+
      scale_color_manual(values=color_selection)+
      scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
      xlab("Generations from present")+ylab("Power to detect selection (|Z-score|>3)")+
    coord_cartesian(ylim=c(0,1))+
      facet_grid(rows=vars(POP),cols=vars(COEF_SELECTION))+
      theme(panel.spacing=unit(0,'pt'))+theme_yann()

      pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Power_overTime_ZscoreFit_th4.pdf",EIP))
      print(p)
      dev.off()




#
# p <- ggplot()+
# geom_rect(data=data.table(xmin=496-120,xmax=496-100,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(1,0,0),alpha=0.2)+
# geom_rect(data=data.table(xmin=496-220,xmax=496-200,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,1,0),alpha=0.2)+
# geom_rect(data=data.table(xmin=496-320,xmax=496-300,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,0,1),alpha=0.2)+
#   geom_line(data=Power_th3,mapping=aes(generation_ago,Zcore_avg,group=paste(RUN,reps_group),color=select_group))+
#   scale_color_manual(values=color_selection)+
# #  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
#   xlab("Generations from present")+ylab("mean Z-score (100 replicates)")+
# coord_cartesian(ylim=c(-6,6))+
#   facet_grid(rows=vars(POP))+
#   theme(panel.spacing=unit(0,'pt'))
  #
  # pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Frequency_overTime_avg_100replicates_v2.pdf",EIP))
  # print(p)
  # dev.off()
  #
  #
  # p <- ggplot()+
  # geom_rect(data=data.table(xmin=496-120,xmax=496-100,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(1,0,0),alpha=0.2)+
  # geom_rect(data=data.table(xmin=496-220,xmax=496-200,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,1,0),alpha=0.2)+
  # geom_rect(data=data.table(xmin=496-320,xmax=496-300,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,0,1),alpha=0.2)+
  #   geom_line(data=averaged_signal,mapping=aes(generation_ago,ZcoreFit_avg,group=select_group,color=select_group))+
  #   scale_color_manual(values=color_selection)+
  # #  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
  #   xlab("Generations from present")+ylab("mean Z-score (100 replicates)")+
  # coord_cartesian(ylim=c(-6,6))+
  #   facet_grid(rows=vars(POP))+
  #   theme(panel.spacing=unit(0,'pt'))
  #
  #   pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Frequency_overTime_avg_100replicates_Fit_v2.pdf",EIP))
  #   print(p)
  #   dev.off()



      averaged_signal_all=all_results_SLiM[,.(Zcore_avg=mean(Zscore),ZcoreFit_avg=mean(ZscoreFit),ZcoreFit2_avg=mean(ZscoreFit2)),by=.(reps_group=(replicate-1)%/%100,generation,generation_ago,COEF_SELECTION,ONSET_SELECTION,END_SELECTION,SELECTED,POP)]
      averaged_signal_all[,select_group:=paste0('s',abs(COEF_SELECTION),'_',ifelse(COEF_SELECTION>0,'+','-'),'_t',ONSET_SELECTION,'-',END_SELECTION)]


p <- ggplot()+
geom_rect(data=data.table(xmin=496-120,xmax=496-100,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(1,0,0),alpha=0.2)+
geom_rect(data=data.table(xmin=496-220,xmax=496-200,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,1,0),alpha=0.2)+
geom_rect(data=data.table(xmin=496-320,xmax=496-300,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,0,1),alpha=0.2)+
geom_rect(data=data.table(xmin=496-420,xmax=496-400,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0.7,0.7,0),alpha=0.2)+
  geom_line(data=averaged_signal_all,mapping=aes(generation_ago,Zcore_avg,group=paste(reps_group,select_group),color=select_group))+
  geom_hline(aes(yintercept=c(-3,3)),col='lightgrey',type=2)+
  scale_color_manual(values=color_selection)+
#  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
  xlab("Generations from present")+ylab("mean Z-score (100 replicates)")+
coord_cartesian(ylim=c(-10,10))+
  facet_grid(rows=vars(POP))+theme_yann()+
  theme(panel.spacing=unit(0,'pt'))

  pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Frequency_overTime_avg_100replicates_v3.pdf",EIP))
  print(p)
dev.off()

      p <- ggplot()+
      geom_rect(data=data.table(xmin=496-120,xmax=496-100,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(1,0,0),alpha=0.2)+
      geom_rect(data=data.table(xmin=496-220,xmax=496-200,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,1,0),alpha=0.2)+
      geom_rect(data=data.table(xmin=496-320,xmax=496-300,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,0,1),alpha=0.2)+
      geom_rect(data=data.table(xmin=496-420,xmax=496-400,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0.7,0.7,0),alpha=0.2)+
        geom_line(data=averaged_signal_all,mapping=aes(generation_ago,ZcoreFit_avg,group=paste(reps_group,select_group),color=select_group))+
        geom_hline(aes(yintercept=c(-3,3)),col='lightgrey',type=2)+
        scale_color_manual(values=color_selection)+
      #  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
        xlab("Generations from present")+ylab("mean Z-score (100 replicates)")+
      coord_cartesian(ylim=c(-10,10))+
        facet_grid(rows=vars(POP))+theme_yann()+
        theme(panel.spacing=unit(0,'pt'))

        pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Frequency_overTime_avg_100replicates_Fit_v3.pdf",EIP))
        print(p)
        dev.off()

        p <- ggplot()+
        geom_rect(data=data.table(xmin=496-120,xmax=496-100,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(1,0,0),alpha=0.2)+
        geom_rect(data=data.table(xmin=496-220,xmax=496-200,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,1,0),alpha=0.2)+
        geom_rect(data=data.table(xmin=496-320,xmax=496-300,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0,0,1),alpha=0.2)+
        geom_rect(data=data.table(xmin=496-420,xmax=496-400,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=rgb(0.7,0.7,0),alpha=0.2)+
          geom_line(data=averaged_signal_all,mapping=aes(generation_ago,ZcoreFit2_avg,group=paste(reps_group,select_group),color=select_group))+
          geom_hline(aes(yintercept=c(-3,3)),col='lightgrey',type=2)+
          scale_color_manual(values=color_selection)+
        #  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0.0","0.25","0.5","0.75","1.0"))+
          xlab("Generations from present")+ylab("mean Z-score (100 replicates)")+
        coord_cartesian(ylim=c(-10,10))+
          facet_grid(rows=vars(POP))+theme_yann()+
          theme(panel.spacing=unit(0,'pt'))

          pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/Frequency_overTime_avg_100replicates_Fit2_v3.pdf",EIP))
          print(p)
          dev.off()


figs8d_data=clues[pop=="CHS"&rsID%in%c(rand_snps,int_snps),]
figs8d_data[,pop:=factor(pop,rev(c('YRI','CEU','CHS')))]
figs8d_data[,random:=ifelse(rsID%chin%int_snps,F,T)]

figs8d_plot=ggplot()+
geom_rect(data=data.table(xmin=770,xmax=970,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="#71458d",alpha=0.2)+
  geom_line(data=figs8d_data[random==T,],mapping=aes(epoch,allele_pp_max_smooth,group=rsID,color=pop),alpha=0.2,color="gray")+
  geom_line(data=figs8d_data[random==F,],mapping=aes(epoch,allele_pp_max_smooth,group=rsID,color=pop))+
  scale_color_manual(values=color_populations_1kg)+
  scale_y_continuous(breaks=c(0,0.5,1),labels=c("0.0","0.5","1.0"))+
  xlab("Generations from present")+ylab("Frequency (f)")+
coord_cartesian(ylim=c(0,1))+
  facet_grid(rows=vars(pop))+
  theme(panel.spacing=unit(0,'pt'))

################################################################################
# Fig. S8e

figs8e_plot=ggplot()+
geom_rect(data=data.table(xmin=770,xmax=970,ymin=-20,ymax=+20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="#71458d",alpha=0.2)+
  geom_line(data=figs8d_data[random==T,],mapping=aes(epoch,dotallele_norm_smooth*1000,group=rsID,color=pop),alpha=0.2,color="gray")+
  geom_line(data=figs8d_data[random==F,],mapping=aes(epoch,dotallele_norm_smooth*1000,group=rsID,color=pop))+
  scale_color_manual(values=color_populations_1kg)+
coord_cartesian(ylim=c(-3.5,4.5))+
  xlab("Generations from present")+ylab("d/dt(f)*1000")+
  facet_grid(rows=vars(pop))+
  theme(panel.spacing=unit(0,'pt'))

################################################################################
# Fig. S8f

figs8f_plot=ggplot()+
geom_rect(data=data.table(xmin=770,xmax=970,ymin=-20,ymax=20),mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="#71458d",alpha=0.2)+
  geom_line(data=figs8d_data[random==T,],mapping=aes(epoch,z_smooth,color=pop,group=rsID),alpha=0.2,color="gray")+
  geom_line(data=figs8d_data[random==F,],mapping=aes(epoch,z_smooth,color=pop,group=rsID))+
  geom_hline(yintercept=3,linetype='dashed',size=0.1)+
  geom_hline(yintercept=-3,linetype='dashed',size=0.1)+
geom_segment(data=data.table(x=c(721,1203),y=c(3,3),yend=c(-20,-20)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
geom_segment(data=data.table(x=c(951,1634),y=c(3,3),yend=c(-20,-20)),mapping=aes(x=x,xend=x,y=y,yend=yend),size=0.1,linetype='dashed')+
  scale_color_manual(values=color_populations_1kg)+
  scale_y_continuous(breaks=c(-3,0,3))+
coord_cartesian(ylim=c(-5,10))+
  xlab("Generations from present")+ylab("Z-score")+
  facet_grid(rows=vars(pop))+
  theme(panel.spacing=unit(0,'pt'))

################################################################################
pname=sprintf("%s/FigS8/FigS8.pdf",FIG_DIR)
pdf(pname,width=7.2,height=6.7)
  grid.arrange(
    grobs=list(
      ggplotGrob(figs8a_plot+theme(legend.position="none")),
      ggplotGrob(figs8b_plot+theme(legend.position="none")),
      ggplotGrob(figs8c_plot+theme(legend.position="none")),
      ggplotGrob(figs8d_plot+theme(legend.position="none",axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())),
      ggplotGrob(figs8e_plot+theme(legend.position="none",axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())),
      ggplotGrob(figs8f_plot+theme(legend.position="none",axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())),
      grid.rect(gp=gpar(col="white"))
    ),
    layout_matrix=rbind(
      c(4,4,1,1),
      c(5,5,1,1),
      c(6,6,1,1),
      c(3,3,2,2),
      c(3,3,2,2),
      c(3,3,2,2),
      c(7,7,7,7)
    ), heights=c(1.5,1.5,1.5,1.5,1.5,1.5,1), widths=c(1,3.5,1.7,2.8)
  )
dev.off()
