################################################################################
################################################################################
# File name: 4c2_SLiMValidation_aggregate.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: aggregate simulated data from various scenarios of natural selection
# effector script
################################################################################

# Setup

# load required packages
LIB_DIR="../../../LIBRARY"
source(sprintf("%s/00_one_lib_to_rule_them_all.R",LIB_DIR))

# declare shortcuts
MISC_DIR="../../../MISC"
source(sprintf("./shortcuts.R",MISC_DIR))

# declare SLIM block functions
source(sprintf(MISC_DIR,'%s/misc_SLiM.R'))

PARAMS=fread(sprintf("%s/Simulation_parameters.tsv",DAT_SLIM_DIR))
PARAMS[,RUN:=1:.N]

result_dir=sprintf("%s/results/",DAT_SLIM_DIR)
result_files=dir(result_dir,pattern='results')

# declare fiunction/object for plotting
source(sprintf("%s/misc_plots.R",MISC_DIR))

##########@@ plot population size and selection
MAX_GENERATION=4962
PARAMS_SIM=data.table(generation_ago=0:MAX_GENERATION)
for(POP in c('CEU','CHS','YRI')){
  Demography=fread(sprintf("%s/spiedel_popsize_%s.coal",DAT_SLIM_DIR,POP,POP))
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
  for(COEF in c(0.01,0.025,0.05)){
    select_group=sprintf('s%s_+_t%s-%s',COEF,round(ONSET),round(ONSET)+20)
    PARAMS_SEL[[select_group]]=data.table(generation_ago=0:MAX_GENERATION,
                                      select = ifelse(MAX_GENERATION:0>=ONSET*10 & MAX_GENERATION:0<=(10*ONSET+200),COEF,0),
                                      onset=ONSET,
                                      select_coeff=COEF)
    }
  }
PARAMS_SEL=rbindlist(PARAMS_SEL,idcol='select_group')


all_results_SLiM=list()
for (task in PARAMS){
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

# compute Zscores
all_results_SLiM[,Freq_smooth:=pmin(1,pmax(0,predict(loess(Freq~generation,span=0.1),generation))),keyby=.(RUN,replicate)]
all_results_SLiM[,DeltaFreq_smooth:=c(NA,diff(Freq_smooth)),keyby=.(RUN,replicate)]
all_results_SLiM[,DeltaFreq_smooth_Norm:=DeltaFreq_smooth/sqrt(Freq*(1-Freq)),keyby=.(RUN,replicate)]

NormFactors_smooth=all_results_SLiM[COEF_SELECTION==0 & paste(replicate,RUN)%chin%notfixed,.(SD_smooth=sd(DeltaFreq_smooth_Norm,na.rm=TRUE)),keyby=.(POP,generation,generation_ago)]

all_results_SLiM=merge(all_results_SLiM,NormFactors_smooth,by=c('POP','generation','generation_ago'))
all_results_SLiM[,Zscore_smooth:=DeltaFreq_smooth_Norm/SD_smooth]

all_results_SLiM[,SELECTED:=generation>=ONSET_SELECTION & generation<=END_SELECTION]
all_results_SLiM[,simID:=paste(RUN,replicate)]

fwrite(all_results_SLiM,file=sprintf("%s/results/0_all_results_SLiM.tsv.gz",DAT_SLIM_DIR),sep='\t')

#### find optimal threshold
Zth=unique(c(seq(0,2,by=.1),seq(2,4,by=.01),seq(4,6,by=.1)))

TH_dt=data.table(expand.grid(simID=all_results_SLiM[generation==310,simID], TH=Zth))
Power_thChoice_310gen=merge(all_results_SLiM[generation==310,],TH_dt,by='simID')
Power_thChoice_310gen=Power_thChoice_310gen[,.(Zcore_power=mean(abs(Zscore)>TH), Zcore_smooth_power=mean(abs(Zscore_smooth)>TH)), by=.(TH, generation, generation_ago, COEF_SELECTION, ONSET_SELECTION, END_SELECTION, SELECTED, POP)]
Power_thChoice_310gen[,select_group:=paste0('s',abs(COEF_SELECTION),'_',ifelse(COEF_SELECTION>0,'+','-'),'_t',ONSET_SELECTION,'-',END_SELECTION)]

Power_thChoice2=merge(Power_thChoice_310gen[ONSET_SELECTION==0,.(TH,generation,generation_ago,POP,Zcore_smooth_power)],Power_thChoice_310gen[ONSET_SELECTION>0,],by=c('generation','POP','TH','generation_ago'),suffix=c('.neutral',''))
fwrite(Power_thChoice2,file=sprintf("%s/results/1_Power_thChoice.tsv.gz",EIP),sep='\t')

# estimate Power at threshold of |Z|>3 for various scenarios
Power_th3=all_results_SLiM[,.(Zcore_smooth_power=mean(abs(Zscore_smooth)>3)), by=.(generation, generation_ago, COEF_SELECTION, ONSET_SELECTION, END_SELECTION, SELECTED, POP)]
Power_th3[,select_group:=paste0('s',abs(COEF_SELECTION),'_',ifelse(COEF_SELECTION>0,'+','-'),'_t',ONSET_SELECTION,'-',END_SELECTION)]
Power_th3[,TH:=3]

fwrite(Power_th3,file=sprintf("%s/results/2_Power_Zth3.tsv.gz",DAT_SLIM_DIR),sep='\t')
