# SCRIPT_DIR="/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/template_scripts/processing_pipeline"
# sbatch --array=1-6675 --parsable --mem=4G --qos=fast -o ${SCRIPT_DIR}/log/logfile_slim_%A_%a.log -J SLiM ${SCRIPT_DIR}/00_Rscript.sh ./18_SLiMValidation_v3.R --set
# sbatch --array=1-101 --parsable --mem=8G -o ${SCRIPT_DIR}/log/logfile_slim_%A_%a.log -J SLiM ${SCRIPT_DIR}/00_Rscript.sh ./18_SLiMValidation_v3.R --set
# sbatch --array=102-2225,2327-4450,4552-6675 --parsable --mem=4G --qos=fast -o ${SCRIPT_DIR}/log/logfile_slim_%A_%a.log -J SLiM ${SCRIPT_DIR}/00_Rscript.sh ./18_SLiMValidation_v3.R --set
# sbatch --array=1920-2225,4145-4450,6370-6675 --parsable --mem=4G --qos=fast -o ${SCRIPT_DIR}/log/logfile_slim_%A_%a.log -J SLiM ${SCRIPT_DIR}/00_Rscript.sh ./18_SLiMValidation_v3.R --set
# sbatch --array=6676-8493 --parsable --mem=4G --qos=fast -o ${SCRIPT_DIR}/log/logfile_slim_%A_%a.log -J SLiM ${SCRIPT_DIR}/00_Rscript.sh ./18_SLiMValidation_v3.R --set
# sbatch --array=8494-10917 --parsable --mem=4G --qos=fast -o ${SCRIPT_DIR}/log/logfile_slim_%A_%a.log -J SLiM ${SCRIPT_DIR}/00_Rscript.sh ./18_SLiMValidation_v3.R --set

options(stringsAsFactors=FALSE, max.print=9999, width=300, datatable.fread.input.cmd.message=FALSE)
.libPaths("/pasteur/zeus/projets/p02/evo_immuno_pop/single_cell/resources/R_libs/4.1.0")
library(tictoc)

cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
	if (cmd[i]=='--set' | cmd[i]=='-s' ){TASK_ID = as.numeric(cmd[i+1])} # ID of the set to test
}

TASK_ID=as.numeric(TASK_ID)

library(data.table)

EIP='/pasteur/zeus/projets/p02/evo_immuno_pop'
#
# PARAMS_POP=list()
# for (POP in c('CEU','CHS','YRI')){
#   PARAMS=list()
#   PARAMS[[1]]=data.table(ONSET_SELECTION=0,END_SELECTION=0,COEF_SELECTION=0,NREP=1000,FREQ0=seq(0, 1,by=0.01))
#   PARAMS[[2]]=data.table(ONSET_SELECTION=100,END_SELECTION=120,COEF_SELECTION=0.005,NREP=100,FREQ0=seq(0, 1,by=0.01))
#   PARAMS[[3]]=data.table(ONSET_SELECTION=100,END_SELECTION=120,COEF_SELECTION=0.01,NREP=100,FREQ0=seq(0, 1,by=0.01))
#   PARAMS[[4]]=data.table(ONSET_SELECTION=100,END_SELECTION=120,COEF_SELECTION=0.05,NREP=100,FREQ0=seq(0, 1,by=0.01))
#   PARAMS[[5]]=data.table(ONSET_SELECTION=200,END_SELECTION=220,COEF_SELECTION=0.005,NREP=100,FREQ0=seq(0, 1,by=0.01))
#   PARAMS[[6]]=data.table(ONSET_SELECTION=200,END_SELECTION=220,COEF_SELECTION=0.01,NREP=100,FREQ0=seq(0, 1,by=0.01))
#   PARAMS[[7]]=data.table(ONSET_SELECTION=200,END_SELECTION=220,COEF_SELECTION=0.05,NREP=100,FREQ0=seq(0, 1,by=0.01))
#   PARAMS[[8]]=data.table(ONSET_SELECTION=300,END_SELECTION=320,COEF_SELECTION=0.005,NREP=100,FREQ0=seq(0, 1,by=0.01))
#   PARAMS[[9]]=data.table(ONSET_SELECTION=300,END_SELECTION=320,COEF_SELECTION=0.01,NREP=100,FREQ0=seq(0, 1,by=0.01))
#   PARAMS[[10]]=data.table(ONSET_SELECTION=300,END_SELECTION=320,COEF_SELECTION=0.05,NREP=100,FREQ0=seq(0, 1,by=0.01))
#   PARAMS[[11]]=data.table(ONSET_SELECTION=100,END_SELECTION=120,COEF_SELECTION=-0.005,NREP=100,FREQ0=seq(0, 1,by=0.01))
#   PARAMS[[12]]=data.table(ONSET_SELECTION=100,END_SELECTION=120,COEF_SELECTION=-0.01,NREP=100,FREQ0=seq(0, 1,by=0.01))
#   PARAMS[[13]]=data.table(ONSET_SELECTION=100,END_SELECTION=120,COEF_SELECTION=-0.05,NREP=100,FREQ0=seq(0, 1,by=0.01))
#   PARAMS[[14]]=data.table(ONSET_SELECTION=200,END_SELECTION=220,COEF_SELECTION=-0.005,NREP=100,FREQ0=seq(0, 1,by=0.01))
#   PARAMS[[15]]=data.table(ONSET_SELECTION=200,END_SELECTION=220,COEF_SELECTION=-0.01,NREP=100,FREQ0=seq(0, 1,by=0.01))
#   PARAMS[[16]]=data.table(ONSET_SELECTION=200,END_SELECTION=220,COEF_SELECTION=-0.05,NREP=100,FREQ0=seq(0, 1,by=0.01))
#   PARAMS[[17]]=data.table(ONSET_SELECTION=300,END_SELECTION=320,COEF_SELECTION=-0.005,NREP=100,FREQ0=seq(0, 1,by=0.01))
#   PARAMS[[18]]=data.table(ONSET_SELECTION=300,END_SELECTION=320,COEF_SELECTION=-0.01,NREP=100,FREQ0=seq(0, 1,by=0.01))
#   PARAMS[[19]]=data.table(ONSET_SELECTION=300,END_SELECTION=320,COEF_SELECTION=-0.05,NREP=100,FREQ0=seq(0, 1,by=0.01))
#   PARAMS[[20]]=data.table(ONSET_SELECTION=1,END_SELECTION=496,COEF_SELECTION=0.001,NREP=100,FREQ0=seq(0, 0.5,by=0.01))
#   PARAMS[[21]]=data.table(ONSET_SELECTION=1,END_SELECTION=496,COEF_SELECTION=0.0025,NREP=100,FREQ0=seq(0,0.5,by=0.01))
#   PARAMS[[22]]=data.table(ONSET_SELECTION=1,END_SELECTION=496,COEF_SELECTION=0.005,NREP=100,FREQ0=seq(0,0.5,by=0.01))
#   PARAMS[[23]]=data.table(ONSET_SELECTION=1,END_SELECTION=496,COEF_SELECTION=-0.001,NREP=100,FREQ0=seq(0.5,1,by=0.01))
#   PARAMS[[24]]=data.table(ONSET_SELECTION=1,END_SELECTION=496,COEF_SELECTION=-0.0025,NREP=100,FREQ0=seq(0.5,1,by=0.01))
#   PARAMS[[25]]=data.table(ONSET_SELECTION=1,END_SELECTION=496,COEF_SELECTION=-0.005,NREP=100,FREQ0=seq(0.5, 1,by=0.01))
#   PARAMS=rbindlist(PARAMS)
# PARAMS_POP[[POP]]=PARAMS
# }
# PARAMS_POP=rbindlist(PARAMS_POP,idcol='POP')
#
# PARAMS_POP2=list()
# for (POP in c('CEU','CHS','YRI')){
# PARAMS_POP2[[POP]]=list()
# PARAMS_POP2[[POP]][[1]]=data.table(ONSET_SELECTION=400,END_SELECTION=420,COEF_SELECTION=-0.005,NREP=100,FREQ0=seq(0, 1,by=0.01))
# PARAMS_POP2[[POP]][[2]]=data.table(ONSET_SELECTION=400,END_SELECTION=420,COEF_SELECTION=-0.01,NREP=100,FREQ0=seq(0, 1,by=0.01))
# PARAMS_POP2[[POP]][[3]]=data.table(ONSET_SELECTION=400,END_SELECTION=420,COEF_SELECTION=-0.05,NREP=100,FREQ0=seq(0, 1,by=0.01))
# PARAMS_POP2[[POP]][[4]]=data.table(ONSET_SELECTION=400,END_SELECTION=420,COEF_SELECTION=0.005,NREP=100,FREQ0=seq(0, 1,by=0.01))
# PARAMS_POP2[[POP]][[5]]=data.table(ONSET_SELECTION=400,END_SELECTION=420,COEF_SELECTION=0.01,NREP=100,FREQ0=seq(0, 1,by=0.01))
# PARAMS_POP2[[POP]][[6]]=data.table(ONSET_SELECTION=400,END_SELECTION=420,COEF_SELECTION=0.05,NREP=100,FREQ0=seq(0, 1,by=0.01))
# PARAMS_POP2[[POP]]=rbindlist(PARAMS_POP2[[POP]])
# }
# PARAMS_POP2=rbindlist(PARAMS_POP2,idcol='POP')
#
# PARAMS=rbind(PARAMS_POP,PARAMS_POP2)
# PARAMS_POP3=list()
# for (POP in c('CEU','CHS','YRI')){
# PARAMS_POP3[[POP]]=list()
# PARAMS_POP3[[POP]][[1]]=data.table(ONSET_SELECTION=100,END_SELECTION=120,COEF_SELECTION=-0.025,NREP=100,FREQ0=seq(0, 1,by=0.01))
# PARAMS_POP3[[POP]][[2]]=data.table(ONSET_SELECTION=200,END_SELECTION=220,COEF_SELECTION=-0.025,NREP=100,FREQ0=seq(0, 1,by=0.01))
# PARAMS_POP3[[POP]][[3]]=data.table(ONSET_SELECTION=300,END_SELECTION=320,COEF_SELECTION=-0.025,NREP=100,FREQ0=seq(0, 1,by=0.01))
# PARAMS_POP3[[POP]][[4]]=data.table(ONSET_SELECTION=400,END_SELECTION=420,COEF_SELECTION=-0.025,NREP=100,FREQ0=seq(0, 1,by=0.01))
# PARAMS_POP3[[POP]][[5]]=data.table(ONSET_SELECTION=100,END_SELECTION=120,COEF_SELECTION=0.025,NREP=100,FREQ0=seq(0, 1,by=0.01))
# PARAMS_POP3[[POP]][[6]]=data.table(ONSET_SELECTION=200,END_SELECTION=220,COEF_SELECTION=0.025,NREP=100,FREQ0=seq(0, 1,by=0.01))
# PARAMS_POP3[[POP]][[7]]=data.table(ONSET_SELECTION=300,END_SELECTION=320,COEF_SELECTION=0.025,NREP=100,FREQ0=seq(0, 1,by=0.01))
# PARAMS_POP3[[POP]][[8]]=data.table(ONSET_SELECTION=400,END_SELECTION=420,COEF_SELECTION=0.025,NREP=100,FREQ0=seq(0, 1,by=0.01))
# PARAMS_POP3[[POP]]=rbindlist(PARAMS_POP3[[POP]])
# }
# PARAMS_POP3=rbindlist(PARAMS_POP3,idcol='POP')
#  PARAMS=rbind(PARAMS,PARAMS_POP3)
#  PARAMS[,RUN:=1:.N]
# fwrite(PARAMS,file=sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulation_parameters.tsv",EIP))
PARAMS=fread(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulation_parameters.tsv",EIP))
POP=PARAMS[TASK_ID,POP]

Demography=fread(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/relate/%s/spiedel_popsize_%s.coal",EIP,POP,POP))
Ngen=round(as.numeric(names(Demography)))
Neff=round(0.5/as.numeric(unlist(Demography)))
DEMO=data.table(start=Ngen[2:31],end=Ngen[1:30],Neff=Neff[3:32])

MAX_GENERATION=DEMO[start>4000,][which.min(start),round(start/10)]

######### CHANGE VALUE OF VARIABLE WITHIN A SLIM BLOCK  #########
udpate_block=function(block,...){
  args=match.call(expand.dots = TRUE)
  PARAMS_NAME=as.character(names(args[-(1:2)]))
  PARAMS_VALUE=as.character(sapply(args[-(1:2)],eval))
  for (i in 1: length(PARAMS_NAME)){
    block=gsub(PARAMS_NAME[i],PARAMS_VALUE[i],block)
  }
  block
 }

######### CONCATENATE SLIM BLOCKS #########
generateScript=function(...,quiet=FALSE){
  args=match.call(expand.dots = TRUE)
  args=args[as.character(names(args))!='quiet']
  BLOCKS=sapply(args[-(1)],function(x){as.character(eval(x),envir=parent.frame())})
  temp_script=sprintf('%s.slim',tempfile())
  if(!quiet){
    cat(paste0(BLOCKS,'\n',collapse='\n\n'))
  }
  paste0(BLOCKS,'\n',collapse='\n\n')
}

######### RUN SCRIPT AND READ OUTPUT #########
runScript=function(SCRIPT,
                    scriptFILE=sprintf('%s.slim',tempfile())
                    ){
  SCRIPT_clean=gsub('#','//',SCRIPT)
  SCRIPT_clean=gsub('\"\t\"','\"\\\\t\"',SCRIPT_clean)

  write(SCRIPT_clean,file=scriptFILE)
  SLIM_CMD=paste('slim',scriptFILE)
  out=system(SLIM_CMD,intern=TRUE)
  out
}


######################################################
########### define operatiing blocks #################
######################################################

###1. initialize smulations
block_init <- 'initialize() {
	initializeMutationRate(0);
	initializeRecombinationRate(0);
  #// neutral (not used)
	initializeMutationType("m1", 0.5, "f", 0.0);
  #// the one to be selected (for now, neutral)
	initializeMutationType("m2",  0.5, "f", 0);
	m2.convertToSubstitution = F;
# // define SNP region, with target mutation
	initializeGenomicElementType("g1", m1, 1.0); #
# // define a single SNP region of size 2
	initializeGenomicElement(g1, 0, 2);
#	// Short examples of various selection model
# //initializeMutationType("m2",  0.1, "f", 0.5); //positively selected mutation (codominance)
# //initializeMutationType("m2",  -0.1, "f", 0.5); //negatively selected mutation (codominance)
# //initializeMutationType("m2",  0, "f", -0.5); //negatively selected mutation (only on homozygotes)
 }'

 ###1. define population
block_pop_setup='1 first() {
        #// Expand the African population to 14474
        #// This occurs 148000 years (5920) generations ago
        sim.addSubpop("p1", VAR_NE_T0);
      }'

###1. create target variant
block_mut_setup='1 late() {
          Ne=asInteger(length(p1.genomes));
          NumberHapsWithMutation = asInteger(round(1e-12+min(max((Ne*VAR_FREQ0)+1e-12,1.),Ne-1.)));
           #// define target haplotype with freq VAR_FREQ0
          target = sample(p1.genomes, NumberHapsWithMutation);
          #// add new mutation at position 1
          target.addNewDrawnMutation(m2, 1);
      }'

block_size_change='VAR_GENERATION_CHANGE first() {
                     p1.setSubpopulationSize(VAR_TARGET_POPSIZE);
                     }'


block_print_freq='1:VAR_NGENERATION_MAX late() {
            muts=sim.mutations;
    	      cat(muts);
            cat(sim.mutationFrequencies(p1,muts)+"\t");
            }'


block_selection_start='VAR_ONSET_SELECTION first() {
                            m2muts = sim.mutationsOfType(m2);
                            for (index in seqAlong(m2muts));
                            m2muts[index].setSelectionCoeff(VAR_SELECT_COEF);
                          }'

block_selection_end='VAR_END_SELECTION first() {
                            m2muts = sim.mutationsOfType(m2);
                            for (index in seqAlong(m2muts));
                            m2muts[index].setSelectionCoeff(0.);
                          }'

block_selection_gaspard='VAR_ONSET_SELECTION:VAR_END_SELECTION fitness(m2) {
															if(sim.generation >= VAR_ONSET_SELECTION & sim.generation <=VAR_END_SELECTION) {
																if(VAR_SELECT_COEF>0.1) {
																	ftmp = sim.mutationFrequencies(p2, mut);
																	for(i in 1:10) {
													    				ftmp = sel(ftmp,VAR_SELECT_COEF);
													    			}
													    			h = hCalc(sim.mutationFrequencies(p2, mut),ftmp);
													    			if(ftmp>0){
																			return 1.0 + h * (1.0);
																	}else{
																		return 1.0;
													    			}
																}else{
											          		return 1.0 + 0.5 * (VAR_SELECT_COEF)*10;
													      		}
															}else{
																return 1.0;
															}
													}'


block_finish='VAR_NGENERATION_MAX late() {
                             sim.simulationFinished();
                           }'


# DONE inverse the demography to account for reversed time between relate and SLIM
# DONE run in batch for TASK_ID in 1-6675 (avant le repas)

demography_block=block_pop_setup
for (i in DEMO[,which(start<4000)]){
  block_size_change_i=gsub('VAR_GENERATION_CHANGE',as.character(DEMO[i,round(MAX_GENERATION-start/10)]),block_size_change)
  block_size_change_i=gsub('VAR_TARGET_POPSIZE',as.character(DEMO[i,round(Neff/10)]),block_size_change_i)
	#block_size_change_i=gsub('VAR_GENERATION_CHANGE',as.character(DEMO[i,round(MAX_GENERATION-start)]),block_size_change)
	#block_size_change_i=gsub('VAR_TARGET_POPSIZE',as.character(DEMO[i,round(Neff)]),block_size_change_i)
  demography_block=generateScript(demography_block,block_size_change_i,quiet=TRUE)
}

selection_block=block_mut_setup
if(PARAMS[TASK_ID,COEF_SELECTION]!=0){
  selection_block=generateScript(selection_block,
                            block_selection_start,
                            block_selection_end,
                            quiet=TRUE)
  selection_block=gsub('VAR_ONSET_SELECTION',as.character(max(2,PARAMS[TASK_ID,ONSET_SELECTION])),selection_block)
  selection_block=gsub('VAR_END_SELECTION',as.character(PARAMS[TASK_ID,END_SELECTION]),selection_block)
# selection_block=gsub('VAR_ONSET_SELECTION',as.character(1000),selection_block)
# selection_block=gsub('VAR_END_SELECTION',as.character(1200),selection_block)
  selection_block=gsub('VAR_SELECT_COEF',as.character(PARAMS[TASK_ID,COEF_SELECTION]*10),selection_block)
}
# DONE: run script


myscript=generateScript(block_init,
                          demography_block,
                          selection_block,
                          block_print_freq,
                          block_finish,quiet=TRUE)

myscript=gsub('VAR_NGENERATION_MAX',as.character(MAX_GENERATION),myscript)
myscript=gsub('VAR_NE_T0',as.character(DEMO[start>4000,][which.min(start),round(Neff/10)]),myscript)
#myscript=gsub('VAR_NE_T0',as.character(DEMO[start>4000,][which.min(start),round(Neff)]),myscript)
myscript=gsub('VAR_FREQ0',as.character(PARAMS[TASK_ID,FREQ0]),myscript)

NREP=PARAMS[TASK_ID,NREP]

DT_res=list()
for (rep in 1:NREP){
  tic()
  output=runScript(myscript)
  Freqs=as.numeric(gsub('Mutation<0:.*>','',strsplit(output[15],'\t')[[1]]))
	# pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/testFreq_NEW.pdf",EIP))
	# plot(10*1:496,Freqs,type='l',ylim=c(0,1));abline(v=c(1000,1200),col='grey');
	# dev.off()
	# pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/testFreq_old.pdf",EIP))
	# plot(10*1:496,Freqs,type='l',ylim=c(0,1));abline(v=c(1000,1200),col='grey');
	# dev.off()
	# pdf(sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/plots/testFreq.pdf",EIP))
	# plot(1:4962,Freqs,type='l');abline(v=c(1000,1200),col='grey');
	# dev.off()
  DT_res[[rep]]=data.table(Freq=Freqs,generation=1:length(Freqs),generation_ago=MAX_GENERATION-(1:length(Freqs)),RUN=TASK_ID)
  toc()
}
DT_res=rbindlist(DT_res,idcol='replicate')
PARAMS[,RUN:=1:.N]
DT_res=merge(DT_res,PARAMS,by='RUN')
fwrite(DT_res,sprintf("%s/single_cell/project/pop_eQTL/data/4_natural_selection/clues/Simulations/results_simTask%s.tsv.gz",EIP,TASK_ID),sep='\t')

#DONE: simulate selection
# DONE: read demography and simulate accordingly
# DONE: readoutput from SlIM
# DONE: automate
# DONE run in batch
