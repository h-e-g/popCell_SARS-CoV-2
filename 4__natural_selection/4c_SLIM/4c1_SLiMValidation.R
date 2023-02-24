
################################################################################
################################################################################
# File name: 4c1_SLiMValidation.R
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: simulate data under various scenarios of natural selection
# and date selective events based on frequency trajectory
# Effector script/ self launching (with 00_Rscript)
################################################################################
################################################################################
# commands for launching script
# SCRIPT_DIR="../../../4__natural_selection/4c_SLIM/"
# LOG_DIR="../../../LOG"
# MISC_DIR=""../../../MISC"

# sbatch --array=1-3939 --parsable --mem=4G -o ${LOG_DIR}/logfile_4c1_%A_%a.log -J SLiM ${MISC_DIR}/00_Rscript.sh ${SCRIPT_DIR}/4c1_SLiMValidation.R
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

#### read command line args
cmd=commandArgs()
print(cmd)
for (i in 1:length(cmd)){
	if (cmd[i]=='--set' | cmd[i]=='-s' ){TASK_ID = as.numeric(cmd[i+1])} # ID of the set to test
}

TASK_ID=as.numeric(TASK_ID)

library(data.table)

PARAMS=fread(sprintf("%s/Simulation_parameters.tsv",DAT_SLIM_DIR))
PARAMS[,RUN:=1:.N]

POP=PARAMS[TASK_ID,POP]

Demography=fread(sprintf("%s/spiedel_popsize_%s.coal",DAT_SLIM_DIR,POP,POP))
Ngen=round(as.numeric(names(Demography)))
Neff=round(0.5/as.numeric(unlist(Demography)))
DEMO=data.table(start=Ngen[2:31],end=Ngen[1:30],Neff=Neff[3:32])

MAX_GENERATION=DEMO[start>4000,][which.min(start),round(start/10)]


######################################################
########### define operating blocks #################
######################################################

### 1. initialize simulations
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
 }'

### 2. define population
block_pop_setup='1 first() {
        sim.addSubpop("p1", VAR_NE_T0);
      }'

### 3. create target variant with initial frequency VAR_FREQ0
block_mut_setup='1 late() {
          Ne=asInteger(length(p1.genomes));
          NumberHapsWithMutation = asInteger(round(1e-12+min(max((Ne*VAR_FREQ0)+1e-12,1.),Ne-1.)));
           #// define target haplotype with freq VAR_FREQ0
          target = sample(p1.genomes, NumberHapsWithMutation);
          #// add new mutation at position 1
          target.addNewDrawnMutation(m2, 1);
      }'

### 4. change population size to VAR_TARGET_POPSIZE at a generation VAR_GENERATION_CHANGE
block_size_change='VAR_GENERATION_CHANGE first() {
                     p1.setSubpopulationSize(VAR_TARGET_POPSIZE);
                     }'


### 5. print frequency until VAR_NGENERATION_MAX
block_print_freq='1:VAR_NGENERATION_MAX late() {
            muts=sim.mutations;
    	      cat(muts);
            cat(sim.mutationFrequencies(p1,muts)+"\t");
            }'


### 6. start natural selection at VAR_ONSET_SELECTION with coeff VAR_SELECT_COEF
block_selection_start='VAR_ONSET_SELECTION first() {
                            m2muts = sim.mutationsOfType(m2);
                            for (index in seqAlong(m2muts));
                            m2muts[index].setSelectionCoeff(VAR_SELECT_COEF);
                          }'

### 6. end natural selection at VAR_END_SELECTION (coeff= 0)
block_selection_end='VAR_END_SELECTION first() {
                            m2muts = sim.mutationsOfType(m2);
                            for (index in seqAlong(m2muts));
                            m2muts[index].setSelectionCoeff(0.);
                          }'

block_finish='VAR_NGENERATION_MAX late() {
                             sim.simulationFinished();
                           }'

# initialize demography
demography_block=block_pop_setup
# update demography to fit speidel et al
for (i in DEMO[,which(start<4000)]){
  block_size_change_i=gsub('VAR_GENERATION_CHANGE',as.character(DEMO[i,round(MAX_GENERATION-start/10)]),block_size_change)
  block_size_change_i=gsub('VAR_TARGET_POPSIZE',as.character(DEMO[i,round(Neff/10)]),block_size_change_i)
  demography_block=generateScript(demography_block,block_size_change_i,quiet=TRUE)
}


# initialize mutation
selection_block=block_mut_setup
# update mutation to add selection at the desired time
if(PARAMS[TASK_ID,COEF_SELECTION]!=0){
  selection_block=generateScript(selection_block,
                            block_selection_start,
                            block_selection_end,
                            quiet=TRUE)
  selection_block=gsub('VAR_ONSET_SELECTION',as.character(max(2,PARAMS[TASK_ID,ONSET_SELECTION])),selection_block)
  selection_block=gsub('VAR_END_SELECTION',as.character(PARAMS[TASK_ID,END_SELECTION]),selection_block)
  selection_block=gsub('VAR_SELECT_COEF',as.character(PARAMS[TASK_ID,COEF_SELECTION]*10),selection_block)
}

# define simulation structure
myscript=generateScript(block_init,
                          demography_block,
                          selection_block,
                          block_print_freq,
                          block_finish,quiet=TRUE)

# set parameters
myscript=gsub('VAR_NGENERATION_MAX',as.character(MAX_GENERATION),myscript)
myscript=gsub('VAR_NE_T0',as.character(DEMO[start>4000,][which.min(start),round(Neff/10)]),myscript)
myscript=gsub('VAR_FREQ0',as.character(PARAMS[TASK_ID,FREQ0]),myscript)

NREP=PARAMS[TASK_ID,NREP]

# run simulations NREP times
DT_res=list()
for (rep in 1:NREP){
  tic()
  output=runScript(myscript)
  Freqs=as.numeric(gsub('Mutation<0:.*>','',strsplit(output[15],'\t')[[1]]))
	# pdf(sprintf("%s/testFreq.pdf",EIP))
	# plot(10*1:496,Freqs,type='l',ylim=c(0,1));abline(v=c(1000,1200),col='grey');
	# dev.off()
  DT_res[[rep]]=data.table(Freq=Freqs,generation=1:length(Freqs),generation_ago=MAX_GENERATION-(1:length(Freqs)),RUN=TASK_ID)
  toc()
}
# write results
DT_res=rbindlist(DT_res,idcol='replicate')
PARAMS[,RUN:=1:.N]
DT_res=merge(DT_res,PARAMS,by='RUN')
fwrite(DT_res,sprintf("%s/results/results_simTask%s.tsv.gz",DAT_SLIM_DIR,TASK_ID),sep='\t')
