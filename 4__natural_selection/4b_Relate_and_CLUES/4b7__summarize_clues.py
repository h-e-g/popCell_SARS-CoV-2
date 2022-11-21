################################################################################
################################################################################
# File name: 4b7__summarize_clues.py
# Author: Y.A., M.R., M.ON.
################################################################################
################################################################################
# Step: Calculate posterior allele frequency mean and 95% CI bounds
# Effector script
################################################################################
################################################################################

################################################################################
# Setup

from matplotlib import pyplot as plt
import numpy as np
import argparse
import gzip
import pandas as pd

from os import listdir
from os.path import isfile, isdir

DAT_CLUES_DIR="4__natural_selection/data/CLUES"
#DAT="%s/YRI_CutOff2000/SNPs/eQTL"%DAT_CLUES_DIR
#DAT="%s/CEU_CutOff2000/SNPs/eQTL"%DAT_CLUES_DIR
#DAT="%s/CHS_CutOff2000/SNPs/eQTL"%DAT_CLUES_DIR
#DAT="%s/YRI_CutOff2000/SNPs/reQTL"%DAT_CLUES_DIR
#DAT="%s/CEU_CutOff2000/SNPs/reQTL"%DAT_CLUES_DIR
DAT="%s/CHS_CutOff2000/SNPs/reQTL"%DAT_CLUES_DIR
CLUES_PREFIX="clues_output_"

SNP_list=[f for f in listdir(DAT) if (isdir("%s/%s"%(DAT,f)) and len(listdir("%s/%s"%(DAT,f)))>=5)]

writePath="%s/allele_freq_max__all_SNPs.tsv.gz"%DAT

with gzip.open(writePath, 'at') as f:
    for idx, SNP in enumerate(SNP_list):
        print(idx)
        epochs = np.load('%s/%s/%s%s.epochs.npy'%(DAT,SNP,CLUES_PREFIX,SNP))
        freqs = np.load('%s/%s/%s%s.freqs.npy'%(DAT,SNP,CLUES_PREFIX,SNP))
        logpost = np.load('%s/%s/%s%s.post.npy'%(DAT,SNP,CLUES_PREFIX,SNP))
        all_freq_max = np.dot(freqs,np.exp(logpost))
        lCI = [np.min(freqs[np.cumsum(np.exp(logpost[:,i]))>0.025]) for i in range(0,len(epochs)-1)]
        uCI = [np.max(freqs[np.cumsum(np.exp(logpost[:,i]))<0.975]) for i in range(0,len(epochs)-1)]
        b = np.array([epochs[:-1],all_freq_max,np.asarray([SNP]*all_freq_max.shape[0]),np.asarray(lCI),np.asarray(uCI)])
        b = pd.DataFrame(b.T,columns=['epoch','all_freq_max','SNP','lowerCI','upperCI'])
        bAsString = b.to_string(header=False, index=False)
        f.write(bAsString+'\n')
