import sims
import numpy as np
import msprime
import pandas as pd

import useful
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
from importlib import reload


import json



import sys

pref=sys.argv[1]

dir=sys.argv[2]



# Generation time, mutation rate and recomination rate
RR = 1e-8
MU = 1.29e-8 
GEN_time = 29.0 

# Split Times
T_NEAND_migration = 55000 #time of Neanderthal migration into Out_of_africa population
T_NEAND_AMH = 550000 # split time between AMH and Neanderthal
T_OOF_AF = 65700 # Out_of_Africa migration time
T_NEAND_samples = 38000

# Effective population size
N_ANC = 18500 # N_e of common  AMH and NEanderthal population 
N_ND = 3400 # N_e of Neanderthal
N_AMH = 23000 # N_e of AMH
N_OOF = 1861 # N_e of Out_of_Africa population
N_AF = 27600 # N_e of Africans
N_EU = 13377 #N_e of Europeans

N_EU_bottleneck = 1080
N_EU_growth = 1450
T_EU_growth = 31900
gr_rate = 0.00202
Portion_admix = 0.05

len_sequence = 50000000 # DNA sequence length

n = 250 # number of generated   AF samples
n_neand = 6 #number of generated Neanderthals

rand_sd =1234 #random seed
T = np.array([T_NEAND_migration, T_NEAND_AMH, T_OOF_AF])/GEN_time


N_ND = 3400 # N_e of Neanderthal
N_e = np.array([N_ANC, N_ND, N_AMH, N_OOF, N_AF, N_EU])

n_eu=1
ploidy=2
L=1000
N_ref_pop=250
N_neanderthal=6

number_chr=1
ts_mas, ND_true_tracts_mas = sims.sim_history_archaic_many_ts(GEN_time, len_sequence, RR, MU, N_e, T,  n, rand_sd, n_neand,  
                              T_NEAND_samples/GEN_time, n_eu, N_EU_growth, 
                              T_EU_growth/GEN_time, N_EU_bottleneck, gr_rate, Portion_admix, ploidy, number_chr)



#functions that move the j-th set of tracts by j*len_sequence
def many_in_one(nd_tracts_mas, number_chr):
    nd2=[]
    for _ in range(len(ts_mas)):
        for j in range(len(nd_tracts_mas[_])):
            nd2.append([nd_tracts_mas[_][j][0]+_*len_sequence, nd_tracts_mas[_][j][1]+_*len_sequence])

    return [useful.tracts_eu(nd2, len_sequence*len(ts_mas) ), nd2]


true_nd=[ND_true_tracts_mas[_][0] for _ in range(len(ts_mas))]  
real_tracts_in_states = many_in_one(true_nd,number_chr)

np.save('./'+dir+'/'+pref+'.ND_true_tracts.npy', np.array(ND_true_tracts_mas, dtype= object))

cover_mas=[0.25,0.4, 0.5, 0.65, 0.8, 0.9, 0.999]

for cover in cover_mas:
    for _ in range(len(ts_mas)):
        sims.create_obs_txt_coverage('./'+dir+'/'+'obs.chr'+pref+'.cover.'+str(cover)+'.txt', ploidy,  n_eu, ts_mas[_], N_ref_pop, N_neanderthal, n, cover)
    sims.create_bed_smpls_arch_cov('./'+dir+'/samples.txt', './'+dir+'/test.bed','./'+dir+'/arch.cover.'+str(cover)+'.txt', len_sequence, cover, 1, n_eu, ploidy)



#preparation_to_hmmix (ancestral alleles, fasta, vcf)

j_son={"ingroup":['tsk_'+str(_) for _ in range(0, n_eu)], "outgroup" : ['tsk_'+str(_) for _ in range(n_eu, n_eu+N_ref_pop)]}
with open('./'+dir+'/individuals.json', 'w') as f:
    json.dump(j_son, f)

for _ in range(len(ts_mas)):
    with open('./'+dir+'/'+'ts.'+'chr'+pref+'.vcf', "w") as vcf_file:
        ts_mas[_].write_vcf(vcf_file, contig_id=pref)

sims.make_ancestral_fasta('./'+dir+'/'+'ancestral.'+pref+'.', ts_mas, len_sequence, pref)

