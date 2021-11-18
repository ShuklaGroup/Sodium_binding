import numpy as np
import glob
import pickle
import pyemma

#######################################Initialization##############################################

#Macrostates (consist of microstates) with sodium bound in different position 

extracellular_binding_site = [128, 299, 16, 47, 175]
orthosteric_binding_site = [729, 312, 129, 458, 115]
binding_site_I = [33, 361, 7, 758, 777]


#initialization of empty array for mean free passage time between transition 
MFPT_ebs_obs = np.empty(200)
MFPT_obs_bsI = np.empty(200)

#####################################TPT calculations for each bootstrap sample####################

for j in range(200):
    files = pickle.load(open("./../CB2_bootstraping/bt_80_" + str(j) +"_files.pkl",'rb'))
    msm = pickle.load(open("./../CB2_bootstraping/bt_80_" + str(j) +"_msm.pkl",'rb'))
    
    #tpt = pyemma.msm.tpt(msm, extracellular_binding_site, orthosteric_binding_site)
    #MFPT_ebs_obs[j] = tpt.mfpt/10000   #extracellular_binding_site to orthosteric_binding_site

    tpt = pyemma.msm.tpt(msm, orthosteric_binding_site, binding_site_I)
    MFPT_obs_bsI[j] = tpt.mfpt/10000   #orthosteric_binding_site to binding_site_I



print('transition  mean  error \n')

#print("extracellular_binding_site -> orthosteric_binding_site " + str(np.mean(MFPT_ebs_obs)) + ' ' + str(np.std(MFPT_ebs_obs)))

print("orthosteric_binding_site -> binding_site_I " + str(np.mean(MFPT_obs_bsI)) + ' ' + str(np.std(MFPT_obs_bsI)))

