import numpy as np
import glob
import pickle
import pyemma

#######################################Initialization##############################################

#Macrostates (consist of microstates) with sodium bound in different position 

extracellular_binding_site = [98, 405, 627, 30, 566]
orthosteric_binding_site = [56, 255, 589, 639, 125]
binding_site_I = [10, 13, 83, 318, 156]
binding_site_II = [77, 459, 316, 17, 588 ]
intracellular_binding_site = [397, 74, 494, 487, 59]


#initialization of empty array for mean free passage time between transition 
MFPT_ebs_obs = np.empty(200)
MFPT_obs_bsI = np.empty(200)
MFPT_bsI_bsII = np.empty(200)
MFPT_ibs_bsII = np.empty(200)

#####################################TPT calculations for each bootstrap sample####################

for j in range(200):
    files = pickle.load(open("./../CB1_bootstraping/bt_80_" + str(j) +"_files.pkl",'rb'))
    msm = pickle.load(open("./../CB1_bootstraping/bt_80_" + str(j) +"_msm.pkl",'rb'))
    
    tpt = pyemma.msm.tpt(msm, extracellular_binding_site, orthosteric_binding_site)
    MFPT_ebs_obs[j] = tpt.mfpt/10000   #extracellular_binding_site to orthosteric_binding_site

    tpt = pyemma.msm.tpt(msm, orthosteric_binding_site, binding_site_I)
    MFPT_obs_bsI[j] = tpt.mfpt/10000   #orthosteric_binding_site to binding_site_I

    tpt = pyemma.msm.tpt(msm, binding_site_I, binding_site_II)
    MFPT_bsI_bsII[j] = tpt.mfpt/10000  #binding_site_I to binding_site_II

    tpt = pyemma.msm.tpt(msm, intracellular_binding_site, binding_site_II)
    MFPT_ibs_bsII[j] = tpt.mfpt/10000  #intracellular_binding_site to binding_site_II


print('transition  mean  error \n')
print("extracellular_binding_site -> orthosteric_binding_site " + str(np.mean(MFPT_ebs_obs)) + ' ' + str(np.std(MFPT_ebs_obs)))

print("orthosteric_binding_site -> binding_site_I " + str(np.mean(MFPT_obs_bsI)) + ' ' + str(np.std(MFPT_obs_bsI)))

print("binding_site_I -> binding_site_II " + str(np.mean(MFPT_bsI_bsII)) + ' ' + str(np.std(MFPT_bsI_bsII)))

print("intracellular_binding_site -> binding_site_II " + str(np.mean(MFPT_ibs_bsII)) + ' ' + str(np.std(MFPT_ibs_bsII)))


