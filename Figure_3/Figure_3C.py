import numpy as np
import matplotlib.pyplot as plt
import pickle
import matplotlib as mpl
from matplotlib import rc


#######################################Initialization##############################################
fraction_CB1 = np.empty(200)
fraction_CB2 = np.empty(200)

#####################################Fraction calculations for each bootstrap sample####################
for j in range(200):        #CB1
    CB1_fractions = pickle.load(open("../CB1_bootstraping/bt_80_" + str(j) +"_Nterm_fraction.pkl",'rb'))
    msm = pickle.load(open("../CB1_bootstraping/bt_80_" + str(j) +"_msm.pkl",'rb'))
    CB1_fractions_con = np.concatenate(CB1_fractions)
    weights = np.concatenate(msm.trajectory_weights())
    fraction_CB1[j] = np.dot(CB1_fractions_con,weights)

for j in range(200):        #CB2
    CB2_fractions = pickle.load(open("../CB2_bootstraping/bt_80_" + str(j) +"_Nterm_fraction.pkl",'rb'))
    msm = pickle.load(open("../CB2_bootstraping/bt_80_" + str(j) +"_msm.pkl",'rb'))
    CB2_fractions_con = np.concatenate(CB2_fractions)
    weights = np.concatenate(msm.trajectory_weights())
    fraction_CB2[j] = np.dot(CB2_fractions_con,weights)

 
fraction = [fraction_CB1,fraction_CB2]

#######################################Figure Specification##############################################
hfont = {'fontname':'Helvetica','fontweight':'bold'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})   #Figure font definition

fig_wid = 10        #Width of the genarated figure
fig_hig = 7         #length of the genarated figure


#######################################Box plot##############################################
fig, axs = plt.subplots(1,1,figsize=(fig_wid,fig_hig)) 

axs.boxplot(fraction)       #box plot 

plt.xticks(fontsize=22)     #x-axis tick size                                         
plt.yticks(fontsize=22)     #x-axis tick size

plt.xticks([1, 2], ['CB$_1$ (D104$^{N-term}$)', 'CB$_2$ (D24$^{N-term}$)'])                 #x-axis label
plt.ylabel('Sodium bound population (D$^{ECL1}$)',**hfont,fontsize=30)                      #y-axis label 

plt.savefig('Figure_3C.png')

