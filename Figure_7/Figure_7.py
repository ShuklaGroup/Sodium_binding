import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc


#######################################Initialization##############################################
volume_CB1 = pickle.load(open('pocket_volume_CB1.pkl','rb'))
volume_CB2 = pickle.load(open('pocket_volume_CB2.pkl','rb'))



#######################################Figure Specification##############################################
hfont = {'fontname':'Helvetica','fontweight':'bold'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})   #Figure font definition

fig_wid = 10        #Width of the genarated figure
fig_hig = 7         #length of the genarated figure


#######################################Distribution plot##############################################

fig, axs = plt.subplots(1,1,figsize=(fig_wid,fig_hig))

sns.distplot(volume_CB1,ax=axs,hist=False,kde_kws={"color": "yellow", "lw": 4})
sns.distplot(volume_CB2,ax=axs,hist=False,kde_kws={"color": "purple", "lw": 4})

axs.set_xlim([0,int(int(max(volume_CB1))/5 + 1)*5])
axs.set_ylim([0,0.5])

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

axs.set_xlabel('Binding Pocket 2 volume (\AA$^3$)',**hfont,fontsize=24)
axs.set_ylabel('Probability density',**hfont,fontsize=24)

plt.savefig('Figure_7.png')
