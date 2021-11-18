import glob
import numpy as np
import os
import pyemma
import math
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc


#######################################Feature matrix loading##############################################

sodium_z_dist = pickle.load(open('./CB1_sodium_ion_z_dist.pkl','rb'))

orthosteric_residue_movement = pickle.load(open('./CB1_orthosteric_residue_movement.pkl','rb'))

#######################################MSM loading##############################################

msm = pickle.load(open("./../CB1_final_MSM_obj.pkl","rb")) #loading of MSM object

weights=np.concatenate(msm.trajectory_weights())                                           #weights(probability density of each frames)

#######################################Parameters and hyperparameters definition##############################################

x_bins = 100    #number of bins used to divide the landscape in x-direction 
y_bins = 20     #number of bins used to divide the landscape in y-direction 

R = 0.001987    #Boltzman's constant (kcal/mol/K)
T = 300         #Temperature (K)

#######################################Figure Specification##############################################
hfont = {'fontname':'Helvetica','fontweight':'bold'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})   #Figure font definition

fig_wid = 10        #Width of the genarated figure
fig_hig = 7         #length of the genarated figure
cmap = mpl.cm.jet   #color bar used in the figure

Max_energy =6       #maximum energy projected in color bar

#######################################2-D Histogram##############################################

x_data =  np.concatenate(sodium_z_dist)*10                     #assign x_data as 1-D array with x-direction feature (change nm to angtrom)
y_data =  np.concatenate(orthosteric_residue_movement)*10      #assign y_data as 1-D array with y-direction feature (change nm to angtrom)




x_data_min =  np.min(x_data)                        #minimum value of x-direction feature
y_data_min =  np.min(y_data)                        #minimum value of y-direction feature
x_data_max =  np.max(x_data)                        #maximum value of x-direction feature
y_data_max =  np.max(y_data)                        #maximum value of y-direction feature

x_hist_lim_low =  x_data_min -0.5                   #mimimum limit of histogram in x-direction 
y_hist_lim_low =  y_data_min -0.5                   #mimimum limit of histogram in y-direction
x_hist_lim_high = x_data_max +0.5                   #maximum limit of histogram in x-direction
y_hist_lim_high = y_data_max  +0.5                  #maximum limit of histogram in y-direction

hist= np.histogram2d(x_data,y_data, bins=[x_bins,y_bins],
                     range = [[x_hist_lim_low,x_hist_lim_high],[y_hist_lim_low,y_hist_lim_high]],
                     density= True,weights=weights) #2-D histogram 

prob_density = hist[0]                              #probablity density obtained from histogram
xedge = hist[1]                                     #edges of the bins obtained from histogram in x-direction 
yedge = hist[2]                                     #edges of the bins obtained from histogram in y-direction 

x_bin_size = xedge[1]-xedge[0]                      #x-bin size
y_bin_size = yedge[1]-yedge[0]                      #y-bin size

#######################################Free energy calculations##############################################

free_energy = -R*T*np.log(prob_density*x_bin_size*y_bin_size)     #absolute value of the free energy in each bin 
min_free_energy= np.min(free_energy)                              #minimum value of the free energy
delta_free_energy = free_energy - min_free_energy                 #Relative free energy in each bin 


#######################################Contour plot##############################################
fig, axs = plt.subplots(1,1,figsize=(fig_wid,fig_hig))            

xx = [(xedge[i]+xedge[i+1])/2 for i in range(len(xedge)-1)]                                     #average values of x-bins
yy = [(yedge[i]+yedge[i+1])/2 for i in range(len(yedge)-1)]                                     #average values of y-bins

cd =axs.contourf(xx,yy,delta_free_energy.T,
                 np.linspace(0,Max_energy,Max_energy*5+1), vmin=0.0, vmax=Max_energy,cmap=cmap) #contour plot

cbar = fig.colorbar(cd,ticks=range(Max_energy+1))                                               #color bar
cbar.ax.set_yticklabels(range(Max_energy+1),fontsize=22)                                        #ticklabels of color bar 
cbar.ax.set_ylabel('Free Energy (Kcal/mol)', labelpad=15,**hfont,fontsize=30)                   #axis labels of color bar 

axs.set_xlim([-10,30])                                                                          #min and max limit of x-axis
axs.set_ylim([3,13])                                                                            #min and max limit of y-axis

axs.set_xticks(range(int(-10),int(30)+1,10))                                                    #x-axis ticks
axs.set_xticklabels(range(int(-10),int(60)+1,10))                                               #y-axis ticks

axs.set_yticks(range(int(3),int(13)+1,2))                                                       #x-axis tick size
axs.set_yticklabels(range(int(3),int(13)+1,2))                                                  #y-axis tick size

plt.xlabel('z' + ' (\AA)', **hfont,fontsize=30)                                                 #x-axis label
plt.ylabel('F170$^{2.57}$-V196$^{3.32}$ distance' + ' (\AA)', **hfont,fontsize=30)               #y-axis label 

plt.xticks(fontsize=22)
plt.yticks(fontsize=22)

plt.tight_layout()

Filename = 'Figure_4A'
plt.savefig(Filename,transparent=True,dpi =700)
plt.close()

