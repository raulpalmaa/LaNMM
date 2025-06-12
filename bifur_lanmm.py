############################################
#Auto07p script for the bifucation diagram of the LaNMM.
#Plots are handled with Python's matplotlib.
#run with 
#$: auto bifur_lanmm.py
############################################
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import rc
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap, Normalize
'''latex configuration
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})

rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})
rc('text', usetex=True)
'''
###############################################################
###############################################################
#################    Useful functions    ######################
###############################################################
###############################################################
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def pt_vals(f):
	return np.array([f[0][i]['PT'] for i in range(len(f[0]))])

def bifs(f,par):
	exceptions = ['No Label', 'RG', 'EP', 'UZ']  # List of exceptions
	return [[f[0][i]['TY name'],f[0][i][par]] for i in range(len(f[0])) if f[0][i]['TY name'] not in exceptions]

prop_cycle = plt.rcParams['axes.prop_cycle']
colores = prop_cycle.by_key()['color']

cmap = plt.get_cmap('Greys')

greys = truncate_colormap(cmap, 0.25, 0.85)
cblues = [(172/255,129/255,255/255),(95/255,14/255,255/255)]
creds  = [(255/255,148/255,115/255),(255/255,71/255,14/255)]
cblues = LinearSegmentedColormap.from_list('cblues',cblues,N = 50)
creds = LinearSegmentedColormap.from_list('creds',creds,N = 50)

###############################################################
###############################################################
##################    AUTO-07p part    ########################
###############################################################
###############################################################

init = run('lanmm',IPS=-2,NMX=100000,PAR={'p1':-100.0, 'p2' : 90.0})
ic = init(201)

fp=run(ic,IPS=1,NMX=400000,ISW=1,ICP=[1,2,7,8],UZSTOP={'p1' : 500.0})

hb1=fp('HB1')                                                                                   
lc1=run(hb1,IPS=2,ISP=2,ICP=[1,11,2,3,4,5,6],NMX=500000,ISW=1, UZSTOP={'p1' : 500.0,'PERIOD':100.0})#p2_0 UZSTOP={'PERIOD' : 100.0}) 

hb2=fp('HB2')                                                                                   
lc2=run(hb2,IPS=2,ISP=2,ICP=[1,11,2,3,4,5,6],NMX=500000,ISW=1, UZSTOP={'p1' : 500.0,'PERIOD':100.0})#p2_90 UZSTOP={'PERIOD' : 100.0})

###############################################################
###############################################################
##################    Plotting part    ########################
###############################################################
###############################################################

fig, axs = plt.subplots(2,1,figsize= (8,6))

axs[0].scatter(fp['p1'], fp['vp1'], c= (pt_vals(fp)< 0), cmap=greys, s=0.5)
axs[0].scatter(lc1['p1'],lc1['y1_min'],c=2 * (pt_vals(lc1) < 0), cmap=cblues, s=2)
axs[0].scatter(lc1['p1'],lc1['y1_max'],c=2 * (pt_vals(lc1) < 0), cmap=cblues, s=2)
axs[0].scatter(lc2['p1'],lc2['y1_min'],c=2 * (pt_vals(lc2) < 0), cmap=cblues, s=2)
axs[0].scatter(lc2['p1'],lc2['y1_max'],c=2 * (pt_vals(lc2) < 0), cmap=cblues, s=2)


axs[1].scatter(fp['p1'], fp['vp2'], c= (pt_vals(fp)< 0), cmap=greys, s=0.5)
axs[1].scatter(lc1['p1'],lc1['y2_min'],c=2 * (pt_vals(lc1) < 0), cmap=creds, s=2)
axs[1].scatter(lc1['p1'],lc1['y2_max'],c=2 * (pt_vals(lc1) < 0), cmap=creds, s=2)
axs[1].scatter(lc2['p1'],lc2['y2_min'],c=2 * (pt_vals(lc2) < 0), cmap=creds, s=2)
axs[1].scatter(lc2['p1'],lc2['y2_max'],c=2 * (pt_vals(lc2) < 0), cmap=creds, s=2)

axs[1].set_ylabel(r'$v_{P_1}$')
axs[1].set_ylabel(r'$v_{P_2}$')
axs[1].set_xlabel(r'$p_1$')

bfp = bifs(fp,'p1')
blc1 = bifs(lc1,'p1')
blc2 = bifs(lc2,'p1')

for ax in axs:
	for i,b in enumerate(bfp):
		ax.axvline(b[1],color='k', ls="--",alpha=0.7)
	for i,b in enumerate(blc1):
		ax.axvline(b[1],color='k', ls="--",alpha=0.7)
	for i,b in enumerate(blc2):
		ax.axvline(b[1],color='k', ls="--",alpha=0.7)
	ax.set_xlim(fp['p1'][0],fp['p1'][-1])

plt.show()
clean()
