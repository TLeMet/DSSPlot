#!usr/bin/env python

# TLeMet March 28, 2019
# program for plotting secondary structure output of gmx do_dssp

'''usage : $python dssp_avg_plot.py .xvg '''

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

file = sys.argv[1]

content = np.transpose(np.array(pd.read_csv(file, skiprows=33, skipfooter=2, sep='\s+', engine='python'), dtype=float))

def CurveSmoothing(curve, degree=1):
    smoothed = curve*0
    for i in xrange(0,len(smoothed),1):
        if (i<degree):
            smoothed[i] = np.mean(curve[0:i+degree+1])
        elif (i>len(smoothed)-degree-1):
            smoothed[i] = np.mean(curve[i-degree:len(smoothed)])
        else:
            smoothed[i] = np.mean(curve[i-degree:i+degree+1])
    return smoothed 

columns = ['Coil','B-Bridge','Bend','Turn','A-Helix','5-Helix','3-Helix']
count = [0,0,0,0,0,0,0]

structure = np.transpose(content[2:])

for i in range(len(structure)):
    somme = sum(structure[i])
    for j in range(len(structure[i])):
        structure[i][j] = structure[i][j]/somme
    count[list(structure[i]).index(max(structure[i]))]+=1
count = [float(count[i])/sum(count) for i in range(len(count))]

averaged = np.concatenate((content[:1], np.transpose(structure)))

for i in range(1,len(averaged)):
	plt.plot(averaged[0], CurveSmoothing(averaged[i],100), marker=',', label=columns[i-1])
plt.title('DSSP plot for '+file[5:-4])
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.xlabel('time (ps)')
plt.ylabel('ratio of residues')

plt.savefig(file[:-4]+'.png', format='png', frameon=True, bbox_inches='tight')
plt.clf()

y_pos = np.arange(len(columns))+1
plt.bar(y_pos, count, align='center')
plt.title('Average DSSP plot for '+file[5:-4])
plt.xticks(y_pos, columns)
plt.ylim(0,1)
plt.ylabel('steps where structure is most found')

plt.savefig(file[:-4]+'_avg.png', format='png', frameon=True, bbox_inches='tight')
