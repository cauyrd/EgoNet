#!/usr/bin/env python
# a bar plot with errorbars
import numpy as np
import matplotlib.pyplot as plt

N = 4 
egoMeans = (0.95, 0.84, 0.97, 0.82)
egoStd =   (0.1, 0.19, 0.03, 0.16)

ind = np.arange(N)  # the x locations for the groups
width = 0.35       # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind, egoMeans, width, color='r', yerr=egoStd)

chuangMeans = (0.92, 0.79, 0.9, 0.81)
chuangStd =   (0.08, 0.17, 0.09, 0.12)
rects2 = ax.bar(ind+width, chuangMeans, width, color='y', yerr=chuangStd)

# add some
ax.set_ylabel('Area Under ROC Curve (AUC)')
#ax.set_title('Scores by group and gender')
ax.set_xticks(ind+width)
ax.set_xticklabels( ('ego_linear', 'ego_nonlinear', 'no_ego_linear', 'no_ego_nonlinear') )
ax.set_yticks(np.arange(0,1.1,0.25))

ax.legend( (rects1[0], rects2[0]), ('EgoNet', 'Chuang et al.') )

def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%.2f'% height,
                ha='center', va='bottom')

autolabel(rects1)
autolabel(rects2)

plt.show()
