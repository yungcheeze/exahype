import matplotlib.pylab
import matplotlib.pyplot as plt

"""
.. module:: scalingplot
  :platform: Unix, Windows
  :synopsis: Provides functions to create scaling plots.
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Provides functions to create scaling plots.
"""

def plot_scaling(ax,x,y,label_,color_="blue",marker_="s",markerfacecolor_="blue",hyperthreading=False,annotate=False):
    """
    Plots y over x. If 'hyperthreading' is selected, the last value is plotted 
    separately. If 'annotate' is selected, the plot is annotated with the 'y'
    values.
    """
    x2 = x
    
    # In the hyperthreading case, we only want to plot a line
    # for the real threads. The hyperthreading value is plotted separately
    if not hyperthreading: 
        plt.plot(x,y,label=r"%s" % label_,markersize=8,marker=marker_,markevery=1,lw=1.2,color=color_,markerfacecolor=color_,markeredgecolor=color_) 
    else:
        x2 = x[0:-1] # [inclusive:exclusive]
        x2.append(x2[-1]+2)
        
        plt.plot(x2[0:-1],y[0:-1],label=r"%s" % label_,markersize=8,marker=marker_,markevery=1,lw=1.1,color=color_,markerfacecolor=markerfacecolor_,markeredgecolor=color_) 
        plt.plot(x2[-1],y[-1],label="",markersize=8,marker="s",markevery=1,lw=1.1,color=color_,markerfacecolor=markerfacecolor_,markeredgecolor=color_) 
    
    if annotate:
        for i,j in zip(x2,y):
            ax.annotate("%.2f" % j,xy=(i,j+0.5),color="blue")
    return
