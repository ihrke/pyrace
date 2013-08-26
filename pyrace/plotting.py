import numpy as np
import pylab as pl
from itertools import cycle

plot_colors=['blue', 'red', 'green', 'yellow', 'magenta', 'cyan']
linestyles=['-', '--', '-.', ':']
linecycler = cycle(linestyles)
colorcycler=cycle(plot_colors)

def plot_bar_parameters(*pars, **kwargs):
    """
    plot different sets of parameters in a barplot

    The real signature (unsupported in Python 2):

    plot_bar_parameters(*pars, labels=None, yerr=None, title=None, legend=False, **kwargs):
    
    labels : list , label for each parameter set
    yerr   : list, yerr-array (lenght of npar) for each set; can be None for some
    """
    ax=pl.gca()
    labels=kwargs['labels'] if 'labels' in kwargs.keys() else None
    yerr=kwargs['yerr'] if 'yerr' in kwargs.keys() else None
    title=kwargs['title']    if 'title' in kwargs.keys() else None
    legend=kwargs['legend'] if 'legend' in kwargs.keys() else True
        
    nparsets=len(pars)
    npars=len(pars[0])
    cl=pars[0].__class__
    for par in pars:
        if cl!=par.__class__:
            raise ValueError('inconsistent list of parameters')

    if isinstance(labels,list):
        if len(labels)!=nparsets:
            raise ValueError('need one label per dataset')
    else:
        labels=['parameters %i'%i for i in range(nparsets)]

    if isinstance(yerr, list):
        if len(yerr)!=nparsets:
            raise ValueError('need one yerr per dataset')
    else:
        yerr=[None for i in range(nparsets)]

    colors=[plot_colors[i % len(plot_colors)] for i in range(nparsets)]
    pad=.1
    x=np.arange(npars)
    for ip,par in enumerate(pars):
        pl.bar(x+((1.-pad)/nparsets)*ip, par, yerr=yerr[ip], width=(1.-pad)/nparsets, label=labels[ip], color=colors[ip])
        
    ax.set_xticks(x+(1.-pad)/2.)
    ax.set_xticklabels(pars[0].parnames, rotation=0)
    if legend:
        pl.legend()
    if isinstance(title, str):
        pl.title(title)
        
