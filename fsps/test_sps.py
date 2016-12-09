# -*- coding: utf-8 -*-

from __future__ import division, print_function

import pdb
import json
import numpy as np
from numpy.testing import assert_allclose
from fsps import StellarPopulation
#from .__init__ import githashes
from fsps import githashes
try:
    import h5py
except(ImportError):
    pass
try:
    import matplotlib.pyplot as pl
except(ImportError):
    pass

pop = StellarPopulation(zcontinuous=1, logzsol=0, imf_type=3, compute_vega_mags=True)
default_params = dict([(k, pop.params[k]) for k in pop.params.all_params])
libs = pop.libraries

def _reset_default_params():
    for k in pop.params.all_params:
        pop.params[k] = default_params[k]


def compare_reference(ref_in, current):
    """This method compares the values in a previously generated HDF5 file to
    the values in a supplied dictionary, and reports the maximum percentage
    change and where that change occurred.

    :param ref_in:
        h5py.File object with keys matching those in the given dictionary of values.

    :param current:
        Dictionary of calculated values that you want to compare to the old
        version.  Keys are the name of the quentities (e.g. maggies,
        stellar_mass, etc.)
    """
    ref_libs = json.loads(ref_in.attrs['libraries'])
    assert ref_libs == list(libs)

    diffstring = ("The maximum difference from the reference quantities "
                  "is {:4.5f} percent in the quantity {}")
    maxdiff, maxdq = 0, 'ANY'
    for k, v in current.items():
        diff = np.atleast_1d(v /ref_in[k]  - 1)  # percent difference
        maxdind = np.argmax(np.abs(diff))
        if diff[maxdind] > maxdiff:
            maxdiff = diff[maxdind]
            maxdq = k

    return diffstring.format(maxdiff, maxdq)



'''        
def test_dl(make_plots = True, reference_in=None, reference_out=None):
    """This function reproduces Conroy & Gunn (2009) Figure 3. These plots
        evaluate the effect of changing \Delta L, the shift in log(L_bol) of the TP-AGB phase
    
        (Fig 4) \Delta T, the shift in log(T_eff) of the TP-AGB phase
        (Fig 5) f_bhb, the fraction of blue horizontal branch stars
        (Fig 6) S_BS, the specific frequency of blue stragglers
        (Fig 7) m_c, the characteristic mass of the IMF

    :make_plots: (optional)
        Accepts a boolean.  Says whether or not you want to make plots or just do quantitative comparison.

    :param reference_in: (optional)
        The name (and path) to a file to use as a refernce for the quantities
        being calculated.  The currently calculated quantities will be compared
        to the ones in this file.

    :param reference_out: (optional)
        Name of hdf5 file to output containing the calculated quantities

    :returns figs:
        List of matplotlib.pyplot figure objects.
    """
    

    _reset_default_params()

    pop.params['tpagb_norm_type'] = 1

    figure_data = {}
    

    for val in [-0.4,-0.2,0.0,0.2,0.4]:
        figure_data[str(val)] = {}
        pop.params['dell'] = val
        vband = pop.get_mags(bands=['v']).flatten()
        rband = pop.get_mags(bands=['sdss_r']).flatten()
        uband = pop.get_mags(bands=['u']).flatten()
        kband = pop.get_mags(bands=['wfcam_k']).flatten()
        times = 1e-9*10**pop.log_age.flatten()
        mass = pop.stellar_mass.flatten()

        l_v = 10**(-0.4*(vband - 4.79)) #luminosity in units of L_sun. 4.79 is v-band mag of Sun
        logl_v = np.log10(l_v)
    
        figure_data[str(val)] = {'log_mL': np.log10(mass) - logl_v,
                            'VRcolor': vband - rband,
                            'UVcolor': uband - vband,
                            'VKcolor': vband - kband,
                            'times': times}

                        
    times = 1e-9*10**pop.log_age.flatten()
    from astropy.table import Table
    t = Table.read('fig3_test3.hdf5')
    val = [-0.4,-0.2,0.0,0.2,0.4]
    for i in range(len(val)):
        figure_data[str(val[i])]={'log_mL': np.log10(t['m2l'][i]), 
                                  'VRcolor':t['VR'][i], 
                                  'UVcolor':t['UV'][i],
                                  'VKcolor':t['VK'][i], 
                                  'times': times}
           
    #Call the plotting function
    #if make_plots:                 
    #    [fig] = make_plots_3to7(figure_data,r'$\Delta_L$')

    # Make quantitative comparisons to an old set of values stored in an hdf5 file
    if reference_in is not None:
        full_ref = h5py.File(reference_in, 'r')
        for val in np.sort(figure_data.keys()):
            current = figure_data[val]
            ref = full_ref[val]
            diffstring = compare_reference(ref, current)
            print(r'Comparing to reference for parameter $\Delta_L$ = '+str(val)+':')
            print(diffstring)
        full_ref.close()        


    # Here we write out an hdf5 file
    if reference_out is not None:
        ref = h5py.File(reference_out, 'w')
        for key, val in figure_data.items():
            grp = ref.create_group(key)
            grp.attrs['libraries'] = json.dumps(libs)
            grp.attrs['default_params'] = json.dumps(default_params)
            #ref.attrs['git_history'] = json.dumps(githashes)
            for k2, v2 in val.items():
                grp.create_dataset(k2, data=v2)
            
        # now we add useful version things

        ref.close()
    
    if make_plots:
        return [fig]
    else:
        return
'''        


def dotests(which_tests = [True, True, True, True, True], ref_in = [None,None,None,None,None], ref_out = [None,None,None,None,None], mk_plots = [True, True, True, True, True]):
    fig3_params = {'var_name':'dell', 'plot_name': r'$\Delta_L$', 'var_values': [-0.4,-0.2,0.0,0.2,0.4]}
    fig4_params = {'var_name':'delt', 'plot_name': r'$\Delta_T$', 'var_values': [-0.2,-0.1,0.0,0.1,0.2]}
    fig5_params = {'var_name':'fbhb', 'plot_name': r'$f_{bhb}$', 'var_values': [0.0,0.2,0.4,0.6,0.8,1.0]}
    fig6_params = {'var_name':'sbss', 'plot_name': r'$S_{bss}$', 'var_values': [0.0,1.0,2.0,3.0,4.0,5.0]}
    fig7_params = {'var_name':'vdmc', 'plot_name': r'$m_c$', 'var_values': [0.08,0.28,0.48,0.68,1.0,2.0]}
    
    fig_params = [fig3_params,fig4_params,fig5_params,fig6_params,fig7_params]
    
    for i in range(len(which_tests)):
        if which_tests[i]:
            test_fig37(fig_params[i],reference_in = ref_in[i], reference_out = ref_out[i], make_plots = mk_plots[i])
    return

        
def test_fig37(fig_params, make_plots = True, reference_in=None, reference_out=None):
    """This function reproduces Conroy & Gunn (2009) Figures 3-7. These plots
        evaluate the effect of changing:
        (Fig 3) \Delta L, the shift in log(L_bol) of the TP-AGB phase
        (Fig 4) \Delta T, the shift in log(T_eff) of the TP-AGB phase
        (Fig 5) f_bhb, the fraction of blue horizontal branch stars
        (Fig 6) S_BS, the specific frequency of blue stragglers
        (Fig 7) m_c, the characteristic mass of the IMF

    :make_plots: (optional)
        Accepts a boolean.  Says whether or not you want to make plots or just do quantitative comparison.

    :param reference_in: (optional)
        The name (and path) to a file to use as a refernce for the quantities
        being calculated.  The currently calculated quantities will be compared
        to the ones in this file.

    :param reference_out: (optional)
        Name of hdf5 file to output containing the calculated quantities

    :returns figs:
        List of matplotlib.pyplot figure objects.
    """
    
    var_name = fig_params['var_name']
    plot_name = fig_params['plot_name']
    var_values = fig_params['var_values']
    

    _reset_default_params()

    pop.params['tpagb_norm_type'] = 1

    figure_data = {}
    

    for val in var_values:
        figure_data[str(val)] = {}
        pop.params[var_name] = val
        vband = pop.get_mags(bands=['v']).flatten()
        rband = pop.get_mags(bands=['cfht_r']).flatten()
        uband = pop.get_mags(bands=['u']).flatten()
        kband = pop.get_mags(bands=['wfcam_k']).flatten()
        times = 1e-9*10**pop.log_age.flatten()
        mass = pop.stellar_mass.flatten()

        l_v = 10**(-0.4*(vband - 4.79)) #luminosity in units of L_sun. 4.79 is v-band mag of Sun
        logl_v = np.log10(l_v)
    
        figure_data[str(val)] = {'log_mL': np.log10(mass) - logl_v,
                            'VRcolor': vband - rband,
                            'UVcolor': uband - vband,
                            'VKcolor': vband - kband,
                            'times': times}
    
                            
           
    #Call the plotting function
    if make_plots:                 
        [fig] = make_plots_3to7(figure_data,plot_name,var_name)

    
    # Make quantitative comparisons to an old set of values stored in an hdf5 file
    if reference_in is not None:
        full_ref = h5py.File(reference_in, 'r')
        for val in np.sort(figure_data.keys()):
            current = figure_data[val]
            ref = full_ref[val]
            diffstring = compare_reference(ref, current)
            print(r'Comparing to reference for parameter '+plot_name+' = '+str(val)+':')
            print(diffstring)
        full_ref.close()        


    # Here we write out an hdf5 file
    if reference_out is not None:
        ref = h5py.File(reference_out, 'w')
        for key, val in figure_data.items():
            grp = ref.create_group(key)
            grp.attrs['libraries'] = json.dumps(libs)
            grp.attrs['default_params'] = json.dumps(default_params)
            ref.attrs['git_history'] = json.dumps(githashes)
            for k2, v2 in val.items():
                grp.create_dataset(k2, data=v2)    
    
    if make_plots:
        return [fig]
    else:
        return
        



def make_plots_3to7(dict_in,param_varied,var_name):
    from matplotlib.ticker import MultipleLocator
    import matplotlib.cm as cm
    
    fig, ((ax1,ax2),(ax3,ax4)) = pl.subplots(2,2, figsize=[12,12], sharex=True)

    dashes = [(None, None), (2,2), (5,5), (5,3,1,3), (5,3,1,3,1,3), (10,4)]
    colors = cm.gnuplot2(np.union1d([0], np.linspace(0.3,0.8,len(dict_in.keys())-1)))

    i = 0
    ks = [float(val) for val in dict_in.keys()]
    ks = [str(val) for val in np.sort(ks)]
    for val in ks:
        log_mL = dict_in[val]['log_mL']
        VRcolor = dict_in[val]['VRcolor']
        UVcolor = dict_in[val]['UVcolor']
        VKcolor = dict_in[val]['VKcolor']
        times = dict_in[val]['times']
    
        ax1.plot(times, log_mL, marker = '', dashes = dashes[i], color=colors[i], lw=2, label=param_varied+' = '+str(val))
        ax2.plot(times, UVcolor, marker = '', dashes = dashes[i], color=colors[i], lw=2, label=param_varied+' = '+str(val))
        ax3.plot(times, VRcolor, marker = '', dashes = dashes[i], color=colors[i], lw=2, label=param_varied+' = '+str(val))
        ax4.plot(times, VKcolor, marker = '', dashes = dashes[i], color=colors[i], lw=2, label=param_varied+' = '+str(val))
        
        i += 1

    # Make figures look nice
    # Upper Left: M/L vs time
    ax1.set_xlim([0,15])
    ax1.set_ylim([-0.3,1.0])
    ax1.set_xticks([0,5,10,15])
    minor_locator = MultipleLocator(1)
    ax1.xaxis.set(minor_locator=minor_locator)
    ax1.set_yticks(np.arange(-0.2, 1.1, 0.2))
    minor_locator = MultipleLocator(0.2/4)
    ax1.yaxis.set(minor_locator=minor_locator)
    ax1.set_ylabel(r'log(M/$\rm L_V$)')
    ax1.legend(loc='lower right')

    # Upper Right: U - V vs time
    ax2.set_ylim([0.6,1.8])
    ax2.set_yticks(np.arange(0.6, 1.81, 0.2))
    minor_locator = MultipleLocator(0.2/4)
    ax2.yaxis.set(minor_locator=minor_locator)
    ax2.set_ylabel(r'U - V')

    # Lower Left: V - R vs time
    ax3.set_ylim(0.4, 0.7)
    ax3.set_yticks(np.arange(0.4, 0.71, 0.05))
    minor_locator = MultipleLocator(0.05/4)
    ax3.yaxis.set(minor_locator=minor_locator)
    ax3.set_ylabel(r'V - R')
    ax3.set_xlabel(r'Time (Gyr)')

    # Lower Right: V - K vs time
    ax4.set_ylim([2.0,4.0])
    ax4.set_yticks(np.arange(2, 4.1, 0.5))
    minor_locator = MultipleLocator(0.5/4)
    ax4.yaxis.set(minor_locator=minor_locator)
    ax4.set_xlabel('Time (Gyr)')
    ax4.set_ylabel('V - K')
    
    fig.savefig('test_%s.png' %var_name,bbox = 'tight')
    #pl.show()
    pl.close()
    return [fig]



#test_dl(reference_in='dltest.hdf5',make_plots = False)
dotests(which_tests = [True, True, True, True, True], ref_in = [None,None,None,None,None], ref_out = ['t3.hdf5','t4.hdf5','t5.hdf5','t6.hdf5','t7.hdf5'], mk_plots = [True, True, True, True, True])

