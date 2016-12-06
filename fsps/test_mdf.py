from __future__ import division
import numpy as np
import fsps
try:
	import h5py
except(ImportError):
	pass
try:
	import matplotlib.pyplot as plt
except(ImportError):
	pass

################################################################################
# Functions to produce a diagnositic plot similar to Figure 9 of Conroy+09     #
#                                                                              #
# usage:                                                                       #
#                                                                              #
# >>> tests_sps_mdf()                                                          #
#                                                                              #
# will produces a plot with default parameters similar to Conroy+09.           #
# Note that the figure is not exactly the same, because a different form of    #
# of the MDF is used, and different stellar models are used.                   #
# if provided with reference .h5 file to compare, then specify h5_in in params #
# otherwise, if you want to generate a reference .h5 file, specify h5_out      #
# e.g.                                                                         #
#                                                                              #
# >>> test_sps_mdf(h5_in='old_reference.h5',h5_out='new_reference.h5')         #
#                                                                              #
# Kareem El-Badry and David Khatami, 12/01/2016                                #
################################################################################

pop = fsps.StellarPopulation(compute_vega_mags = True, zcontinuous=2,
            sfh=0, logzsol = logzsol)
default_params = dict([(k,pop.params[k]) for k in pop.params.all_params])
libs = pop.libraries

def _reset_default_params():
	for k in pop.params.all_params:
		pop.params[k] = default_params[k]


#specify the test run parameters below
logzsol = -0.75
pmetals = np.array((0.5,1.0,2.0,5.0,-99.0))
ages = np.array((3.2,12.6))

def test_sps_mdf(make_plots=True,h5_in=None,h5_out=None):
    '''
    Perform a unit test similar to Figure 9 in Conroy, Gunn and White (2009).
    For given metallicity distribution functions (MDF) specified by pmetals and ages,
    plot against single-metallicity SSP runs for 3 different colors -- (V-K),(U-V), (V-R).
    If a reference file is specified, compare the run and report maximum percent differences
    of color+metallicity output. If output path is specified, write the test run to file.
    
    parameters
    ----------
    make_plots: bool. set to True to generate diagnostic plots
    h5_in: string. path to the input reference .h5 file. If None, then
                   don't perform comparison
    h5_out: string. path to write the current test run to an .h5 file. If None,
                    then skips the write.
    '''
    _reset_default_params()

    mean_logZs, ssp_logzsols, mdf_colors, ssp_colors = get_sps_colors()
    if make_plots:
        plot_diagnostic_metallicity_figure(mean_logZs, ssp_logzsols,mdf_colors,ssp_colors)
        
    V_K, U_V, V_R = mdf_colors
        
    test_dataset = {'ages' : ages,
                    'pmetals' : pmetals,
                    'logzsol' : logzsol,
                    'ssp_logz' : ssp_logzsols,
                    'ssp_colors' : ssp_colors,
                    'mdf_logz' : mean_logZs,
                    'V_K' : V_K,
                    'U_V' : U_V,
                    'V_R' : V_R}
    
    #perform the comparison if provided with reference .h5 file
    if h5_in is not None:
        ref = h5py.File(h5_in,'r')
        diffstrings = compare_sps_mdf(ref,test_dataset)
        print(diffstrings[0])
        print(diffstrings[1])
        ref.close()
    
    #write the test to file if the h5_out path is specified
    if h5_out is not None:
        ref = h5py.File(h5_out,'w')
        for k, v in test_dataset.items():
            d = ref.create_dataset(k,data=v)
        ref.attrs['default_params'] = json.dumps(default_params)
        ref.attrs['libraries'] = json.dumps(libs)
        ref.attrs['git_history'] = json.dumps(githashes)
        ref.close()
        
    return test_dataset
        

def compare_sps_mdf(ref,test):
    '''
    compares a test MDF+SSP run to a reference file, reporting
    which logZ, age, MDF, and color give the largest percent difference
    
    parameters
    ----------
    ref: h5 file. contains the reference datsets to compare against
    test: dictionary. holds the output of the test run.
    
    
    '''
    
    #here we make sure that the input parameters are the same
    #otherwise the comparison isn't veyr meaningful
    assert np.all(np.equal(ref['ages'],test['ages'])), "Reference ages %s differ from test ages %s" %(ref['ages'],test['ages'])
    assert np.all(np.equal(ref['logzsol'],test['logzsol'])), "Reference logzsol %s differs from test logzsol %s" \
                                                % (ref['logzsol'],test['logzsol'])
    assert np.all(np.equal(ref['pmetals'],test['pmetals'])), "Reference mdf arr %s differs from test mdf arr %s" \
                                                % (ref['pmetals'],test['pmetals'])
        
    #differences in mean logZ
    logz_diff = np.abs(test['mdf_logz']/ref['mdf_logz']-1.)
    logz_max_ind = np.argmax(logz_diff)
    logz_max_diff = logz_diff[logz_max_ind]
    logz_max_pmetals = test['pmetals'][logz_max_ind]
    
    #differences in MDF colors
    max_col = 'ANY'
    max_col_diff = -1.
    max_col_age = -1.
    max_pmetal = -np.inf
    
    for col in ['V_K','U_V','V_R']:
        col_diff = np.reshape(np.abs(np.array(test[col])/np.array(ref[col])-1.),-1)
        max_diff_ind = np.argmax(col_diff)
        max_diff = col_diff[max_diff_ind]
        if max_diff > max_col_diff:
            max_col = col
            max_col_diff = col_diff[max_diff_ind]
            max_pmetal = ref['pmetals'][max_diff_ind%np.size(pmetals)]
            max_col_age = ref['ages'][max_diff_ind/np.size(pmetals)]
        
        
    max_logz_string = ("Max difference in LogZ is %4.5f percent for pmetal of %.2f" %
                       (logz_max_diff*100,logz_max_pmetals))
    max_col_string = ("Max difference in color is %4.5f percent for %s with age %.2f Gyr and pmetal %.2f" %
                      (max_col_diff*100,max_col,max_col_age,max_pmetal))
            
    return max_logz_string, max_col_string
    
        
        
def get_sps_colors():
    '''
    computes U-V, V-K, and V-R colors of both stellar populations with different 
    MDFs and SSPs, all for two ages. 
    
    parameters
    ----------
    pop: fsps StellarPopulation object.
    logzsol: float. this sets the mode of the MDF. 
    pmetals: list of floats. these set the shape of the mdf, with higher 
        positive p values producing a narrower mdf. If negative, uses a simple 
        triangular mdf. 
    ages: list of 2 floats. Age at which to compute colors, in Gyr.  
    '''
    bands = ['u', 'cousins_r', 'v', 'wfcam_k']
    V_Ks, U_Vs, V_Rs, mean_logZs = [], [], [], []
    pop._zcontinuous = 2
    for i, age in enumerate(ages):
        this_age_VK, this_age_UV, this_age_VR = [], [], []
        for j, p in enumerate(pmetals):
            pop.params['pmetals'] = p
            U, R, V, K = pop.get_mags(tage = age, bands = bands)
            V_K, U_V, V_R = V - K, U - V, V - R
            this_age_VK.append(V_K); this_age_UV.append(U_V); this_age_VR.append(V_R)
            if not i:
                mean_logZ = get_mean_logZ_of_mdf(pmetals = p)
                mean_logZs.append(mean_logZ)
        V_Ks.append(this_age_VK); U_Vs.append(this_age_UV); V_Rs.append(this_age_VR)
        
    # Now the SSPs
    pop._zcontinuous = 1
    ssp_logzsols = np.linspace(-2., 0.3, 20)
    V_Ks_ssp, U_Vs_ssp, V_Rs_ssp = [], [], []
    for i, age in enumerate(ages):
        this_age_VK, this_age_UV, this_age_VR = [], [], []
        for j, logz in enumerate(ssp_logzsols):
            pop.params['logzsol'] = logz
            U, R, V, K = pop.get_mags(tage = age, bands = bands)
            V_K, U_V, V_R = V - K, U - V, V - R
            this_age_VK.append(V_K); this_age_UV.append(U_V); this_age_VR.append(V_R)
        V_Ks_ssp.append(this_age_VK); U_Vs_ssp.append(this_age_UV); V_Rs_ssp.append(this_age_VR)
    
    mdf_colors, ssp_colors = [V_Ks, U_Vs, V_Rs], [V_Ks_ssp, U_Vs_ssp, V_Rs_ssp]
    
    return np.array(mean_logZs), np.array(ssp_logzsols), mdf_colors, ssp_colors




def plot_diagnostic_metallicity_figure(mean_logZs, ssp_logzsols, mdf_colors, ssp_colors):
    '''
    This creates a diagnostic plot similar to figure 9 in Conroy, Gunn and White (2009).
    
    The figure is not identical to Figure 9, because the form of the MDF currently 
    supported by python fsps is different from the form used in the paper, and because
    python fsps defaults to the MIST stellar models rather than Padova. 
    
    parameters
    -----------
    mean_logZs: list of floats. array of mean logZ for the MDFs
    ssp_logzsols: 2 lists of floats. one for each age (2). gives logZ of
                  various single metallicity runs.
    mdf_colors: 2x3xN list of floats. one for each (3) color (V-K, U-V, V-R) and age (2).
                gives colors of each MDF run at two different ages.
    ssp_colors: 2x3xN list of floats. one for each (3) color (V-K,U-V,V-R) and age(2).
                gives colors of each single-metallicity SSP at two different ages.
    '''

    V_Ks, U_Vs, V_Rs = mdf_colors
    V_Ks_ssp, U_Vs_ssp, V_Rs_ssp = ssp_colors
    ylabels = [[r'${\rm d}n/{\rm d}\log(Z)$', 'U - V'], ['V-R', 'V-K']]
    fig, ax = plt.subplots(2, 2, figsize = (11, 8))
    
    colors = 'red','blue','green','brown','grey'
    
    for i in range(2):
        for j in range(2):
            ax[i, j].xaxis.set_tick_params(labelsize = 18)
            ax[i, j].yaxis.set_tick_params(labelsize = 18)
            ax[i, j].set_xlim(-2, 0.5)
            ax[i, j].set_ylabel(ylabels[i][j], fontsize = 18)
            if i == 0:
                ax[i, j].set_xticklabels([])
            else:
                ax[i, j].set_xlabel(r'$\log(Z/Z_{\odot})$', fontsize = 18)
            if (i == 0) and (j == 0):
                # show the mdfs
                logZ = np.linspace(-2, 1.5, 500)
                for k,p in enumerate(pmetals):
                    log_Z_over_Z_sun, normed_mdf = normalized_mdf(logZ = logZ,pmetals=p)
                    ax[i, j].plot(log_Z_over_Z_sun, normed_mdf, label = r'$p = %.1f$' % p,color=colors[k])
                ax[i, j].legend(loc = 'best', frameon = False, fontsize = 12)
            else:
                # plot the colors
                if i == 0 and j == 1:
                    color, ssp_color = U_Vs, U_Vs_ssp
                elif i == 1 and j == 0:
                    color, ssp_color = V_Rs, V_Rs_ssp
                elif i == 1 and j == 1:
                    color, ssp_color = V_Ks, V_Ks_ssp
                for k in range(5):
                    ax[i, j].plot(mean_logZs[k], color[0][k], 's', markerfacecolor=colors[k])
                    ax[i, j].plot(mean_logZs[k], color[1][k], 'ko',markerfacecolor=colors[k])
                ax[i, j].plot(ssp_logzsols, ssp_color[0], 'k--')
                ax[i, j].plot(ssp_logzsols, ssp_color[1], 'k-')
            if (i == 0) and (j == 1):
                ax[i, j].plot([], [], 'k--', label = 't = %.1f Gyr, ssp' %ages[0])
                ax[i, j].plot([], [], 'k-', label = 't = %.1f Gyr, ssp' %ages[1])
                ax[i, j].plot([], [], 's', markerfacecolor='None', label = 't = %.1f Gyr, mdf' %ages[0])
                ax[i, j].plot([], [], 'ko', label = 't = %.1f Gyr, mdf' %ages[1])
                ax[i, j].legend(loc = 'upper left', fontsize = 12, frameon = False)

    plt.subplots_adjust(wspace = 0.4, hspace = 0.1)

def not_normalized_mdf(logZ,pmetals):
    '''
    the form of the mdf implemented in fsps if zcontinuous == 2
    and pmetals > 0:
    
    dn/dlogZ = (Z*exp(-Z))**pmetals
    
    This does *not* correspond to a closed box or any other mdf
    predicted by analytical models. 
        
    parameters
    -------------
    logZ: array or float; this is log10(Z/Z0), where Z0 is 
        effectively the mode of the MDF. 
    pmetals: float. 
    '''
    mdf = (10**logZ * np.exp(-10**logZ))**pmetals
    return mdf
    
def get_mdf_triangular_kernel(zmets = 'MIST'):
    '''
    The triangular smoothing kernel implemented in zinterp.f90
    
    The form and width of the kernel depends on the number of 
    different metallicities ('zmets'). For the MIST models, fsps
    currently only has 12 zmets, so the kernel is fairly wide. 
    Increasing the number of zmets would decrease the width
    of the kernel. 
    
    fsps normalizes the kernel so that the weights sum to 1. This
    normalizes it so that it integrates 1. 
    
    parameters
    ------------
    zmets: array of logzsols for which stellar models are available; i.e.
        zmets = log10(pop.zlegend/Z_sun). If set to "MIST", uses the default
        set of isochrones from fsps. 
    '''
    if zmets == 'MIST':
        zmets = np.array([-2.50, -2.00, -1.75, -1.50, -1.25, 
                           -1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50])
    nz = len(zmets)
    w1, w2, w3 = 0.25, 0.5, 0.25
    nearest_ind = np.where(zmets == np.max(zmets[zmets < logzsol]))[0][0]
    zlo = max(min(nearest_ind, nz-1), 0)
    dz = (logzsol - zmets[zlo])/(zmets[zlo + 1] - zmets[zlo])

    mdf = np.zeros(len(zmets))
    mdf[max(zlo-1, 0)] = w1*(1-dz)
    mdf[min(zlo+2, nz-1)] = w3*dz
    mdf[zlo] += w2*(1-dz) + w1*dz
    mdf[zlo+1] += w3*(1-dz) + w2*dz
    
    #renormalize
    mdf /= np.trapz(mdf, x = zmets)
    
    return zmets, mdf
    
def normalized_mdf(logZ, pmetals = 2):
    '''
    Form of the mdf implemented in fsps if zcontinuous = 2
    
    This normalizes the mdf by numerically integrating. 
    
    parameters
    ------------
    logZ: array or float; this is log10(Z/Z0), where Z0 is 
        effectively the mode of the MDF. 
    pmetals: float. 
    '''
    log_Z_over_Z_sun = logZ + logzsol 

    if pmetals < 0:
        zmets, mdf = get_mdf_triangular_kernel()
        interped_mdf = np.interp(log_Z_over_Z_sun, zmets, mdf)
        return log_Z_over_Z_sun, interped_mdf
        
    else:
        # normalization gets wonky outside this range.
        assert (pmetals >= 0.1) and (pmetals <= 10)
        integrate_logz = np.linspace(-10, 5, 10000)
        mdfs = not_normalized_mdf(integrate_logz,pmetals=pmetals)
        intgrl = np.trapz(mdfs, x = integrate_logz)
        mdf = not_normalized_mdf(logZ = logZ,pmetals=pmetals)
        normed_mdf = mdf/intgrl
        return log_Z_over_Z_sun, normed_mdf

def get_mean_logZ_of_mdf(pmetals = 2):
    '''
    computes the mean logZ of the mdf implemented in fsps if 
    zcontinuous = 2, in solar units. 
    
    parameters
    -----------
    pmetals: float. 
    '''        
    logZ = np.linspace(-10, 5, 10000)
    log_Z_over_Z_sun, dn_dlogZ = normalized_mdf(logZ = logZ,pmetals=pmetals)
    mean_logZ = np.trapz(log_Z_over_Z_sun*dn_dlogZ, x = log_Z_over_Z_sun)
    return mean_logZ

