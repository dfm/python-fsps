# We want to run a few basic tests.  These should check that 
#   x1. The code is importable
#   x2. The SP object instantiates and the libraries attrbute is correct
#   x3. The returned spectrum has the right shape
#   x4. The attributes have the right shape
#   x5. The spectrum is not all zeros
#   x6. Nebular emission works for young populations when add_neb_emission=True
# Then we have a few things to look at for different compiler options
#  x1. The SSP color evolution for afe=0 makes sense between C3K_LR, C3K_LR+alpha, and MILES
#  x2. The C3K_LR+alpha spectra are different from the C3K_LR spectra for afe>0
#  x3. The C3K_LR spectra are the same as the C3K_LR+alpha spectra for afe=0

import sys, os, argparse
import numpy as np
from astropy.table import Table
from sedpy.observate import load_filters, getSED


filternames = ['bessell_U', 'bessell_B', 'bessell_V', 'bessell_R', 'twomass_J', 'twomass_Ks']
filternames += ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0']

FILTERS = load_filters(filternames)
#FILTERS = load_filters(['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0'])

def get_githash(dir):
    """
    Get the current git hash for a repo
    """
    import subprocess
    cmd = ["git", "rev-parse", "HEAD"]
    hash = subprocess.check_output(cmd, cwd=dir).decode("utf-8").strip()
    return hash


def test_import_and_instantiate(zcontinuous=1):
    import fsps
    from fsps import StellarPopulation
    sp = StellarPopulation(zcontinuous=zcontinuous)
    print(sp.libraries)
    assert "afeindx" in sp.params.csp_params
    assert "afe" in sp.params.csp_params
    print(f"Number of afe gridpoints: {sp.n_afe}")
    print(f'afe index: {sp.params["afeindx"]}')
    print(f'afe value: {sp.params["afe"]}')
    sp.params["use_lw_tpagb"] = True

    default_params = dict([(k, sp.params[k]) for k in sp.params.all_params])
    return sp, default_params


def _reset_default_params(pop, params, zcontinuous=1):
    pop._zcontinuous = zcontinuous
    for k in pop.params.all_params:
        pop.params[k] = params[k]


def test_spectrum_shape(sp):
    wave, spec = sp.get_spectrum(tage=1.0, peraa=False)
    print(len(wave))
    assert spec.shape == (len(wave),)
    assert not np.all(spec == 0)
    wave, spec = sp.get_spectrum(peraa=True)
    assert spec.shape == (len(sp.ssp_ages), len(wave),)
    assert spec.shape[0]  > 1
    return None


def test_attributes(sp, sfh=3):

    sp.params["sfh"] = sfh
    if sfh == 3:
        sp.set_tabular_sfh(np.linspace(0, 10, 10),
                            np.random.uniform(0, 1, 10))
    wave, spec = sp.get_spectrum(tage=1.0, peraa=True)
    young, old = sp._csp_young_old
    assert young.shape == old.shape
    assert young.shape == (len(wave),)

    w = sp.wavelengths
    res = sp.resolutions
    assert w.shape == res.shape
    assert w.shape == (len(wave),)

    zleg = sp.zlegend
    n_afe = sp.n_afe
    ages = sp.ssp_ages
    weights = sp._ssp_weights
    assert weights.shape == (len(ages), len(zleg))
    assert np.any(weights > 0)
    assert n_afe >= 1

    age = sp.log_age
    mstar = sp.stellar_mass
    lbol = sp.log_lbol
    sfr = sp.sfr
    mdust = sp.dust_mass
    mform = sp.formed_mass
    assert np.isscalar(age)
    assert np.isscalar(mstar)
    assert np.isscalar(lbol)
    assert np.isscalar(sfr)
    assert np.isscalar(mdust)
    assert np.isscalar(mform)
    assert age > 0
    assert mstar > 0
    assert mform > 0
    assert lbol > -33
    assert mstar/mform < 1

    ewave = sp.emline_wavelengths
    elum = sp.emline_luminosity
    assert ewave.shape == elum.shape
    assert ewave.size > 0

    return None


def test_sfhs(sp):

    for sfh in [0, 1, 2, 3, 4, 5]:
        try:
            stat = test_attributes(sp, sfh=sfh)
        except Exception as e:
            print(f"Error testing sfh={sfh}: {e}")
            raise
    return None


def test_neb(sp, tage=0.003, gas_logz=-0.5, gas_logu=-2.0):
    current = sp.params["add_neb_emission"]
    sp.params["add_neb_emission"] = True
    sp.params["logzsol"] = gas_logz
    sp.params["gas_logu"] = gas_logu
    sp.params["gas_logz"] = gas_logz
    sp.params["sfh"] = 0
    wave, spec = sp.get_spectrum(tage=tage, peraa=True)
    eline_waves = sp.emline_wavelengths
    eline_lums = sp.emline_luminosity
    assert eline_waves.shape == eline_lums.shape
    assert np.any(eline_lums > 0)
    sp.params["add_neb_emission"] = current
    return None


def color_evol(sp, afe=0, plot=False, filterlist=FILTERS):
    fn = (1, 2) # B-V
    sp.params["sfh"] = 0
    sp.params["afe"] = afe
    ssp_ages = sp.ssp_ages.copy() + 0.02
    mags = np.zeros((len(ssp_ages), len(filterlist)))
    for it, log_age in enumerate(ssp_ages):
        age = 10**(log_age-9)
        wave, spec = sp.get_spectrum(tage=age, peraa=True)
        mags[it] = getSED(wave, spec, filterlist=filterlist, linear_flux=False)
        #print(f"Age: {age:.2f} Gyr, AFE: {afe}, mags: {mags}")

    # make an astropy table of the mangitudes and ages
    colnames = [f"{filt.name}" for filt in filterlist]
    colnames.insert(0, "log_age")
    t = Table(data=np.column_stack((ssp_ages, mags)), names=colnames)

    if plot:
        import matplotlib.pyplot as pl
        fig, ax = pl.subplots()
        ax.plot(ssp_ages, mags[:, fn[0]] - mags[:, fn[1]], "-o")
        ax.set_xlabel("SSP Age (log yr)")
        ax.set_ylabel(f"Color ({filterlist[fn[0]].name} - {filterlist[fn[1]].name})")
        ax.set_ylim(-1.4, 1.8)
        ax.set_xlim(4.8, 10.4)
        return t, fig
    else:
        return t, None


def get_parser():
    parser = argparse.ArgumentParser(description="Run pyFSPS AFE tests")
    parser.add_argument("--outdir", type=str, default="./lib_tests",
                        help="Output directory for test results")
    parser.add_argument("--show_hashes", action="store_true",
                        help="Only show git hashes for FSPS and pyFSPS")
    parser.add_argument("--logzsol", type=float, default=0.0,
                        help="Metallicity for the tests")
    parser.add_argument("--plot", type=int, default=1,
                        help="Plot the color evolution (1) or not (0)")
    parser.add_argument("--write", type=int, default=1,
                        help="Write the results to disk (1) or not (0)")
    return parser


if __name__ == "__main__":

    parser = get_parser()
    args = parser.parse_args()
    outdir = args.outdir

    pyfsps_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    fsps_hash = get_githash(os.environ["SPS_HOME"])
    pyfsps_hash = get_githash(os.path.dirname(__file__))
    pyfsps_lib_hash = get_githash(os.path.join(pyfsps_dir, "src", "fsps", "libfsps"))
    print(f"FSPS git hash: {fsps_hash}")
    print(f"pyFSPS git hash: {pyfsps_hash}")
    print(f"pyFSPS lib git hash: {pyfsps_lib_hash}")

    hashtag = f"sps_home-{fsps_hash[:6]}_pyfsps-{pyfsps_hash[:6]}_pyfslib-{pyfsps_lib_hash[:6]}"
    if args.show_hashes:
        sys.exit()

    print(f"writing results to {outdir}/{hashtag}")
    if args.write | args.plot:
        os.makedirs(f"{outdir}/{hashtag}", exist_ok=True)

    sp, defaults = test_import_and_instantiate(zcontinuous=1)
    defaults["logzsol"] = args.logzsol

    libs, nafe = [s.decode("utf-8") for s in sp.libraries], sp.n_afe
    tag = f"{libs[0]}+{libs[1]}+afe{int(nafe > 1)}"

    test_spectrum_shape(sp)
    _reset_default_params(sp, defaults, zcontinuous=1)
    test_attributes(sp, sfh=0)
    _reset_default_params(sp, defaults, zcontinuous=1)
    test_sfhs(sp)
    _reset_default_params(sp, defaults, zcontinuous=1)
    test_neb(sp)
    _reset_default_params(sp, defaults, zcontinuous=1)

    # get all the color data
    table, fig = color_evol(sp, afe=0, plot=True)

    if args.write:
        table.write(f"{outdir}/{hashtag}/color_evol_{tag}.csv", format="csv", overwrite=True)
    if args.plot:
        fig.axes[0].grid(True)
        fig.axes[0].set_title(f"{tag}")
        fig.savefig(f"{outdir}/{hashtag}/color_evol_{tag}.png", dpi=300)
    if nafe > 1:
        # just run one alpha_enhancement
        t3, fig3 = color_evol(sp, afe=0.3, plot=args.plot)
        if args.write:
            t3.write(f"{outdir}/{hashtag}/color_evol_{tag}_afe0.3.csv", format="csv", overwrite=True)
        if args.plot:
            fig3.axes[0].set_title(f"{tag}, afe=0.3")
            fig3.axes[0].grid(True)
            fig3.axes[0].plot(table["log_age"], table["bessell_B"] - table["bessell_V"], "-o", alpha=0.5, label="[a/Fe]=0.0")
            fig3.axes[0].legend()
            fig3.savefig(f"{outdir}/{hashtag}/color_evol_{tag}_afe0.3.png", dpi=300)
    _reset_default_params(sp, defaults, zcontinuous=1)

    # -- make spectra for different SFHs --
    sp.params["logzsol"] = args.logzsol
    sp.params["sfh"] = 1
    sp.params["const"] = 1.0
    wave, spec = sp.get_spectrum(tage=1, peraa=True)
    tbl = Table(data=np.column_stack((wave, spec)), names=["wave", "spec"])
    if args.write:
        tbl.write(f"{outdir}/{hashtag}/spectrum_{tag}_const.csv", format="csv", overwrite=True)
    _reset_default_params(sp, defaults, zcontinuous=1)

    sp.params["logzsol"] = args.logzsol
    sp.params["sfh"] = 1
    sp.params["const"] = 1.0
    sp.params["add_neb_emission"] = True
    wave, spec = sp.get_spectrum(tage=1, peraa=True)
    tbl = Table(data=np.column_stack((wave, spec)), names=["wave", "spec"])
    if args.write:
        tbl.write(f"{outdir}/{hashtag}/spectrum_{tag}_constneb.csv", format="csv", overwrite=True)
    _reset_default_params(sp, defaults, zcontinuous=1)

    sp.params["logzsol"] = args.logzsol
    sp.params["sfh"] = 1
    sp.params["const"] = 0.0
    sp.params["tau"] = 1.0
    wave, spec = sp.get_spectrum(tage=1, peraa=True)
    tbl = Table(data=np.column_stack((wave, spec)), names=["wave", "spec"])
    if args.write:
        tbl.write(f"{outdir}/{hashtag}/spectrum_{tag}_tau1.csv", format="csv", overwrite=True)
    _reset_default_params(sp, defaults, zcontinuous=1)

