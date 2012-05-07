import numpy as np
import scipy.optimize as op
import scipy.interpolate as interp
import fsps
import sys
import h5py
import emcee
import matplotlib.pyplot as pl


fsps.ssps()


class ModelGalaxy(object):
    def __init__(self, redshift, logz, **kwargs):
        self.redshift = redshift
        mets = fsps.logz
        mi1, mi2 = sorted(np.argsort(np.abs(mets - logz))[:2])
        self.met1, self.met2 = mets[mi1], mets[mi2]
        self.zmet = logz

        fsps.compute(mi1, **kwargs)
        self._mags_1 = fsps.get_mags(redshift=redshift)
        s1 = fsps.get_stats()

        fsps.compute(mi2, **kwargs)
        self._mags_2 = fsps.get_mags(redshift=redshift)
        s2 = fsps.get_stats()

        self.age = self._interp_met(s1["age"], s2["age"])
        self.mass = self._interp_met(s1["mass"], s2["mass"])
        self.lbol = self._interp_met(s1["lbol"], s2["lbol"])
        self.sfr = self._interp_met(s1["sfr"], s2["sfr"])
        self.mdust = self._interp_met(s1["mdust"], s2["mdust"])

    def _interp_met(self, y1, y2):
        m = (y2 - y1) / (self.met2 - self.met1)
        b = y1 - m * self.met1
        return m * self.zmet + b

    def mag(self, band, age):
        m1 = self._mags_1[band]
        m2 = self._mags_2[band]
        m = self._interp_met(m1, m2)
        s = interp.interp1d(self.age, m)
        return s(age)


if __name__ == "__main__":
    age = 9  # Gyr
    gal = ModelGalaxy(0.1, -0.09, sfh=1, const=0.4)

    bands = ["sdss_" + b for b in "ugriz"]
    error = 0.1
    data = np.array([gal.mag(b, age) for b in bands])
    data += error * np.random.randn(len(data))

    opfn = "opt.dat"

    def chi2(p, opt):
        if p[1] < 0 or not -1.98 < p[0] < 0.2:
            return 1e10
        try:
            m = ModelGalaxy(0.1, p[0], sfh=1, const=0.4)
        except Exception:
            return 1e10
        try:
            diff = np.array([m.mag(b, p[1]) for b in bands]) - data
        except ValueError:
            return 1e10
        c2 = np.sum(diff / error) ** 2
        if opt:
            f = open(opfn, "a")
            f.write("{0} {1} {2}\n".format(p[0], p[1], c2))
            f.close()
        print c2
        return c2

    def lnlike(p):
        ll = - 0.5 * chi2(p, False)
        return ll

    truth = np.array([gal.zmet, age])
    if "--mcmc" in sys.argv:
        ndim, nwalkers = len(truth), 10
        p0 = truth[None, :] \
                + 0.01 * np.random.randn(nwalkers * ndim)\
                .reshape((nwalkers, ndim))
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlike)

        fn = "mcmc.h5"
        f = h5py.File(fn, "w")
        f.create_dataset("truth", data=truth)
        f.close()
        for i, (pos, lnprob, state) in \
                                enumerate(sampler.sample(p0, iterations=10)):
            f = h5py.File(fn, "a")
            g = f.create_group(str(i))
            g.create_dataset("pos", data=pos)
            g.create_dataset("lnprob", data=lnprob)
            f.close()
    else:
        f = open(opfn, "w")
        f.write("# {0}\n".format(truth))
        f.close()
        p = op.fmin_bfgs(chi2, [-1, 8], args=(True,))
        f = open(opfn, "a")
        f.write("# optimal = {0}\n".format(p))
        f.close()
