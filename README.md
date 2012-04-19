Work in progress...

See [FSPS](https://www.cfa.harvard.edu/~cconroy/FSPS.html) by Charlie Conroy.

To get started:

```bash
svn checkout http://fsps.googlecode.com/svn/trunk/ fsps/fsps
export SPS_HOME=`pwd`/fsps/fsps
make
python -c 'import fsps'
```

Basic usage:

```python
import fsps
fsps.compute(zmet=10) # Other options coming soon.
g = fsps.get_mag("sdss_g")
r = fsps.get_mag("sdss_r")
lam, flux = fsps.get_spec()
```

