*This is a work in progress and it's probably broken. Use at your own risk.*

See [FSPS](http://www.ucolick.org/~cconroy/FSPS.html) by Charlie Conroy.

To get started:

1. [Install FSPS](http://www.ucolick.org/~cconroy/FSPS.html).
2. Make sure that your `$SPS_HOME` environment variable is set.
3. Run `python setup.py build_fsps` in the `python-fsps` root directory.

Basic usage:

```python
import fsps
fsps.ssps() # Compute all the SSPs
fsps.compute(5) # Compute the CSP at metallicity=5
lam = fsps.get_lambda()
spec = fsps.get_spec()

import numpy as np
import matplotlib.pyplot as pl
pl.plot(np.log10(lam), spec)
```

