*This is a work in progress and it's probably broken. Use at your own risk.*

See [FSPS](https://www.cfa.harvard.edu/~cconroy/FSPS.html) by Charlie Conroy.

To get started:

1. [Install FSPS](https://www.cfa.harvard.edu/~cconroy/FSPS.html).
2. Make sure that your `$SPS_HOME` environment variable is set.
3. Run `python setup.py build_fsps` in the `python-fsps` root directory.

Basic usage:

```python
import fsps
fsps.ssps() # Compute all the SSPs
fsps.compute(5) # Compute the CSP at metallicity=5
lambda = fsps.get_lambda()
spec = fsps.get_spec()

import numpy as np
import matplotlib.pyplot as pl
pl.plot(np.log10(lambda), spec)
```

