import fsps

sps = fsps.StellarPopulation()
assert sps.libraries[0] == b"pdva", "Incorrect isochrone library"
assert sps.libraries[1] == b"basel", "Incorrect spectral library"
