F90 = gfortran
F90FLAGS = -O3 -cpp -march=native

SRC_DIR = fsps/src

COMMON = sps_vars.o nrtype.o sps_utils.o nr.o nrutil.o compsp.o	 \
ssp_gen.o getmags.o ran.o locate.o qromb.o polint.o trapzd.o tridag.o \
sps_setup.o pz_convol.o get_tuniv.o imf.o imf_weight.o add_dust.o getspec.o \
add_bs.o mod_hb.o add_remnants.o getindx.o velbroad.o mod_gb.o \
write_isochrone.o sfhstat.o linterp.o

MODULE = pyfsps

%.o : $(SRC_DIR)/%.f90
	$(F90) $(F90FLAGS) -o $@ -c $<

%.o : $(SRC_DIR)/nr/%.f90
	$(F90) $(F90FLAGS) -o $@ -c $<

$(MODULE): $(MODULE).f95 $(COMMON)
	rm -f $(MODULE).so
	f2py -c -m $(MODULE) --fcompiler=gnu95 $(MODULE).f95 $(COMMON)

clean:
	rm -rf *.o *.mod *.so

# --f90flags="-fopenmp -lm " --opt="-O3 --ffast-math" -lgomp --fcompiler=gnu95
