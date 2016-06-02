F90 = gfortran
F90FLAGS = -O -cpp -fPIC

PROGS = simple lesssimple autosps

COMMON = sps_vars.o sps_utils.o compsp.o csp_gen.o ssp_gen.o getmags.o locate.o funcint.o \
sps_setup.o pz_convol.o get_tuniv.o intsfwght.o imf.o imf_weight.o add_dust.o getspec.o sbf.o \
add_bs.o mod_hb.o add_remnants.o getindx.o smoothspec.o mod_gb.o add_nebular.o \
write_isochrone.o sfhstat.o linterp.o tsum.o add_agb_dust.o linterparr.o \
ztinterp.o vacairconv.o igm_absorb.o get_lumdist.o attn_curve.o \
sfh_weight.o sfhlimit.o sfhinfo.o setup_tabular_sfh.o agn_dust.o

all : $(PROGS)

clean :
	rm -rf *.o *.mod *.MOD *~

autosps : autosps.o $(COMMON)
	$(F90) -o autosps.exe autosps.o $(COMMON)

simple : simple.o $(COMMON)
	$(F90) -o simple.exe simple.o $(COMMON)

lesssimple : lesssimple.o $(COMMON)
	$(F90) -o lesssimple.exe lesssimple.o $(COMMON)

autosps.o : sps_vars.o sps_utils.o

simple.o : sps_vars.o  sps_utils.o

lesssimple.o : sps_vars.o sps_utils.o

sps_utils.o :  sps_vars.o

%.o : %.f90
	$(F90) $(F90FLAGS) -o $@ -c $<

%.o : nr/%.f90
	$(F90) $(F90FLAGS) -o $@ -c $<

