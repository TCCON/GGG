# include file contains: chosen fortran compiler and link options
include $(GGGPATH)/install/.compiler_ggg

TARGET=gfit
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn

SOURCES := \
  $(TSRC)/$(TARGET).f \
  $(TSRC)/fetch.f \
  $(TSRC)/solar_gas_shift.f \
  $(TSRC)/write_col.f \
  $(TSRC)/compute_residual.f \
  $(TSRC)/compute_vertical_paths.f \
  $(TSRC)/spectrum_loop.f \
  $(TSRC)/do_retrieval.f \
  $(TSRC)/do_retrieval3.f \
  $(TSRC)/subtract_cf.f \
  $(TSRC)/fm.f \
  $(TSRC)/fm2.f \
  $(TSRC)/fm3.f \
  $(TSRC)/fringes.f \
  $(TSRC)/jetspe.f \
  $(TSRC)/write_spt.f \
  $(TSRC)/read_mav.f \
  $(COMN)/get_ggg_environment.f \
  $(COMN)/abscoi.f \
  $(COMN)/gasdev.f \
  $(COMN)/gct_humlik.f \
  $(COMN)/getendian.f \
  $(COMN)/gindfile.f \
  $(COMN)/file_size_in_bytes.f \
  $(COMN)/hitran_to_atmos_gas_numbering.f \
  $(COMN)/fbc.f \
  $(COMN)/fnbc.f \
  $(COMN)/hfti.f \
  $(COMN)/julian.f \
  $(COMN)/lloc.f \
  $(COMN)/lnbc.f \
  $(COMN)/lowercase.f \
  $(COMN)/modfft.f \
  $(COMN)/newdec.f \
  $(COMN)/posnall.f \
  $(COMN)/profzl.f \
  $(COMN)/ran1.f \
  $(COMN)/rbyte.f \
  $(COMN)/read_runlog.f \
  $(COMN)/skiprec.f \
  $(COMN)/solar_pseudo_trans_spec.f \
  $(COMN)/subsolar.f \
  $(COMN)/substr.f \
  $(COMN)/tlpath.f \
  $(COMN)/vsubs.f \
  $(COMN)/read_isotop.f \
  $(COMN)/wlimit.f  

OBJECTS := $(patsubst %.f, %.o, $(SOURCES))

#  Compile the .f files that are newer than the corresponding .o files
.f.o:
	$(FC)  $(CFLAGS) -o $@ $<

#  Link the object files to create the executable
$(TARGET) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^
