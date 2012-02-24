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
  $(TSRC)/fm3.f \
  $(TSRC)/fringes.f \
  $(TSRC)/jetspe.f \
  $(TSRC)/write_spt.f \
  $(COMN)/read_mavfile_body.f \
  $(COMN)/get_ggg_environment.f \
  $(COMN)/abscoj.f \
  $(COMN)/gasdev.f \
  $(COMN)/gct2_humlik.f \
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
  $(COMN)/solar_pts.f \
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

#  Add dependencies to object files
$(TSRC)/do_retrieval.o:  $(GGGPATH)/src/ggg_const_params.f \
                         $(TSRC)/const_params.f $(TSRC)/int_params.f
$(TSRC)/do_retrieval3.o: $(GGGPATH)/src/ggg_const_params.f \
                         $(TSRC)/const_params.f $(TSRC)/int_params.f
$(TSRC)/fm.o:            $(GGGPATH)/src/ggg_const_params.f \
                         $(GGGPATH)/src/ggg_int_params.f $(TSRC)/const_params.f
$(TSRC)/fm3.o:           $(GGGPATH)/src/ggg_const_params.f $(TSRC)/const_params.f
$(TSRC)/fringes.o:       $(GGGPATH)/src/ggg_const_params.f
$(TSRC)/gfit.o:          $(GGGPATH)/src/ggg_int_params.f $(TSRC)/int_params.f
$(TSRC)/jetspe.o:        $(GGGPATH)/src/ggg_const_params.f
$(TSRC)/spectrum_loop.o: $(GGGPATH)/src/ggg_const_params.f \
                         $(GGGPATH)/src/ggg_int_params.f   \
                         $(TSRC)/const_params.f $(TSRC)/int_params.f
$(TSRC)/subtract_cf.o:   $(GGGPATH)/src/ggg_const_params.f
$(COMN)/abscoj.o:        $(GGGPATH)/src/ggg_const_params.f $(GGGPATH)/src/ggg_int_params.f
$(COMN)/gindfile.o:      $(GGGPATH)/src/ggg_int_params.f
$(COMN)/profzl.o:        $(GGGPATH)/src/ggg_const_params.f
$(COMN)/tlpath.o:        $(GGGPATH)/src/ggg_const_params.f $(GGGPATH)/src/ggg_int_params.f

clean:
	rm -f $(OBJECTS) $(GGGPATH)/bin/$(TARGET)

