FC = g77
TARGET=gfit
TSRC = $(GGGPATH)/src/$(TARGET)
COMN = $(GGGPATH)/src/comn
CFLAGS = -c -O -fbounds-check -Wunused -Wuninitialized -W -Wall
LDFLAGS =

$(TARGET): \
  $(TSRC)/$(TARGET).o \
  $(TSRC)/fetch.o \
  $(TSRC)/solar_gas_shift.o \
  $(TSRC)/write_col.o \
  $(TSRC)/compute_residual.o \
  $(TSRC)/compute_vertical_paths.o \
  $(TSRC)/spectrum_loop.o \
  $(TSRC)/do_retrieval.o \
  $(TSRC)/do_retrieval2.o \
  $(TSRC)/do_retrieval3.o \
  $(TSRC)/subtract_cf.o \
  $(TSRC)/fm.o \
  $(TSRC)/fm2.o \
  $(TSRC)/fm3.o \
  $(TSRC)/fringes.o \
  $(TSRC)/jetspe.o \
  $(TSRC)/write_spt.o \
  $(TSRC)/read_mav.o \
  $(COMN)/abscoi.o \
  $(COMN)/gasdev.o \
  $(COMN)/gct_humlik.o \
  $(COMN)/getendian.o \
  $(COMN)/gindfile.o \
  $(COMN)/file_size_in_bytes.o \
  $(COMN)/hitran_to_atmos_gas_numbering.o \
  $(COMN)/fbc.o \
  $(COMN)/fnbc.o \
  $(COMN)/hfti.o \
  $(COMN)/julian.o \
  $(COMN)/lloc.o \
  $(COMN)/lnbc.o \
  $(COMN)/lowercase.o \
  $(COMN)/modfft.o \
  $(COMN)/newdec.o \
  $(COMN)/posnall.o \
  $(COMN)/profzl.o \
  $(COMN)/ran1.o \
  $(COMN)/rbyte.o \
  $(COMN)/read_runlog.o \
  $(COMN)/skiprec.o \
  $(COMN)/solar_pseudo_trans_spec.o \
  $(COMN)/subsolar.o \
  $(COMN)/substr.o \
  $(COMN)/tlpath.o \
  $(COMN)/vsubs.o \
  $(COMN)/read_isotop.o \
  $(COMN)/wlimit.o \
  $(COMN)/zenaz.o

#  Link the object files to create the executable
	$(FC) $(LDFLAGS) -o $(GGGPATH)/bin/$@ $^

#  Compile the .f files that are newer than the corresponding .o files
%.o : %.f
	$(FC)  $(CFLAGS) -o $@ $<
