BINDIR = $(GGGPATH)/bin
CSRC = $(GGGPATH)/src/comn
TSRC = $(GGGPATH)/src/create_sunrun/izana

#sun solaris f77, compiler flags, linker flags
#FC = f77
#FFLAGS = -c -Nl50 -O
#LDFLAGS =

#linux gnu g77, compiler flags, linker flags
FC = g77
FFLAGS = -c -O -fbounds-check -Wunused -Wuninitialized -W -Wall
LDFLAGS =

#linux absoft f95, compiler flags, linker flags
#FC = f95
#FFLAGS = -c -ea -en -g -g77
#LDFLAGS = -g -g77 -lU77

SRCS = \
	$(TSRC)/create_sunrun_from_izana.f

CSRCS = \
	$(CSRC)/gindfile.f \
	$(CSRC)/lloc.f \
	$(CSRC)/lnbc.f \
	$(CSRC)/fnbc.f \
	$(CSRC)/fbc.f \
	$(CSRC)/rbyte.f \
	$(CSRC)/write_sunrun.f \
	$(CSRC)/get_ggg_environment.f

OBJS = $(SRCS:.f=.o)
COBJS = $(CSRCS:.f=.o)

.SUFFIXES : .f
.f.o :
	$(FC) $(FFLAGS) $<

create_sunrun_from_izana : $(OBJS) comn_objs
	$(FC) $(LDFLAGS) -o $(BINDIR)/$@ $(OBJS) $(COBJS)

comn_objs :
	cd $(CSRC); make create_sunrun

#  Add dependencies to object files
$(TSRC)/create_sunrun_from_izana.o:  $(GGGPATH)/src/ggg_int_params.f $(TSRC)/../params.f
$(CSRC)/gindfile.o:          $(GGGPATH)/src/ggg_int_params.f
$(CSRC)/write_sunrun.o:      $(GGGPATH)/src/ggg_int_params.f

clean :
	/bin/rm -f $(OBJS) $(COBJS)
