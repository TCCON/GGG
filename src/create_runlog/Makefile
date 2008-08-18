BINDIR = $(GGGPATH)/bin
CSRC = $(GGGPATH)/src/comn

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
	rdsphead.f \
	create_runlog.f

CSRCS = \
	$(CSRC)/read_opus_header.f \
	$(CSRC)/getopusparval.f \
	$(CSRC)/read_mkiv_header.f \
	$(CSRC)/getendian.f \
	$(CSRC)/gindfile.f \
	$(CSRC)/lloc.f \
	$(CSRC)/julian.f \
	$(CSRC)/lnbc.f \
	$(CSRC)/fnbc.f \
	$(CSRC)/fbc.f \
	$(CSRC)/rbyte.f \
	$(CSRC)/subsolar.f \
	$(CSRC)/read_sunrun.f \
	$(CSRC)/write_runlog.f \
	$(CSRC)/wlimit.f \
	$(CSRC)/zenaz.f

OBJS = $(SRCS:.f=.o)
COBJS = $(CSRCS:.f=.o)

.SUFFIXES : .f
.f.o :
	$(FC) $(FFLAGS) $<

create_runlog : $(OBJS) comn_objs
	$(FC) $(LDFLAGS) -o $(BINDIR)/$@ $(OBJS) $(COBJS)

comn_objs :
	cd $(CSRC); make create_runlog

clean :
	/bin/rm -f $(OBJS) $(COBJS)
