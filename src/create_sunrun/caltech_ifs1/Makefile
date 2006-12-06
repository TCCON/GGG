BINDIR = $(GGGPATH)/bin
CSRC = $(GGGPATH)/src/comn

#sun solaris f77, compiler flags, linker flags
#FC = f77
#FFLAGS = -c -Nl50 -O
#LDFLAGS =

#linux gnu g77, compiler flags, linker flags
FC = g77
FFLAGS = -c -O -pedantic -Wimplicit -Wunused -Wuninitialized -Wsurprising -W
LDFLAGS =

#linux absoft f95, compiler flags, linker flags
#FC = f95
#FFLAGS = -c -ea -en -g -g77
#LDFLAGS = -g -g77 -lU77

SRCS = \
	create_sunrun_from_opus_ifs125.f

CSRCS = \
	$(CSRC)/gindfile.f \
	$(CSRC)/lnbc.f \
	$(CSRC)/fnbc.f \
	$(CSRC)/fbc.f \
	$(CSRC)/rbyte.f \
	$(CSRC)/write_sunrun.f \

OBJS = $(SRCS:.f=.o)
COBJS = $(CSRCS:.f=.o)

.SUFFIXES : .f
.f.o :
	$(FC) $(FFLAGS) $<

create_sunrun_from_opus_ifs125 : $(OBJS) comn_objs
	$(FC) $(LDFLAGS) -o $(BINDIR)/$@ $(OBJS) $(COBJS)

comn_objs :
	cd $(CSRC); make create_sunrun

clean :
	/bin/rm -f $(OBJS) $(COBJS)
