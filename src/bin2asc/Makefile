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
#FFLAGS = -c -en -g
#LDFLAGS = -g -lU77

SRCS = \
	bin2asc.f

CSRCS = \
	$(CSRC)/gindfile.f \
	$(CSRC)/lnbc.f \
	$(CSRC)/lloc.f \
	$(CSRC)/rbyte.f \
	$(CSRC)/getendian.f \
        $(CSRC)/write_runlog.f\
        $(CSRC)/wlimit.f\
	$(CSRC)/read_runlog.f

OBJ = $(SRCS:.f=.o)
COBJ = $(CSRCS:.f=.o)

.SUFFIXES : .f

.f.o :
	$(FC) $(FFLAGS) $<

bin2asc : $(OBJ) comn_objs
	$(FC) $(LDFLAGS) -o $(BINDIR)/$@ $(OBJ) $(COBJ)

comn_objs :
	cd $(CSRC); make bin2asc

clean :
	/bin/rm -f $(OBJ) $(COBJ)
