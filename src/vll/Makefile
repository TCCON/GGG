BINDIR = $(GGGPATH)/bin
CSRC = $(GGGPATH)/src/comn

#sun solaris f77, compiler flags, linker flags
#FC = f77
#FFLAGS = -c -Nl50 -O
#LDFLAGS =

#linux gnu g77, compiler flags, linker flags
FC = g77
FFLAGS = -c -C -pedantic -Wimplicit -Wunused -Wuninitialized -Wsurprising -W
LDFLAGS =

#linux absoft f95, compiler flags, linker flags
#FC = f95
#FFLAGS = -c -en -g
#LDFLAGS = -g -lU77

SRCS = \
	vll.f

CSRCS = \
	$(CSRC)/lnbc.f \
	$(CSRC)/lloc.f \
	$(CSRC)/hitran_to_atmos_gas_numbering.f \
	$(CSRC)/file_size_in_bytes.f \
	$(CSRC)/posnall.f

OBJ = $(SRCS:.f=.o)
COBJ = $(CSRCS:.f=.o)

.SUFFIXES : .f

.f.o :
	$(FC) $(FFLAGS) $<

vll : $(OBJ) comn_objs
	$(FC) $(LDFLAGS) -o $(BINDIR)/$@ $(OBJ) $(COBJ)

comn_objs :
	cd $(CSRC); make vll

clean :
	/bin/rm -f $(OBJ) $(COBJ)
