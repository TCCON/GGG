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
	ak.f \
	fetch.f \
	fm.f \
	fringes.f \
	gfit.f \
	jetspe.f \
	write_spt.f \
	read_mav.f

CSRCS = \
	$(CSRC)/abscoh.f \
	$(CSRC)/gasdev.f \
	$(CSRC)/gct_humlik.f \
	$(CSRC)/getendian.f \
	$(CSRC)/gindfile.f \
	$(CSRC)/file_size_in_bytes.f \
	$(CSRC)/hitran_to_atmos_gas_numbering.f \
	$(CSRC)/fbc.f \
	$(CSRC)/fnbc.f \
	$(CSRC)/hfti.f \
	$(CSRC)/julian.f \
	$(CSRC)/lloc.f \
	$(CSRC)/lnbc.f \
	$(CSRC)/lowercase.f \
	$(CSRC)/modfft.f \
	$(CSRC)/newdec.f \
	$(CSRC)/posnall.f \
	$(CSRC)/profzl.f \
	$(CSRC)/ran1.f \
	$(CSRC)/rbyte.f \
	$(CSRC)/read_runlog.f \
	$(CSRC)/skiprec.f \
	$(CSRC)/solar_spectrum.f \
	$(CSRC)/subsolar.f \
	$(CSRC)/substr.f \
	$(CSRC)/vsubs.f \
	$(CSRC)/read_isotop.f \
	$(CSRC)/zenaz.f

OBJ = $(SRCS:.f=.o)
COBJ = $(CSRCS:.f=.o)

.SUFFIXES : .f

.f.o :
	$(FC) $(FFLAGS) $<

gfit : $(OBJ) comn_objs
	$(FC) $(LDFLAGS) -o $(BINDIR)/$@ $(OBJ) $(COBJ)

comn_objs :
	cd $(CSRC); make gfit

clean :
	/bin/rm -f $(OBJ) $(COBJ)
