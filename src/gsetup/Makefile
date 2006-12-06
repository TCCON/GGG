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
#FFLAGS = -c -en -g
#LDFLAGS = -g -lU77

SRCS = \
	gsetup.f \
	write_mav.f \
	setvmr.f \
	gravity.f \
	readmenu.f \
	readmodFC.f \
	readvmrFC.f

CSRCS = \
	$(CSRC)/clistindex.f \
	$(CSRC)/fbc.f \
	$(CSRC)/fnbc.f \
	$(CSRC)/height.f \
	$(CSRC)/lnbc.f \
	$(CSRC)/lowercase.f \
	$(CSRC)/read_runlog.f \
	$(CSRC)/read_isotop.f \
	$(CSRC)/log1pxox.f \
	$(CSRC)/tlpath.f \
	$(CSRC)/substr.f

OBJ = $(SRCS:.f=.o)
COBJ = $(CSRCS:.f=.o)

.SUFFIXES : .f

.f.o :
	$(FC) $(FFLAGS) $<
    
gsetup : $(OBJ) comn_objs
	$(FC) $(LDFLAGS) -o $(BINDIR)/$@ $(OBJ) $(COBJ)

comn_objs :
	cd $(CSRC); make gsetup
    
clean :
	/bin/rm -f $(OBJ) $(COBJ)
