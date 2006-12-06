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
	ss2ames1001.f

CSRCS = \
	$(CSRC)/fbc.f \
	$(CSRC)/fnbc.f \
	$(CSRC)/lnbc.f \
	$(CSRC)/skiprec.f \
	$(CSRC)/substr.f

OBJ = $(SRCS:.f=.o)
COBJ = $(CSRCS:.f=.o)

.SUFFIXES : .f

.f.o :
	$(FC) $(FFLAGS) $<
    
ss2ames1001 : $(OBJ) comn_objs
	$(FC) $(LDFLAGS) -o $(BINDIR)/$@ $(OBJ) $(COBJ)

comn_objs :
	cd $(CSRC); make ss2ames1001
    
clean :
	/bin/rm -f $(OBJ) $(COBJ)
