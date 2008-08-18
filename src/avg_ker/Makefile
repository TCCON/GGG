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
	avg_ker.f \

CSRCS = \
	$(CSRC)/fbc.f \
	$(CSRC)/fnbc.f \
	$(CSRC)/lnbc.f \
	$(CSRC)/vsubs.f \
	$(CSRC)/hfti.f

OBJS = $(SRCS:.f=.o)
COBJS = $(CSRCS:.f=.o)

.SUFFIXES : .f
.f.o :
	$(FC) $(FFLAGS) $<
 
avg_ker: $(OBJS) comn_objs
	$(FC) $(LDFLAGS) -o $(BINDIR)/$@ $(OBJS) $(COBJS)
 
comn_objs:
	cd $(CSRC); make avg_ker
 
clean:
	/bin/rm -f $(OBJS) $(COBJS)
