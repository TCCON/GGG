G77FLAGS=-c -O -fbounds-check -Wunused -Wuninitialized -W -Wall
AVG_KER=fbc.o fnbc.o lnbc.o hfti.o vsubs.o
BIN2ASC=rbyte.o getendian.o gindfile.o lloc.o lnbc.o read_runlog.o write_runlog.o wlimit.o
CREATE_RUNLOG=read_opus_header.o getopusparval.o getendian.o gindfile.o lnbc.o \
	subsolar.o zenaz.o read_sunrun.o write_runlog.o julian.o \
	read_mkiv_header.o rbyte.o fnbc.o fbc.o lloc.o wlimit.o
DAYAV=fbc.o fnbc.o lnbc.o substr.o
EXTRACT_PTH=read_runlog.o lnbc.o log1pxox.o julian.o caldat.o
GFIT=abscoi.o fbc.o fnbc.o lloc.o gct_humlik.o getendian.o gindfile.o \
	julian.o lloc.o lnbc.o lowercase.o modfft.o newdec.o posnall.o \
	skiprec.o solar_pseudo_trans_spec.o read_runlog.o \
	subsolar.o substr.o vsubs.o tlpath.o wlimit.o \
	gasdev.o ran1.o rbyte.o zenaz.o profzl.o hfti.o read_isotop.o \
	hitran_to_atmos_gas_numbering.o file_size_in_bytes.o
GGGAVG=fbc.o fnbc.o lnbc.o skiprec.o substr.o read_runlog.o
GSETUP=fbc.o fnbc.o height.o lnbc.o lowercase.o tlpath.o \
	read_sunrun.o read_runlog.o substr.o log1pxox.o read_isotop.o
LIST_MAKER=lnbc.o fnbc.o gindfile.o lloc.o julian.o caldat.o
SPEC_AVG=rbyte.o getendian.o gindfile.o lloc.o lnbc.o read_runlog.o write_runlog.o
SS2AMES1001=fbc.o fnbc.o lnbc.o skiprec.o substr.o
SUMMARIZE_LINELIST=read_isotop.o hitran_to_atmos_gas_numbering.o
VLL=hitran_to_atmos_gas_numbering.o file_size_in_bytes.o posnall.o lnbc.o lloc.o

.SUFFIXES : .o .f
.f.o :
	g77 $(G77FLAGS) $<

avg_ker: $(AVG_KER)
bin2asc : $(BIN2ASC)
dayav : $(DAYAV)
extract_pth : $(EXTRACT_PTH)
gfit: $(GFIT)
gggavg : $(GGGAVG)
gsetup : $(GSETUP)
list_maker: $(LIST_MAKER)
create_runlog : $(CREATE_RUNLOG)
spec_avg: $(SPEC_AVG)
ss2ames1001 : $(SS2AMES1001)
summarize_linelist : $(SUMMARIZE_LINELIST)
vll: $(VLL)

clean :
	rm -f $(ALL)
