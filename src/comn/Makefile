G77FLAGS=-c -O -fbounds-check -Wunused -Wuninitialized -W -Wall
F90FLAGS=-c 
ACE_IPP=lnbc.o lloc.o getendian.o write_runlog.o file_size_in_bytes.o rbyte.o
ADJUST_SOLAR_LINELIST=lnbc.o hfti.o posnall.o vsubs.o
AVG_KER=fbc.o fnbc.o lnbc.o hfti.o vsubs.o
BIN2ASC=rbyte.o getendian.o gindfile.o lloc.o lnbc.o read_runlog.o write_runlog.o wlimit.o
CREATE_RUNLOG=read_opus_header.o getopusparval.o getendian.o gindfile.o lnbc.o \
	subsolar.o zenaz.o read_sunrun.o write_runlog.o julian.o \
	read_mkiv_header.o rbyte.o fnbc.o fbc.o lloc.o wlimit.o
CREATE_SUNRUN=read_mkiv_header.o read_opus_header.o getopusparval.o \
	getendian.o gindfile.o lloc.o lnbc.o write_sunrun.o fnbc.o fbc.o rbyte.o
DAYAV=fbc.o fnbc.o lnbc.o substr.o
DIURNRET=fbc.o fnbc.o julian.o lnbc.o lowercase.o skiprec.o \
	subsolar.o substr.o twodint.o zenaz.o
EXTRACT_PTH=read_runlog.o lnbc.o log1pxox.o julian.o caldat.o
GFIT=abscoi.o fbc.o fnbc.o lloc.o gct_humlik.o getendian.o gindfile.o \
	julian.o lloc.o lnbc.o lowercase.o modfft.o newdec.o posnall.o \
	skiprec.o solar_pseudo_trans_spec.o read_runlog.o \
	subsolar.o substr.o vsubs.o tlpath.o wlimit.o \
	gasdev.o ran1.o rbyte.o zenaz.o profzl.o hfti.o read_isotop.o \
	hitran_to_atmos_gas_numbering.o file_size_in_bytes.o
GFIT_BJC=abscoh.o fbc.o fnbc.o lloc.o gct_humlik.o getendian.o gindfile.o \
	julian.o lloc.o lnbc.o lowercase.o modfft.o newdec.o posnall.o \
	skiprec.o solar_pseudo_trans_spec.o read_runlog.o \
	subsolar.o substr.o vsubs.o tlpath.o wlimit.o \
	gasdev.o ran1.o rbyte.o zenaz.o profzl.o hfti.o dhfti.o read_isotop.o \
	hitran_to_atmos_gas_numbering.o file_size_in_bytes.o
GGGAVG=fbc.o fnbc.o lnbc.o skiprec.o substr.o read_runlog.o
GSETUP=fbc.o fnbc.o height.o lnbc.o lowercase.o tlpath.o \
	read_sunrun.o read_runlog.o substr.o log1pxox.o read_isotop.o
LIST_MAKER=lnbc.o fnbc.o gindfile.o lloc.o julian.o caldat.o
NONGRAY=rbyte.o getendian.o gindfile.o lloc.o lnbc.o read_runlog.o
SPEC_AVG=rbyte.o getendian.o gindfile.o lloc.o lnbc.o read_runlog.o write_runlog.o
SS2AMES1001=fbc.o fnbc.o lnbc.o skiprec.o substr.o
SUMMARIZE_LINELIST=read_isotop.o hitran_to_atmos_gas_numbering.o
TEXO=getendian.o gindfile.o gasdev.o ran1.o fnbc.o lloc.o lnbc.o profzl.o rbyte.o \
        read_runlog.o write_runlog.o vsubs.o
VLL=hitran_to_atmos_gas_numbering.o file_size_in_bytes.o posnall.o lnbc.o lloc.o
XCO2_CALIBRATE=lnbc.o fnbc.o fbc.o substr.o
ZENANG=fbc.o fnbc.o lnbc.o substr.o
ZENCORR=lnbc.o read_runlog.o write_runlog.o

.SUFFIXES : .o .f
.f.o :
	gfortran $(G77FLAGS) $<
#	f90 $(F90FLAGS) $<

ace_ipp: $(ACE_IPP)
avg_ker: $(AVG_KER)
bin2asc : $(BIN2ASC)
dayav : $(DAYAV)
diurnret : $(DIURNRET)
extract_pth : $(EXTRACT_PTH)
gfit: $(GFIT)
gfit_bjc: $(GFIT_BJC)
gggavg : $(GGGAVG)
gpsdads : $(GPSDADS)
gsetup : $(GSETUP)
list_maker: $(LIST_MAKER)
create_runlog : $(CREATE_RUNLOG)
create_sunrun: $(CREATE_SUNRUN)
spec_avg: $(SPEC_AVG)
ss2ames1001 : $(SS2AMES1001)
summarize_linelist : $(SUMMARIZE_LINELIST)
texo: $(TEXO)
vll: $(VLL)
xco2_calibrate: $(XCO2_CALIBRATE)
zenang : $(ZENANG)
zencorr: $(ZENCORR)

clean :
	rm -f $(ALL)
