F77FLAGS=-Nl65 -O -c -u
G77FLAGS=-c -C -pedantic -Wimplicit -Wunused -Wuninitialized -Wsurprising -W
ALL=abscog.o abscoh.o byterev.o clistindex.o endian.o fbc.o findfile.o fnbc.o \
	gasdev.o gct_humlik.o getendian.o getopusparval.o gindfile.o height.o \
	hfti.o hitran2atmos.o julian.o lloc.o lnbc.o lowercase.o modfft.o\
	molnum.o newdec.o plint.o posnall.o profzl.o ran1.o rbyte.o \
	read_isotop.o read_mkiv_header.o read_opus_header.o readnew.o \
	read_runlog.o read_sunrun.o skiprec.o slpath.o file_size_in_bytes.o \
	solar_spectrum.o subsolar.o substr.o sunspect.o tlpath.o twodint.o \
	vsubs.o write_runlog.o write_sunrun.o xloc.o zenaz.o
ADJUST_SOLAR_LINELIST=lnbc.o hfti.o posnall.o vsubs.o
AKER=fbc.o fnbc.o lnbc.o hfti.o vsubs.o
AVG_KER=fbc.o fnbc.o lnbc.o hfti.o vsubs.o
ALTCORR=lnbc.o
ALTADJ=fbc.o fnbc.o lnbc.o molnum.o substr.o
AMAL=lnbc.o
ASZA=subsolar.o zenaz.o julian.o
BIN2ASC=rbyte.o getendian.o gindfile.o lnbc.o read_runlog.o
BIRDFIX=byterev.o endian.o findfile.o lnbc.o readnew.o vsubs.o
CHECKINT=findfile.o lnbc.o
CONVERT=lnbc.o lowercase.o
DAYAV=fbc.o fnbc.o lnbc.o substr.o
DIURNRET=fbc.o fnbc.o julian.o lnbc.o lowercase.o skiprec.o \
	subsolar.o substr.o twodint.o zenaz.o
ENG=byterev.o fbc.o findfile.o fnbc.o lnbc.o substr.o
EXTRACT_TPH=read_runlog.o log1pxox.o
GFIT=abscoh.o fbc.o fnbc.o gct_humlik.o getendian.o gindfile.o \
	julian.o lloc.o lnbc.o lowercase.o modfft.o newdec.o posnall.o \
	skiprec.o solar_spectrum.o read_runlog.o subsolar.o substr.o vsubs.o \
	gasdev.o ran1.o rbyte.o zenaz.o profzl.o hfti.o read_isotop.o \
	hitran_to_atmos_gas_numbering.o file_size_in_bytes.o
GFIT.OLD=endian.o fbc.o findfile.o fnbc.o lnbc.o modfft.o newdec.o posnall.o \
	profzl.o skiprec.o  substr.o vsubs.o
GGGAVG=fbc.o fnbc.o lnbc.o skiprec.o substr.o read_runlog.o
GGGAVG.OLD=fbc.o fnbc.o lnbc.o skiprec.o substr.o
GPSDADS=plint.o skiprec.o
GSETUP=clistindex.o fbc.o fnbc.o height.o lnbc.o lowercase.o tlpath.o \
	tlpath.o read_runlog.o substr.o log1pxox.o read_isotop.o
GSETUP.OLD=fbc.o fnbc.o height.o lnbc.o molnum.o slpath.o substr.o
HALOE=julian.o lnbc.o skiprec.o
IPP=byterev.o endian.o findfile.o lnbc.o julian.o modfft.o readnew.o \
	subsolar.o write_runlog.o vsubs.o xloc.o zenaz.o
MAPTELEM=byterev.o endian.o findfile.o lnbc.o
PTDAO=skiprec.o xloc.o
RATIO=findfile.o lnbc.o newdec.o profzl.o vsubs.o zenaz.o subsolar.o
RH94=lloc.o lnbc.o posnall.o
CREATE_RUNLOG=getopusparval.o byterev.o getendian.o gindfile.o lnbc.o \
	subsolar.o zenaz.o read_sunrun.o write_runlog.o julian.o \
	read_mkiv_header.o rbyte.o fnbc.o fbc.o
CREATE_SUNRUN=read_mkiv_header.o read_opus_header.o getopusparval.o \
	byterev.o getendian.o \
	gindfile.o lnbc.o write_sunrun.o fnbc.o fbc.o rbyte.o
SOLAR=abscog.o fbc.o fnbc.o gct_humlik.o getendian.o gindfile.o lloc.o \
	lnbc.o newdec.o posnall.o profzl.o rbyte.o substr.o vsubs.o
SOLAR.OLD=findfile.o lnbc.o newdec.o profzl.o vsubs.o
SPAVG=getendian.o write_runlog.o gindfile.o lnbc.o
SS2AMES1001=fbc.o fnbc.o lnbc.o skiprec.o substr.o
SUMMARIZE_LINELIST=read_isotop.o hitran_to_atmos_gas_numbering.o
VLL=hitran2atmos.o posnall.o lnbc.o lloc.o
WINDOWS=posnall.o
ZENANG=fbc.o fnbc.o lnbc.o molnum.o substr.o
ZENCORR=lnbc.o read_runlog.o write_runlog.o

.SUFFIXES : .o .f
.f.o :
	g77 $(G77FLAGS) $<

adjust_solar_linlist: $(ADJUST_SOLAR_LINELIST)
aker: $(AKER)
avg_ker: $(AVG_KER)
all: $(ALL)
altcorr: $(ZENCORR)
altadj : $(ALTADJ)
amal : $(AMAL)
asza : $(ASZA)
bin2asc : $(BIN2ASC)
birdfix : $(BIRDFIX)
checkint : $(CHECKINT)
convert : $(CONVERT)
dayav : $(DAYAV)
diurnret : $(DIURNRET)
diurnret.old : $(DIURNRET)
eng : $(ENG)
extract_pth : $(EXTRACT_PTH)
gfit: $(GFIT)
gfit.dc8 : $(GFIT.OLD)
gfit.old : $(GFIT.OLD)
gggavg : $(GGGAVG)
gggavg.old : $(GGGAVG.OLD)
gpsdads : $(GPSDADS)
gsetup : $(GSETUP)
gsetup.dc8 : $(GSETUP.OLD)
gsetup.old : $(GSETUP.OLD)
haloe : $(HALOE)
ipp: $(IPP)
maptelem : $(MAPTELEM)
newret : $(NEWRET)
ptdao : $(PTDAO)
ratio : $(RATIO)
rh94 : $(RH94)
create_runlog : $(CREATE_RUNLOG)
create_sunrun: $(CREATE_SUNRUN)
solar : $(SOLAR)
solar.dc8 : $(SOLAR.OLD)
spavg: $(SPAVG)
ss2ames1001 : $(SS2AMES1001)
summarize_linelist : $(SUMMARIZE_LINELIST)
vll: $(VLL)
windows : $(WINDOWS)
zenang : $(ZENANG)
zencorr: $(ZENCORR)

clean :
	rm -f $(ALL)
