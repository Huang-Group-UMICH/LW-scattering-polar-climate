
FC = ifort
FCOPTS = -w -O2 -r8 -i8
#FCOPTK = -O0 
EXENAME = test111

###########################################
# Source, module and object names and paths
###########################################

MAIN = zz-mc6

SHRK = shr_kind_mod
SHRC = shr_const_mod
SHRG = ppgrid
SHRP = physconst

MPAR = parrrtm
MCLD = rrlw_cld
MCON = rrlw_con
MK01 = rrlw_kg01
MK02 = rrlw_kg02
MK03 = rrlw_kg03
MK04 = rrlw_kg04
MK05 = rrlw_kg05
MK06 = rrlw_kg06
MK07 = rrlw_kg07
MK08 = rrlw_kg08
MK09 = rrlw_kg09
MK10 = rrlw_kg10
MK11 = rrlw_kg11
MK12 = rrlw_kg12
MK13 = rrlw_kg13
MK14 = rrlw_kg14
MK15 = rrlw_kg15
MK16 = rrlw_kg16
MREF = rrlw_ref
MTBL = rrlw_tbl
MVSN = rrlw_vsn
MWVN = rrlw_wvn

CRTM = rrtmg_lw_rad
CCLD = rrtmg_lw_cldprop
CCMC = rrtmg_lw_cldprmc
CINI = rrtmg_lw_init
CKGS = rrtmg_lw_k_g
CRTR = rrtmg_lw_rtrn
CRMR = rrtmg_lw_rtrnmr
CRMC = rrtmg_lw_rtrnmc
CRM2 = rrtmg_lw_rtr2function
CSET = rrtmg_lw_setcoef
CTAU = rrtmg_lw_taumol

OBJMAN = $(MAIN).o

OBJSHR = $(SHRK).o \
         $(SHRC).o \
         $(SHRG).o \
         $(SHRP).o 

OBJPTH =  \
         $(CCLD).o \
         $(CCMC).o \
         $(CKGS).o \
         $(CRTR).o \
         $(CRMR).o \
         $(CRMC).o \
         $(CRM2).o \
         $(CSET).o \
         $(CTAU).o \
         $(CINI).o \
         $(CRTM).o 

OBMPTH =  \
         $(MPAR).o \
         $(MCLD).o \
         $(MCON).o \
         $(MK01).o \
         $(MK02).o \
         $(MK03).o \
         $(MK04).o \
         $(MK05).o \
         $(MK06).o \
         $(MK07).o \
         $(MK08).o \
         $(MK09).o \
         $(MK10).o \
         $(MK11).o \
         $(MK12).o \
         $(MK13).o \
         $(MK14).o \
         $(MK15).o \
         $(MK16).o \
         $(MREF).o \
         $(MTBL).o \
         $(MVSN).o \
         $(MWVN).o

MODSHR = $(SHRK).mod \
         $(SHRC).mod \
         $(SHRG).mod \
         $(SHRP).mod 

MODPTH1 =  \
         $(CCLD).mod \
         $(CCMC).mod \
         $(CRTR).mod \
         $(CRMR).mod \
         $(CRMC).mod \
         $(CRM2).mod \
         $(CSET).mod \
         $(CTAU).mod \
         $(CINI).mod \
         $(CRTM).mod 

MODPTH2 =  \
         $(MPAR).mod \
         $(MCLD).mod \
         $(MCON).mod \
         $(MK01).mod \
         $(MK02).mod \
         $(MK03).mod \
         $(MK04).mod \
         $(MK05).mod \
         $(MK06).mod \
         $(MK07).mod \
         $(MK08).mod \
         $(MK09).mod \
         $(MK10).mod \
         $(MK11).mod \
         $(MK12).mod \
         $(MK13).mod \
         $(MK14).mod \
         $(MK15).mod \
         $(MK16).mod \
         $(MREF).mod \
         $(MTBL).mod \
         $(MVSN).mod \
         $(MWVN).mod


#################
# Compile modules
#################

$(EXENAME): $(OBJSHR) $(OBMPTH) $(OBJPTH) $(OBJMAN)
	$(FC) -o $(EXENAME) $(OBJSHR) $(OBMPTH) $(OBJPTH) $(OBJMAN)
#test111: $(OBJSHR) $(OBMPTH) $(OBJPTH)
#	$(FC) -o test111 $(OBJSHR) $(OBMPTH) $(OBJPTH)

#$(MPAR).o : $(MPAR).f90
#	$(FC) -O -c $(MPAR).f90

$(SHRK).o : $(SHRK).F90
	$(FC) $(FCOPTS) -c $(SHRK).F90

$(SHRC).o : $(SHRC).F90
	$(FC) $(FCOPTS) -c $(SHRC).F90

$(SHRG).o : $(SHRG).F90
	$(FC) $(FCOPTS) -c $(SHRG).F90

$(SHRP).o : $(SHRP).F90
	$(FC) $(FCOPTS) -c $(SHRP).F90

$(MPAR).o : $(MPAR).f90
	$(FC) $(FCOPTS) -c $(MPAR).f90

$(MCLD).o : $(MCLD).f90
	$(FC) $(FCOPTS) -c $(MCLD).f90

$(MCON).o : $(MCON).f90
	$(FC) $(FCOPTS) -c $(MCON).f90

$(MK01).o : $(MK01).f90
	$(FC) $(FCOPTS)  -c $(MK01).f90


$(MK02).o : $(MK02).f90
	$(FC) $(FCOPTS)  -c $(MK02).f90


$(MK03).o : $(MK03).f90
	$(FC) $(FCOPTS)  -c $(MK03).f90


$(MK04).o : $(MK04).f90
	$(FC) $(FCOPTS)  -c $(MK04).f90


$(MK05).o : $(MK05).f90
	$(FC) $(FCOPTS)  -c $(MK05).f90


$(MK06).o : $(MK06).f90
	$(FC) $(FCOPTS)  -c $(MK06).f90


$(MK07).o : $(MK07).f90
	$(FC) $(FCOPTS)  -c $(MK07).f90


$(MK08).o : $(MK08).f90
	$(FC) $(FCOPTS)  -c $(MK08).f90


$(MK09).o : $(MK09).f90
	$(FC) $(FCOPTS)  -c $(MK09).f90


$(MK10).o : $(MK10).f90
	$(FC) $(FCOPTS)  -c $(MK10).f90


$(MK11).o : $(MK11).f90
	$(FC) $(FCOPTS)  -c $(MK11).f90


$(MK12).o : $(MK12).f90
	$(FC) $(FCOPTS)  -c $(MK12).f90


$(MK13).o : $(MK13).f90
	$(FC) $(FCOPTS)  -c $(MK13).f90


$(MK14).o : $(MK14).f90
	$(FC) $(FCOPTS)  -c $(MK14).f90


$(MK15).o : $(MK15).f90
	$(FC) $(FCOPTS)  -c $(MK15).f90


$(MK16).o : $(MK16).f90
	$(FC) $(FCOPTS)  -c $(MK16).f90


$(MREF).o : $(MREF).f90
	$(FC) $(FCOPTS)  -c $(MREF).f90


$(MTBL).o : $(MTBL).f90
	$(FC) $(FCOPTS)  -c $(MTBL).f90


$(MVSN).o : $(MVSN).f90
	$(FC) $(FCOPTS)  -c $(MVSN).f90


$(MWVN).o : $(MWVN).f90
	$(FC) $(FCOPTS)  -c $(MWVN).f90


#################
# Compile source
##################

$(CCLD).o : $(CCLD).f90
	$(FC) $(FCOPTS)  -c $(CCLD).f90


$(CCMC).o : $(CCMC).f90
	$(FC) $(FCOPTS)  -c $(CCMC).f90


$(CKGS).o : $(CKGS).f90
	$(FC) $(FCOPTS)  -c $(CKGS).f90


$(CRTR).o : $(CRTR).f90
	$(FC) $(FCOPTS)  -c $(CRTR).f90


$(CRMR).o : $(CRMR).f90
	$(FC) $(FCOPTS)  -c $(CRMR).f90


$(CRMC).o : $(CRMC).f90
	$(FC) $(FCOPTS)  -c $(CRMC).f90

$(CRM2).o : $(CRM2).f90
	$(FC) $(FCOPTS)  -c $(CRM2).f90

$(CSET).o : $(CSET).f90
	$(FC) $(FCOPTS)  -c $(CSET).f90


$(CTAU).o : $(CTAU).f90
	$(FC) $(FCOPTS)  -c $(CTAU).f90


$(CINI).o : $(CINI).f90
	$(FC) $(FCOPTS)  -c $(CINI).f90


$(CRTM).o : $(CRTM).f90
	$(FC) $(FCOPTS)  -c $(CRTM).f90

$(MAIN).o : $(MAIN).f90
	$(FC) $(FCOPTS)  -c $(MAIN).f90


#########
# clean
#########

clean:
	rm $(EXENAME) $(OBJSHR) $(OBJPTH) $(OBMPTH) $(MODSHR) $(MODPTH1) $(MODPTH2) $(OBLMAN) 

###########################################
# Create object file directory if necessary
###########################################

#makedir :
#	if [ ! -d $(CDOTO) ] ; \
#	then \
#	   mkdir $(CDOTO) ; \
#	fi

