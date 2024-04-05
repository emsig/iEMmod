# Makefile

include ../Make_include
#OPTC:=
#OPTF:=
#PGI
#LIBS    = -L$L -lgenfft -lpgf90 -lpgf90_rpm1 -lpgf90rtl -lpgf902 -lpghpf -lpgftnrtl -lm
#Intel
#LIBS    = -L$L -lgenfft -lifcore -lm
#pathscale
#LIBS    = -L$L -lgenfft -lm
# GNU Cray
#LIBS    = -L$L -lgenfft -lgfortran 
#AMD openf95
#LIBS    += -L$L -lgenfft -lacml_mv -lfortran -lm -lffio -lc -lrt

ALL: iemmod 

SRC	= iemmod.c \
		emmod.c \
		kxwmod.c \
		bessel.c \
		verbosepkg.c \
		docpkge.c \
                atopkge.c \
                getpars.c \
		wallclock_time.c

SRCF77 = zqk61n.f \
                j4save.f \
                fdump.f \
                i1mach.f \
                xermsg.f \
                xerprn.f \
                xersve.f \
                xerhlt.f \
                xercnt.f \
                xgetua.f

SRCF = Gamma.F90 \
		dGamma.F90 \
		Wprop.F90 \
		Wpropderiv.F90 \
		Rplusderiv.F90 \
		Rplus.F90 \
		Rminderiv.F90 \
		Rmin.F90 \
		Pupplus.F90 \
		Pupplusderiv.F90 \
		Pupmin.F90 \
		Pupminderiv.F90 \
		Pdownplus.F90 \
		Pdownplusderiv.F90 \
		Pdownmin.F90 \
		Pdownminderiv.F90 \
		Ptotalx.F90 \
		dPtotalx.F90 \
		Ptotaly.F90 \
		dPtotaly.F90 \
		Ptotalz.F90 \
		dPtotalz.F90 \
		Ptotalref.F90 \
		dPtotalref.F90 \
		Gin11.F90 \
		dGin11.F90 \
		Gin12.F90 \
		dGin12.F90 \
		Gin13.F90 \
		dGin13.F90 \
		Gin22.F90 \
		dGin22.F90 \
		Gin23.F90 \
		dGin23.F90 \
		Gin31.F90 \
		dGin31.F90 \
		Gin32.F90 \
		dGin32.F90 \
		Gin33.F90 \
		dGin33.F90 \
		getcoords.F90 \
		hankeltrans_xx.F90 \
		hankeltrans_yy.F90 \
		hankeltrans_zz.F90 \
		hankeltrans_ref.F90 \
		evalpoints.F90 \
		gridit_xx.F90 \
		dgridit_xx.F90 \
		gridit_yy.F90 \
		dgridit_yy.F90 \
		gridit_zz.F90 \
		dgridit_zz.F90 \
		gridit_ref.F90 \
		zqk61_setup_grid.F90

OBJ	= $(SRC:%.c=%.o)
OBJF	= $(SRCF:%.F90=%.o)
OBJF77	= $(SRCF77:%.f=%.o)

iemmod:	$(OBJ) $(OBJF) $(OBJF77)
	$(CC) $(OPTC) $(CFLAGS) -o iemmod $(OBJ) $(OBJF) $(OBJF77) $(LIBS)

install: iemmod
	cp iemmod $B

clean:
		rm -f core $(OBJ) $(OBJF) $(OBJF77) iemmod

realclean: clean
		rm -f $B/iemmod




