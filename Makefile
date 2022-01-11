# Makefile 
# Nom du compilateur
FC = f95

# Options de compilation: optimisation, debug etc...
OPT = -fdefault-real-8 -fbounds-check   
# nom de l'executable
EXE = rayleigh
# Options de l'edition de lien..
LINKOPT =  

# Defining the objects (OBJS) variables

OBJS =  \
       Initialisation.o \
       main.o \
       matgen.o \
       Solveur.o \
       reprise.o \
       calc_dt.o \
       calc_forces_volumiques.o \
       calc_visc.o \
       calc_vort_uv_centre.o \
       conv_uv.o \
       bc_uv.o \
       bc_temp.o \
       calc_visc_temp.o \
       conv_temp.o \
       ensight.o 

# Linking object files
exe :   $(OBJS)
	$(FC) $(LINKOPT) $(MODS) $(OBJS)  -o $(EXE) 
	
Initialisation.o : Initialisation.f90
	$(FC) -c $(OPT)  Initialisation.f90

main.o : main.f90
	$(FC) -c $(OPT)  main.f90

matgen.o : matgen.f90
	$(FC) -c $(OPT) matgen.f90

Solveur.o : Solveur.f90
	$(FC) -c $(OPT) Solveur.f90
	
reprise.o : reprise.f90
	$(FC) -c $(OPT) reprise.f90
	
calc_dt.o : calc_dt.f90
	$(FC) -c $(OPT) calc_dt.f90
	
calc_forces_volumiques.o : calc_forces_volumiques.f90
	$(FC) -c $(OPT) calc_forces_volumiques.f90
	
calc_visc.o : calc_visc.f90
	$(FC) -c $(OPT) calc_visc.f90
	
calc_vort_uv_centre.o : calc_vort_uv_centre.f90
	$(FC) -c $(OPT) calc_vort_uv_centre.f90
	
conv_uv.o : conv_uv.f90
	$(FC) -c $(OPT) conv_uv.f90
	
bc_uv.o : bc_uv.f90
	$(FC) -c $(OPT) bc_uv.f90

bc_temp.o : bc_temp.f90
	$(FC) -c $(OPT) bc_temp.f90

calc_visc_temp.o : calc_visc_temp.f90
	$(FC) -c $(OPT) calc_visc_temp.f90

conv_temp.o : conv_temp.f90
	$(FC) -c $(OPT) conv_temp.f90

ensight.o : ensight.f90
	$(FC) -c $(OPT) ensight.f90
	

# Removing object files
clean :
	/bin/rm -f $(OBJS) $(EXE)  *.mod

config :
	if [ ! -d obj ] ; then mkdir obj ; fi
	if [ ! -d run ] ; then mkdir run ; fi

.SUFFIXES: .f90

