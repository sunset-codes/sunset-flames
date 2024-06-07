# Makefile for sunset_code
#
#
# ========================== OPTIONS ==============================================================
# -------------------------------------------------------------------------------------------------
# mpi        Shared only (0) or distributed-shared (1) acceleration                    (default: 1)          
# dim3       Two (0) or three (1) dimensional simulation                               (default: 0)
# react      Reacting (1) or inert (0) flow                                            (default: 1)
# thermo     Isothermal (0) or thermal (1) flow                                        (default: 1)
# restart    Start from initial conditions (0) or restart file (1)                     (default: 0)
# pgrad      Drive the flow with a pressure gradient and P.I.D control                 (default: 0)
# flout      Output the flame structure (1) or don't (0)                               (default: 0)
# allout     If 3D, output the entire domain (1) or just a slice (0)                   (default: 1)
# morder     m (order) value = 4,6,8,10                                                (default: 8)
# tarout     Compress output files (1) or don't (0)                                    (default: 0)
# -------------------------------------------------------------------------------------------------
#
# EXAMPLE USAGE:
#
# Default make for reacting 2D flows:
# make

# For isothermal flows:
# make thermo=0 dim3=X mpi=X pgrad=X
#
# For standard combustion problems, react=1 will overwrite any thermo flags:
# make react=1 dim3=X mpi=X flout=X      <---------- standard combustion make
#
# For thermal flows with real gas properties:
# make thermo=1 dim3=X mpi=X pgrad=X 
#

#
# Choose compiler depending on whether mpi
ifneq ($(mpi),0)
FC := mpif90
LD := mpif90
else
FC := gfortran
LD := gfortran
endif

# Set compiler flags based on make options
CFLAGS := -Wall -O3 -g -m64
FFLAGS := -fopenmp -fbounds-check -ffpe-trap=zero -O3 -Wall -g -J./obj -I./obj -m64

# Order of the numerical scheme (even, from 4 to 10, default 8)
ifeq ($(morder),4)
FFLAGS += -Dorder=$(morder)
else
ifeq ($(morder),6)
FFLAGS += -Dorder=$(morder)
else
ifeq ($(morder),10)
FFLAGS += -Dorder=$(morder)
else
FFLAGS += -Dorder=8
endif
endif
endif

# Isothermal or not. isoT can only happen if not(react)
ifeq ($(thermo), 0)
ifeq ($(react), 0)
FFLAGS += -DisoT
endif
endif

# Reacting? set react if so
ifneq ($(react), 0)
FFLAGS += -Dreact
endif

# Restart from dump file.
ifeq ($(restart), 1)
FFLAGS += -Drestart
endif

# Compress output files
ifeq ($(tarout), 1)
FFLAGS += -Dtarout
endif

# Multiprocessor? (use mpi?)
ifneq ($(mpi),0)
FFLAGS += -Dmp
endif

# Three dimensional?
ifeq ($(dim3),1)
FFLAGS += -Ddim3
endif

# Flow driven by pressure gradient?
ifeq ($(pgrad),1)
FFLAGS += -Dpgrad
endif

# Output full chemical composition?
ifeq ($(flout),1)
FFLAGS += -Doutput_flame
endif

# Output entire domain?
ifneq ($(allout),0)
FFLAGS += -Dallout
endif

LDFLAGS := -fopenmp -m64
###### -lopenblas 

# Identify directories
SUB_DIRS := common base
SRC_DIR  := $(addprefix source/,$(SUB_DIRS))

# identify object files
#parameters come first, as almost everything depends on them.
OBJ_FILES := obj/kind_parameters.o obj/common_parameter.o obj/common_vars.o
OBJ_FILES += obj/rbfs.o obj/mirror_boundaries.o obj/derivatives.o 
OBJ_FILES += obj/mpi_transfers.o obj/transport.o obj/thermodynamics.o
OBJ_FILES += obj/neighbours.o obj/output.o obj/statistics.o 
OBJ_FILES += obj/turbulence.o obj/svdlib.o
OBJ_FILES += obj/load_data.o obj/setup_domain.o obj/setup_flow.o
OBJ_FILES += obj/labf.o obj/fd.o obj/chemistry.o
OBJ_FILES += obj/characteristic_boundaries.o obj/rhs.o
OBJ_FILES += obj/step.o
OBJ_FILES += $(foreach sdir,$(SRC_DIR),$(patsubst $(sdir)/%.F90,obj/%.o,$(wildcard $(sdir)/*.F90)))

vpath %.F90 $(SRC_DIR)

#-------
default: sunset
sunset: $(OBJ_FILES)
	$(LD) -o $@ $^ $(LDFLAGS)

obj/%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -vf ./obj/*.o
	rm -vf ./obj/*.mod
	rm -vf ./sunset
	rm -rfv fort.*	
	rm -vf ./data_out/fields*
	rm -vf ./data_out/nodes*
	rm -vf ./data_out/flame*
	rm -vf ./data_out/time.out
	rm -vf ./data_out/statistics/*.out
	rm -vf ./paraview_files/LAYER*
	rm -vf ./data_out/*.tar.gz

