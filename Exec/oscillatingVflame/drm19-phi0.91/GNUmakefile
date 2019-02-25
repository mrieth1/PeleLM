#TOP = ../../../..
#AMREX_HOME ?= ${TOP}/amrex
#IAMR_HOME ?= ${TOP}/IAMR
#PELELM_HOME ?= ${TOP}/PeleLM

#
# Build configuration
#

DIM             = 2
COMP            = gnu
DEBUG           = FALSE
USE_MPI         = TRUE
USE_OMP         = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE

USE_CUDA        = FALSE

# Enable flctin namelist
USE_FLCTS       = TRUE

USE_SUNDIALS_3x4x=TRUE
CVODE_LIB_DIR = /home/amfelden/CODES/PeleLM/Source/ChemDriver/cvode-4.0.2/instdir/lib/

USE_KLU = TRUE
ifeq ($(USE_KLU), TRUE)
    DEFINES  += -DUSE_KLU
endif

# Enable ht.velocity_plotfile option
#USE_VELOCITY    = TRUE

# Enable the restart-from-plotfile option (ht.pltfile)
#  ... default this off because when on, code expects it to be used
# USE_NEWMECH    = TRUE

# Enable time-dependent forcing term
#USE_TDEPFORCING = TRUE

# Choose model (from list below), and pmf file
Chemistry_Model = DRM19
fEXE_sources += drm19-phi0.91.f

#Chemistry_Model = LIDRYER
#fEXE_sources += LiDryer_040_01ATM.f

#Chemistry_Model = LUDME
#fEXE_sources += LuDME_0700.f

#Chemistry_Model = GRI30
#fEXE_sources += gri30_070.f

Blocs   := .. .

CEXE_sources += actual_Creactor.cpp  actual_Creactor_unit.cpp 
FEXE_sources += Prob_$(DIM)d.F Prob_F.F
CEXE_headers += 
FEXE_headers += probdata.H Prob_F.H

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
include Make.CVODE
# Available chemistry models:
# CanuH2 chem-CH4-2step chem-H dme glar gri Hai H-CW header inert Konnov
# LiDryer Marinov prf_ethanol Roehl sandiego smooke usc
