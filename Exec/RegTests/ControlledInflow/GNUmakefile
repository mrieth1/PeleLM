TOP = ../../../..
AMREX_HOME ?= ${TOP}/amrex
IAMR_HOME ?= ${TOP}/IAMR
PELELM_HOME ?= ${TOP}/PeleLM

#
# Build configuration
#

DIM             = 2
COMP            = gnu
DEBUG           = FALSE
USE_MPI         = TRUE
USE_OMP         = TRUE
PRECISION       = DOUBLE
VERBOSE         = FALSE

TINY_PROFILE = FALSE
#USE_WBAR = TRUE

DEFINES+=-DALWAYS_NEW_J

# Enable ht.velocity_plotfile option
#USE_VELOCITY    = FALSE

# Enable the restart-from-plotfile option (ht.pltfile)
#  ... default this off because when on, code expects it to be used
# USE_NEWMECH    = FALSE


# Choose model (from list below), and pmf file
#Chemistry_Model = drm19
#fEXE_sources += drm19Soln_seed_0.50.f

Chemistry_Model = LiDryer
fEXE_sources += LiDryer_040_01ATM.f

#Chemistry_Model = LUDME
#fEXE_sources += LuDME_0700.f

#Chemistry_Model = GRI30
#fEXE_sources += gri30_070.f

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := EGLib

Blocs   := .

CEXE_sources +=
F90EXE_sources += Prob_$(DIM)d.F90 probdata.F90
CEXE_headers += 
FEXE_headers += 

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
# Available chemistry models:
# CanuH2 chem-CH4-2step chem-H dme glar gri Hai H-CW header inert Konnov
# LiDryer Marinov prf_ethanol Roehl sandiego smooke usc
