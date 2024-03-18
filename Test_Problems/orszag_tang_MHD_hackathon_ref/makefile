# *********************************************************
# PLUTO 4.4-patch3  Makefile  
#
# Last modified: May 22, 2023
# *********************************************************

plutino:                              # Default target

#ARCH         = Linux.gcc.defs

CC      = nvc++
CFLAGS  = -c -O3 -cudalib=nvtx3 -D_USE_NVTX -acc -Minfo=accel
LDFLAGS = -lm -cudalib=nvtx3 -acc

#PLUTINO_DIR    = /home/truzzi/gplutino
SRC          = $(PLUTINO_DIR)/Src
INCLUDE_DIRS = -I. -I$(SRC)
VPATH        = ./:$(SRC) 

#include $(PLUTINO_DIR)/Config/$(ARCH)

# ---------------------------------------------------------
#         Set headers and object files 
# ---------------------------------------------------------

HEADERS = plutino.hpp prototypes.hpp macro.hpp definition.hpp structs.hpp
OBJ = main.o initial_conditions.o arrays.o boundary.o \
			cons_to_prim.o debug.o flux.o input_output.o \
			reconstruct.o riemann.o rungekutta.o set_indexes.o \
			source.o update.o wave_speed.o nvtx.o
			
#include $(SRC)/Math_Tools/makefile

# ---------------------------------------------------------
#  Define macros by adding -D<name> where <name> has been
#  set to TRUE in the system configuration file (.defs) 
# ---------------------------------------------------------

#ifeq ($(strip $(PARALLEL)), TRUE)
# CFLAGS += -I$(SRC)/Parallel -DPARALLEL
# include $(SRC)/Parallel/makefile
# ifeq ($(strip $(USE_ASYNC_IO)), TRUE)
#  CFLAGS += -DUSE_ASYNC_IO
# endif
#endif

#ifeq ($(strip $(USE_HDF5)), TRUE)
# CFLAGS += -DUSE_HDF5
# OBJ    += hdf5_io.o
#endif
      
#ifeq ($($strip $(USE_PNG)), TRUE)
# CFLAGS += -DUSE_PNG
#endif

-include local_make

# ---------------------------------------------------------
#   Additional_CFLAGS_here   ! dont change this line
# ---------------------------------------------------------


# ---------------------------------------------------------
#   Additional_header_files_here   ! dont change this line
# ---------------------------------------------------------


# ---------------------------------------------------------
#   Additional_object_files_here   ! dont change this line
# ---------------------------------------------------------

#OBJ += plm_states.o
#OBJ += vec_pot_diff.o
#OBJ += vec_pot_update.o
#OBJ += ctu_step.o
#OBJ += hancock.o
#include $(SRC)/MHD/makefile
#include $(SRC)/EOS/Ideal/makefile

# ---------------------------------------------------------
#    PLUTINO target rule
# ---------------------------------------------------------

plutino: $(OBJ) 
	$(CC) $(OBJ) $(LDFLAGS) -o $@

# ---------------------------------------------------------
#                    Suffix rule
# ---------------------------------------------------------

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDE_DIRS) $<

clean:
	@rm -f	*.o
	@echo make clean: done

# ---------------------------------------------------------
#          Dependencies for object files
# ---------------------------------------------------------

$(OBJ):  $(HEADERS)

