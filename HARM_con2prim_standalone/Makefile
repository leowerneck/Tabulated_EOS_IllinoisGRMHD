SRC       = standalone.c EOS_functions.c nuc_eos_press.c readtable.c utoprim_3d_eos.c utoprim_2d_eos.c utoprim_2d_eos_safe_guess.c utoprim_2d_eos_backup.c
INC       = decs.h driver.h globals.h helpers.h macros.h nuc_eos.h options.h u2p_defs.h units.h
OBJ       = $(SRC:.c=.o)
EXE       = standalone

CC        = gcc
CFLAGS    = -march=native -Ofast
CLIBS     =

HDF5DIR   = /usr/local/hdf5
HDF5INC   = $(HDF5DIR)/include
HDF5LIB   = $(HDF5DIR)/lib

EXTRAINCS = -I$(HDF5INC) -L$(HDF5LIB)
EXTRALIBS = -lm -lhdf5 -lz

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(CFLAGS) $(EXTRAINCS) $(OBJ) -o $@ $(EXTRALIBS)

$(OBJ): %.o : %.c $(INC)
	$(CC) $(CFLAGS) $(EXTRAINCS) -c $< -o $@

clean:
	rm -rf $(EXE) $(OBJ) *.asc
