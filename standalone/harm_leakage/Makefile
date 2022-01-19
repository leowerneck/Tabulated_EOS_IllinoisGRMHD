# HDF5 directory
# HDF_DIR = /usr/local/opt/hdf5
HDF_DIR = /usr/lib/x86_64-linux-gnu/hdf5/serial
NRPYEOS_DIR = ../NRPyEOS/

CC      = gcc
CFLAGS  = -Wall -I../ -I../NRPyEOS -I$(HDF_DIR)/include
LDFLAGS = -L$(HDF_DIR)/lib -lhdf5 -lm

all: harm_neutrinos_standalone

NRPyEOS_readtable_set_EOS_params.o: $(NRPYEOS_DIR)/NRPyEOS_readtable_set_EOS_params.c
	$(CC) $(CFLAGS) -c $< -o $@

NRPyEOS_free_memory.o: $(NRPYEOS_DIR)/NRPyEOS_free_memory.c
	$(CC) $(CFLAGS) -c $< -o $@

NRPyEOS_Tabulated_general_interpolators.o: $(NRPYEOS_DIR)/NRPyEOS_Tabulated_general_interpolators.c
	$(CC) $(CFLAGS) -c $< -o $@

NRPyEOS_Tabulated_known_T.o: $(NRPYEOS_DIR)/NRPyEOS_Tabulated_known_T.c
	$(CC) $(CFLAGS) -c $< -o $@

NRPyEOS_Tabulated_unknown_T.o: $(NRPYEOS_DIR)/NRPyEOS_Tabulated_unknown_T.c
	$(CC) $(CFLAGS) -c $< -o $@

NRPyEOS_validation.o: $(NRPYEOS_DIR)/NRPyEOS_validation.c
	$(CC) $(CFLAGS) -c $< -o $@

harm_neutrinos.o: harm_neutrinos.c harm_neutrinos.h
	$(CC) $(CFLAGS) -c $< -o $@

harm_neutrinos_standalone.o: harm_neutrinos_standalone.c
	$(CC) $(CFLAGS) -c $< -o $@

harm_neutrinos_standalone: harm_neutrinos_standalone.o harm_neutrinos.o NRPyEOS_readtable_set_EOS_params.o NRPyEOS_free_memory.o NRPyEOS_Tabulated_general_interpolators.o NRPyEOS_Tabulated_known_T.o NRPyEOS_Tabulated_unknown_T.o
	$(CC) $(CFLAGS) harm_neutrinos_standalone.o harm_neutrinos.o NRPyEOS_readtable_set_EOS_params.o NRPyEOS_free_memory.o NRPyEOS_Tabulated_general_interpolators.o NRPyEOS_Tabulated_known_T.o NRPyEOS_Tabulated_unknown_T.o -o harm_neutrinos_standalone $(LDFLAGS)


clean:
	rm -f *.o harm_neutrinos_standalone

veryclean: clean
	rm -f *.txt
