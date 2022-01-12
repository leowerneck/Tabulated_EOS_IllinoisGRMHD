# HDF5 directory
# HDF_DIR = /usr/local/opt/hdf5
HDF_DIR = /usr/lib/x86_64-linux-gnu/hdf5/serial

CC      = gcc
CFLAGS  = -Wall -I../ -I$(HDF_DIR)/include
LDFLAGS = -L$(HDF_DIR)/lib -lhdf5 -lm

all: NRPyEOS_validation

NRPyEOS_readtable_set_EOS_params.o: NRPyEOS_readtable_set_EOS_params.c
	$(CC) $(CFLAGS) -c $< -o $@

NRPyEOS_free_memory.o: NRPyEOS_free_memory.c
	$(CC) $(CFLAGS) -c $< -o $@

NRPyEOS_Tabulated_general_interpolators.o: NRPyEOS_Tabulated_general_interpolators.c
	$(CC) $(CFLAGS) -c $< -o $@

NRPyEOS_Tabulated_known_T.o: NRPyEOS_Tabulated_known_T.c
	$(CC) $(CFLAGS) -c $< -o $@

NRPyEOS_Tabulated_unknown_T.o: NRPyEOS_Tabulated_unknown_T.c
	$(CC) $(CFLAGS) -c $< -o $@

NRPyEOS_validation.o: NRPyEOS_validation.c
	$(CC) $(CFLAGS) -c $< -o $@

NRPyEOS_validation: NRPyEOS_readtable_set_EOS_params.o NRPyEOS_free_memory.o NRPyEOS_Tabulated_general_interpolators.o NRPyEOS_Tabulated_known_T.o NRPyEOS_Tabulated_unknown_T.o NRPyEOS_validation.o
	$(CC) $(CFLAGS) NRPyEOS_readtable_set_EOS_params.o NRPyEOS_free_memory.o NRPyEOS_Tabulated_general_interpolators.o NRPyEOS_Tabulated_known_T.o NRPyEOS_Tabulated_unknown_T.o NRPyEOS_validation.o -o NRPyEOS_validation $(LDFLAGS)

clean:
	rm -f *.o NRPyEOS_validation

veryclean: clean
	rm -f *.txt
