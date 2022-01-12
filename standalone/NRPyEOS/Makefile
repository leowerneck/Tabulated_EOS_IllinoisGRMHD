CC     = gcc
CFLAGS = -Wall -I../ -I/usr/lib/x86_64-linux-gnu/hdf5/serial/include

all: NRPyEOS_free_memory.o NRPyEOS_readtable.o NRPyEOS_Tabulated_general_interpolators.o NRPyEOS_Tabulated_known_T.o

NRPyEOS_readtable.o: NRPyEOS_readtable.c
	$(CC) $(CFLAGS) -c $< -o $@

NRPyEOS_free_memory.o: NRPyEOS_free_memory.c
	$(CC) $(CFLAGS) -c $< -o $@

NRPyEOS_Tabulated_general_interpolators.o: NRPyEOS_Tabulated_general_interpolators.c
	$(CC) $(CFLAGS) -c $< -o $@

NRPyEOS_Tabulated_known_T.o: NRPyEOS_Tabulated_known_T.c
	$(CC) $(CFLAGS) -c $< -o $@

NRPyEOS_Tabulated_unknown_T.o: NRPyEOS_Tabulated_unknown_T.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o
