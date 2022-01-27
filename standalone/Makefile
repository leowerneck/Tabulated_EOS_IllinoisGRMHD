BASEDIR = .
NRPyEOS_DIR = $(BASEDIR)/NRPyEOS
LEAKAGE_DIR = $(BASEDIR)/NRPyLeakage
TESTS_DIR   = $(BASEDIR)/Leakage_tests
HARM_DIR    = $(BASEDIR)/harm_leakage

HDF5_DIR = /usr/local/opt/hdf5
# HDF5_DIR = /usr/lib/x86_64-linux-gnu/hdf5/serial

CC       = gcc
CFLAGS   = -Wall -O2 -fopenmp -march=native -I$(BASEDIR) -I$(NRPyEOS_DIR) -I$(LEAKAGE_DIR) -I$(HARM_DIR) -I$(TESTS_DIR) -I$(HDF5_DIR)/include

CXX      = g++
CXXFLAGS = -Wall -O2 -fopenmp -march=native -I$(BASEDIR) -I$(NRPyEOS_DIR) -I$(LEAKAGE_DIR) -I$(HARM_DIR) -I$(TESTS_DIR) -I$(HDF5_DIR)/include

LD_FLAGS = -lm -L$(HDF5_DIR)/lib -lhdf5

SRC  = standalone.c
EXE  = $(SRC:.c=)

NRPyEOS_SRC = NRPyEOS_Tabulated_general_interpolators.c \
              NRPyEOS_Tabulated_known_T.c \
              NRPyEOS_Tabulated_unknown_T.c \
              NRPyEOS_readtable_set_EOS_params.c \
              NRPyEOS_free_memory.c

NRPyEOS_INC = $(NRPyEOS_DIR)/NRPyEOS_Tabulated_headers.h \
              $(NRPyEOS_DIR)/NRPyEOS_Tabulated_helpers.h

LEAKAGE_SRC = NRPyLeakage_Fermi_Dirac_integrals.c \
              NRPyLeakage_compute_GRMHD_source_terms_nrpy_constants.c \
              NRPyLeakage_compute_GRMHD_source_terms_harm_constants.c \
              NRPyLeakage_compute_GRMHD_source_terms.c \
              NRPyLeakage_compute_optical_depths.c

LEAKAGE_INC = $(LEAKAGE_DIR)/NRPyLeakage.h

TESTS_SRC   = OpticallyThinGas_NRPyLeakage.c OpticallyThinGas_harm_leakage.c 
TESTS_INC   = $(TESTS_DIR)/Leakage_tests.h

HARM_SRC    = harm_neutrinos.c #harm_optical_depth.c
HARM_INC    = $(HARM_DIR)/harm_neutrinos.h $(HARM_DIR)/harm_units.h 

OBJ = $(addprefix obj/,$(SRC:.c=.o)) \
      $(addprefix obj/,$(NRPyEOS_SRC:.c=.o)) \
      $(addprefix obj/,$(LEAKAGE_SRC:.c=.o)) \
      $(addprefix obj/,$(TESTS_SRC:.c=.o)) \
      $(addprefix obj/,$(HARM_SRC:.c=.o))

INC = Basic_defines.h $(NRPyEOS_INC) $(LEAKAGE_INC) $(TESTS_INC) $(HARM_INC)

all: obj_dir $(EXE)

obj_dir:
	mkdir -p obj

$(EXE): $(OBJ) $(COBJ)
	$(CC) $(CFLAGS) $(OBJ) $(COBJ) -o $(EXE) $(LD_FLAGS)

obj/standalone.o: standalone.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/NRPyEOS_Tabulated_general_interpolators.o: $(NRPyEOS_DIR)/NRPyEOS_Tabulated_general_interpolators.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/NRPyEOS_Tabulated_known_T.o: $(NRPyEOS_DIR)/NRPyEOS_Tabulated_known_T.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/NRPyEOS_Tabulated_unknown_T.o: $(NRPyEOS_DIR)/NRPyEOS_Tabulated_unknown_T.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/NRPyEOS_readtable_set_EOS_params.o: $(NRPyEOS_DIR)/NRPyEOS_readtable_set_EOS_params.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/NRPyEOS_free_memory.o: $(NRPyEOS_DIR)/NRPyEOS_free_memory.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/NRPyLeakage_Fermi_Dirac_integrals.o: $(LEAKAGE_DIR)/NRPyLeakage_Fermi_Dirac_integrals.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/NRPyLeakage_compute_GRMHD_source_terms_nrpy_constants.o: $(LEAKAGE_DIR)/NRPyLeakage_compute_GRMHD_source_terms_nrpy_constants.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/NRPyLeakage_compute_GRMHD_source_terms_harm_constants.o: $(LEAKAGE_DIR)/NRPyLeakage_compute_GRMHD_source_terms_harm_constants.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/NRPyLeakage_compute_GRMHD_source_terms.o: $(LEAKAGE_DIR)/NRPyLeakage_compute_GRMHD_source_terms.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/NRPyLeakage_compute_optical_depths.o: $(LEAKAGE_DIR)/NRPyLeakage_compute_optical_depths.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/OpticallyThinGas_NRPyLeakage.o: $(TESTS_DIR)/optically_thin/OpticallyThinGas_NRPyLeakage.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/OpticallyThinGas_harm_leakage.o: $(TESTS_DIR)/optically_thin/OpticallyThinGas_harm_leakage.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/harm_neutrinos.o: $(HARM_DIR)/harm_neutrinos.c
	$(CC) $(CFLAGS) -c $< -o $@

obj/harm_optical_depth.o: $(HARM_DIR)/harm_optical_depth.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(EXE) $(OBJ) $(COBJ) obj

veryclean: clean
	rm -rf *.txt *.png
