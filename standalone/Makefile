BASEDIR = .
NRPyEOS_DIR = $(BASEDIR)/NRPyEOS
LEAKAGE_DIR = $(BASEDIR)/NRPyLeakage
TESTS_DIR   = $(BASEDIR)/Leakage_tests
HARM_DIR    = $(BASEDIR)/harm_leakage

HDF5_DIR = /usr/local/opt/hdf5
# HDF5_DIR = /usr/lib/x86_64-linux-gnu/hdf5/serial

CC       = gcc
CFLAGS   = -Wall -O2 -fopenmp -march=native -I$(BASEDIR) -I$(NRPyEOS_DIR) -I$(LEAKAGE_DIR) -I$(HARM_DIR) -I$(TESTS_DIR) -I$(HDF5_DIR)/include

FC       = gfortran
FCFLAGS  = -Wall -O2 -fopenmp -march=native -J./obj

CXX      = g++
CXXFLAGS = -Wall -O2 -fopenmp -march=native -I$(BASEDIR) -I$(NRPyEOS_DIR) -I$(LEAKAGE_DIR) -I$(HARM_DIR) -I$(TESTS_DIR) -I$(HDF5_DIR)/include

LD_FLAGS = -lm -L$(HDF5_DIR)/lib -lhdf5 -lgfortran

SRC  = standalone.c
EXE  = $(SRC:.c=)

NRPyEOS_SRC = NRPyEOS_Tabulated_general_interpolators.c \
              NRPyEOS_Tabulated_known_T.c \
              NRPyEOS_Tabulated_unknown_T.c \
              NRPyEOS_readtable_set_EOS_params.c \
              NRPyEOS_free_memory.c \
              NRPyEOS_constant_entropy_slice.c \
	      NRPyEOS_module.F90 \
              NRPyEOS_fortran_interface.F90

NRPyEOS_INC = $(NRPyEOS_DIR)/NRPyEOS_Tabulated_headers.h \
              $(NRPyEOS_DIR)/NRPyEOS_Tabulated_helpers.h

LEAKAGE_SRC = NRPyLeakage_Fermi_Dirac_integrals.c \
              NRPyLeakage_compute_GRMHD_source_terms_and_opacities_nrpy_constants.c \
              NRPyLeakage_compute_GRMHD_source_terms_and_opacities_harm_constants.c \
              NRPyLeakage_compute_GRMHD_source_terms_and_opacities.c \
              NRPyLeakage_compute_optical_depths.c \
              NRPyLeakage_compute_opacities_nrpy_constants.c \
              NRPyLeakage_compute_opacities_harm_constants.c \
              NRPyLeakage_compute_opacities.c

LEAKAGE_INC = $(LEAKAGE_DIR)/NRPyLeakage.h

ZELMANI_SRC = ZelmaniLeak_calc_taus.F90 ZelmaniLeak_calc_leak.F90

TESTS_SRC   = OpticallyThinGas_NRPyLeakage.c \
              OpticallyThinGas_harm_leakage.c \
              ConstantDensitySphere_NRPyLeakage.c \
              ConstantDensitySphere_harm_leakage.c \
              NRPyEOS_benchmark.c

TESTS_INC   = $(TESTS_DIR)/Leakage_tests.h

HARM_SRC    = harm_neutrinos.c #harm_optical_depth.c
HARM_INC    = $(HARM_DIR)/harm_neutrinos.h $(HARM_DIR)/harm_units.h

OBJ = $(addprefix obj/,$(SRC:.c=.o)) \
      $(addprefix obj/,$(patsubst %.F90,%.o, $(patsubst %.c,%.o,$(NRPyEOS_SRC)))) \
      $(addprefix obj/,$(LEAKAGE_SRC:.c=.o)) \
      $(addprefix obj/,$(TESTS_SRC:.c=.o)) \
      $(addprefix obj/,$(HARM_SRC:.c=.o)) \
      $(addprefix obj/,$(ZELMANI_SRC:.F90=.o))

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

obj/NRPyEOS_constant_entropy_slice.o: $(NRPyEOS_DIR)/NRPyEOS_constant_entropy_slice.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/NRPyEOS_module.o: $(NRPyEOS_DIR)/NRPyEOS_module.F90
	$(FC) $(FCFLAGS) -c $< -o $@

obj/NRPyEOS_fortran_interface.o: $(NRPyEOS_DIR)/NRPyEOS_fortran_interface.F90
	$(FC) $(FCFLAGS) -c $< -o $@

obj/NRPyLeakage_Fermi_Dirac_integrals.o: $(LEAKAGE_DIR)/NRPyLeakage_Fermi_Dirac_integrals.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/NRPyLeakage_compute_GRMHD_source_terms_and_opacities_nrpy_constants.o: $(LEAKAGE_DIR)/NRPyLeakage_compute_GRMHD_source_terms_and_opacities_nrpy_constants.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/NRPyLeakage_compute_GRMHD_source_terms_and_opacities_harm_constants.o: $(LEAKAGE_DIR)/NRPyLeakage_compute_GRMHD_source_terms_and_opacities_harm_constants.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/NRPyLeakage_compute_GRMHD_source_terms_and_opacities.o: $(LEAKAGE_DIR)/NRPyLeakage_compute_GRMHD_source_terms_and_opacities.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/NRPyLeakage_compute_opacities_nrpy_constants.o: $(LEAKAGE_DIR)/NRPyLeakage_compute_opacities_nrpy_constants.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/NRPyLeakage_compute_opacities_harm_constants.o: $(LEAKAGE_DIR)/NRPyLeakage_compute_opacities_harm_constants.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/NRPyLeakage_compute_opacities.o: $(LEAKAGE_DIR)/NRPyLeakage_compute_opacities.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/NRPyLeakage_compute_optical_depths.o: $(LEAKAGE_DIR)/NRPyLeakage_compute_optical_depths.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/OpticallyThinGas_NRPyLeakage.o: $(TESTS_DIR)/optically_thin/OpticallyThinGas_NRPyLeakage.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/OpticallyThinGas_harm_leakage.o: $(TESTS_DIR)/optically_thin/OpticallyThinGas_harm_leakage.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/ConstantDensitySphere_NRPyLeakage.o: $(TESTS_DIR)/constant_density_sphere/ConstantDensitySphere_NRPyLeakage.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/ConstantDensitySphere_harm_leakage.o: $(TESTS_DIR)/constant_density_sphere/ConstantDensitySphere_harm_leakage.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/harm_neutrinos.o: $(HARM_DIR)/harm_neutrinos.c
	$(CC) $(CFLAGS) -c $< -o $@

obj/harm_optical_depth.o: $(HARM_DIR)/harm_optical_depth.c
	$(CC) $(CFLAGS) -c $< -o $@

obj/NRPyEOS_benchmark.o: Leakage_tests/eos_benchmarks/NRPyEOS_benchmark.c $(INC)
	$(CC) $(CFLAGS) -c $< -o $@

obj/ZelmaniLeak_calc_taus.o: ZelmaniLeak/ZelmaniLeak_calc_taus.F90
	$(FC) $(FCFLAGS) -c $< -o $@

obj/ZelmaniLeak_calc_leak.o: ZelmaniLeak/ZelmaniLeak_calc_leak.F90
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	rm -rf $(EXE) $(OBJ) $(COBJ) obj

veryclean: clean
	rm -rf *.txt *.png *.asc
