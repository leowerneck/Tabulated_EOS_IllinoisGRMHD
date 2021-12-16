BASEDIR = .
WVU_EOS_DIR = $(BASEDIR)/WVU_EOS
LEAKAGE_DIR = $(BASEDIR)/NRPyLeakage

HDF5_DIR = /usr/local/opt/hdf5
#HDF5_DIR = /usr/lib/x86_64-linux-gnu/hdf5/serial

CC       = gcc
CFLAGS   = -Wall -O2 -fopenmp -march=native -I$(BASEDIR) -I$(WVU_EOS_DIR) -I$(LEAKAGE_DIR) -I$(HDF5_DIR)/include

CXX      = g++
CXXFLAGS = -Wall -O2 -fopenmp -march=native -I$(BASEDIR) -I$(WVU_EOS_DIR) -I$(LEAKAGE_DIR) -I$(HDF5_DIR)/include

LD_FLAGS = -L$(HDF5_DIR)/lib -lhdf5

SRC  = standalone.cc
EXE  = $(SRC:.cc=)

WVU_EOS_SRC = WVU_EOS_Tabulated_general_interpolators.cc \
              WVU_EOS_Tabulated_known_T.cc \
              WVU_EOS_Tabulated_unknown_T.cc \
              WVU_EOS_readtable.cc \
              WVU_EOS_free_memory.cc

WVU_EOS_INC = $(WVU_EOS_DIR)/WVU_EOS_Tabulated_headers.hh \
              $(WVU_EOS_DIR)/WVU_EOS_Tabulated_helpers.hh

LEAKAGE_SRC = Leakage_Fermi_Dirac_integrals.cc Leakage_compute_free_rates.cc

LEAKAGE_INC = $(LEAKAGE_DIR)/Leakage.hh

OBJ = $(addprefix obj/,$(SRC:.cc=.o)) \
      $(addprefix obj/,$(WVU_EOS_SRC:.cc=.o)) \
      $(addprefix obj/,$(LEAKAGE_SRC:.cc=.o))

INC = Basic_defines.hh $(WVU_EOS_INC) $(LEAKAGE_INC)

all: obj_dir $(EXE)

obj_dir:
	mkdir -p obj

$(EXE): $(OBJ) $(COBJ)
	$(CXX) $(CXXLAGS) $(OBJ) $(COBJ) -o $(EXE) $(LD_FLAGS)

obj/standalone.o: standalone.cc $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj/WVU_EOS_Tabulated_general_interpolators.o: $(WVU_EOS_DIR)/WVU_EOS_Tabulated_general_interpolators.cc $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj/WVU_EOS_Tabulated_known_T.o: $(WVU_EOS_DIR)/WVU_EOS_Tabulated_known_T.cc $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj/WVU_EOS_Tabulated_unknown_T.o: $(WVU_EOS_DIR)/WVU_EOS_Tabulated_unknown_T.cc $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj/WVU_EOS_readtable.o: $(WVU_EOS_DIR)/WVU_EOS_readtable.cc $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj/WVU_EOS_free_memory.o: $(WVU_EOS_DIR)/WVU_EOS_free_memory.cc $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj/Leakage_Fermi_Dirac_integrals.o: $(LEAKAGE_DIR)/Leakage_Fermi_Dirac_integrals.cc $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj/Leakage_compute_free_rates.o: $(LEAKAGE_DIR)/Leakage_compute_free_rates.cc $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(EXE) $(OBJ) $(COBJ) obj
