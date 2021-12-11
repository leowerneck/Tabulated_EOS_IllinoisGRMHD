BASEDIR = $(shell pwd)
WVU_EOS_DIR = $(BASEDIR)/WVU_EOS

CC       = gcc
CFLAGS   = -Wall -O2 -fopenmp -march=native -I$(WVU_EOS_DIR) -I/usr/local/opt/hdf5/include

CXX      = g++
CXXFLAGS = -Wall -O2 -fopenmp -march=native -I$(WVU_EOS_DIR) -I/usr/local/opt/hdf5/include

LD_FLAGS = -L/usr/local/opt/hdf5/lib -lhdf5

EXE  = standalone
CSRC =

WVU_EOS_SRC = WVU_EOS_Tabulated_general_interpolators.cc \
              WVU_EOS_Tabulated_known_T.cc \
              WVU_EOS_Tabulated_unknown_T.cc \
              WVU_EOS_readtable.cc \
              WVU_EOS_free_memory.cc

WVU_EOS_INC = $(WVU_EOS_DIR)/WVU_EOS_Tabulated_headers.hh \
              $(WVU_EOS_DIR)/WVU_EOS_Tabulated_helpers.hh

SRC = standalone.cc


OBJ = $(addprefix obj/,$(SRC:.cc=.o)) \
      $(addprefix obj/,$(WVU_EOS_SRC:.cc=.o))

INC = Basic_defines.hh $(WVU_EOS_INC)

all: obj_dir $(EXE)

obj_dir:
	mkdir -p obj

$(EXE): $(OBJ) $(COBJ)
	$(CXX) $(CXXLAGS) $(OBJ) $(COBJ) -o $(EXE) $(LD_FLAGS)

obj/standalone.o: standalone.cc $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj/WVU_EOS_Tabulated_general_interpolators.o: WVU_EOS/WVU_EOS_Tabulated_general_interpolators.cc $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj/WVU_EOS_Tabulated_known_T.o: WVU_EOS/WVU_EOS_Tabulated_known_T.cc $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj/WVU_EOS_Tabulated_unknown_T.o: WVU_EOS/WVU_EOS_Tabulated_unknown_T.cc $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj/WVU_EOS_readtable.o: WVU_EOS/WVU_EOS_readtable.cc $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj/WVU_EOS_free_memory.o: WVU_EOS/WVU_EOS_free_memory.cc $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(EXE) $(OBJ) $(COBJ) obj
