###############################################################################
#  Makefile: OpenMP  +  CUDA
###############################################################################

TARGETS_CPU  := time_smh build_sketch selection
TARGETS_CUDA := selection_cuda
TARGETS      := $(TARGETS_CPU) $(TARGETS_CUDA)

SRC_DIRS := experiments/src src

# ------------------------ Compiladores y flags ----------------------------
CXX       := c++
CXXFLAGS  := -Wall -Wextra -std=c++14 -fopenmp -O3 -march=native \
             -DSEQAN_HAS_ZLIB=1 -DNDEBUG
LDFLAGS   := -lstdc++ -lm -lz -pthread

NVCC      := nvcc
CUDAFLAGS := -O3 --std=c++14 -Xcompiler "-Wall -Wextra -fopenmp -march=native" \
             -arch=sm_86 -DNDEBUG -lineinfo
LDFLAGS_CUDA := -lcudart -lm -lz

INCLUDE := -I. -Isketch/ -Isketch/include -Isketch/include/blaze \
           -Iseqan-library-2.4.0/include -Iinclude

BUILD    := build
OBJ_DIR  := $(BUILD)/objects
BIN_DIR  := $(BUILD)

###############################################################################
# Fuentes individuales por ejecutable
###############################################################################
# CPU
time_smh_src      := experiments/src/time_smh.cpp
build_sketch_src  := src/build_sketch.cpp
selection_src     := src/selection.cpp

# CUDA
# (Si luego agregas build_sketch_cuda o time_smh_cuda, solo los agregas aquí)
selection_cuda_src := src/selection_cuda.cu
selection_main_src := src/selection_main.cpp

###############################################################################
# Listas de fuentes
###############################################################################
SRCS_CPU  := $(time_smh_src) $(build_sketch_src) $(selection_src)
SRCS_CUDA := $(selection_cuda_src) $(selection_main_src)
SRCS      := $(SRCS_CPU) $(SRCS_CUDA)

###############################################################################
# Listas de objetos
###############################################################################
OBJECTS_CPU  := $(SRCS_CPU:%.cpp=$(OBJ_DIR)/%.o)
OBJECTS_CUDA := $(selection_main_src:src/%.cpp=$(OBJ_DIR)/src/%.o) \
                $(selection_cuda_src:src/%.cu=$(OBJ_DIR)/src/%.o)
OBJECTS      := $(OBJECTS_CPU) $(OBJECTS_CUDA)

###############################################################################
# Binarios a generar
###############################################################################
BINARIES_CPU  := $(addprefix $(BIN_DIR)/,$(TARGETS_CPU))
BINARIES_CUDA := $(BIN_DIR)/selection_cuda
BINARIES      := $(BINARIES_CPU) $(BINARIES_CUDA)

###############################################################################
# 7. Regla por defecto
###############################################################################
all: build $(BINARIES)

###############################################################################
# 8. Compilar objetos
###############################################################################
$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

$(OBJ_DIR)/%.o: %.cu
	@mkdir -p $(@D)
	$(NVCC) $(CUDAFLAGS) $(INCLUDE) -c $< -o $@

###############################################################################
# 9. Enlazar ejecutables
###############################################################################
# CPU targets (uno a uno, limpio y claro)
$(BIN_DIR)/time_smh: $(OBJ_DIR)/experiments/src/time_smh.o
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ $(LDFLAGS) -o $@

$(BIN_DIR)/build_sketch: $(OBJ_DIR)/src/build_sketch.o
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ $(LDFLAGS) -o $@

$(BIN_DIR)/selection: $(OBJ_DIR)/src/selection.o
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ $(LDFLAGS) -o $@

# CUDA target (host+device juntos)
$(BIN_DIR)/selection_cuda: $(OBJ_DIR)/src/selection_main.o $(OBJ_DIR)/src/selection_cuda.o
	@mkdir -p $(@D)
	$(NVCC) $(CUDAFLAGS) $(INCLUDE) $^ $(LDFLAGS_CUDA) -o $@

build:
	@mkdir -p $(OBJ_DIR)

clean:
	@echo "Removing all build directory"
	@rm -rf $(BUILD)

.PHONY: all clean build
