###############################################################################
#  Makefile: OpenMP  +  CUDA
###############################################################################

TARGETS_CPU  := metrics_hll metrics_smh time_hll time_smh build_sketch selection
TARGETS_CUDA := $(addsuffix _cuda,$(TARGETS_CPU))
TARGETS      := $(TARGETS_CPU) $(TARGETS_CUDA)

SRC_DIRS := experiments/src src

# ------------------------ Compiladores y flags ----------------------------
CXX       := c++
CXXFLAGS  := -Wall -Wextra -std=c++14 -fopenmp -O3 -march=native \
             -DSEQAN_HAS_ZLIB=1 -DNDEBUG
LDFLAGS   := -lstdc++ -lm -lz -pthread

NVCC      := nvcc
CUDAFLAGS := -O3 --std=c++14 -Xcompiler "-Wall -Wextra -fopenmp -march=native" \
             -arch=sm_75 -DNDEBUG -lineinfo
LDFLAGS_CUDA := -lcudart -lm -lz -pthread

INCLUDE    := -I. -Isketch/ -Isketch/include -Isketch/include/blaze \
              -I/home/centos/dw/seqan/include -Iinclude

BUILD    := build
OBJ_DIR  := $(BUILD)/objects
BIN_DIR  := $(BUILD)

# ------------------------ Mapear ejecutables → archivo fuente -------------
metrics_hll_src       := experiments/src/metrics_hll.cpp
metrics_smh_src       := experiments/src/metrics_smh.cpp
time_hll_src          := experiments/src/time_hll.cpp
time_smh_src          := experiments/src/time_smh.cpp
build_sketch_src      := src/build_sketch.cpp
selection_src         := src/selection.cpp

metrics_hll_cuda_src  := experiments/src/metrics_hll_cuda.cu
metrics_smh_cuda_src  := experiments/src/metrics_smh_cuda.cu
time_hll_cuda_src     := experiments/src/time_hll_cuda.cu
time_smh_cuda_src     := experiments/src/time_smh_cuda.cu
build_sketch_cuda_src := src/build_sketch_cuda.cu
selection_cuda_src    := src/selection_cuda.cu

SRCS_CPU   := $(foreach t,$(TARGETS_CPU),  $($(t)_src))
SRCS_CUDA  := $(foreach t,$(TARGETS_CUDA), $($(t)_src))
SRCS       := $(SRCS_CPU) $(SRCS_CUDA)

OBJECTS_CPU  := $(SRCS_CPU:%.cpp=$(OBJ_DIR)/%.o)
OBJECTS_CUDA := $(SRCS_CUDA:%.cu=$(OBJ_DIR)/%.o)
OBJECTS      := $(OBJECTS_CPU) $(OBJECTS_CUDA)

BINARIES_CPU  := $(addprefix $(BIN_DIR)/,$(TARGETS_CPU))
BINARIES_CUDA := $(addprefix $(BIN_DIR)/,$(TARGETS_CUDA))
BINARIES      := $(BINARIES_CPU) $(BINARIES_CUDA)

# ------------------------ 7. Regla por defecto -------------------------------
all: build $(BINARIES)

# ------------------------ 8. Compilar objetos --------------------------------
$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

$(OBJ_DIR)/%.o: %.cu
	@mkdir -p $(@D)
	$(NVCC) $(CUDAFLAGS) $(INCLUDE) -c $< -o $@

# ------------------------ 9. Enlazar ejecutables -----------------------------
# CPU
$(foreach t,$(TARGETS_CPU),\
  $(eval $(BIN_DIR)/$(t): $(OBJ_DIR)/$($(t)_src:.cpp=.o))\
)

$(BIN_DIR)/%:
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ $(LDFLAGS) -o $@

# CUDA
$(foreach t,$(TARGETS_CUDA),\
  $(eval $(BIN_DIR)/$(t): $(OBJ_DIR)/$($(t)_src:.cu=.o))\
)

$(BIN_DIR)/%_cuda:
	@mkdir -p $(@D)
	$(NVCC) $(CUDAFLAGS) $(INCLUDE) $^ $(LDFLAGS_CUDA) -o $@

build:
	@mkdir -p $(OBJ_DIR)

clean:
	@echo "Removing all build directory"
	@rm -rf $(BUILD)

.PHONY: all clean build
