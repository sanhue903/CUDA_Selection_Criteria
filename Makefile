###############################################################################
#  Makefile  —  C++ (OpenMP) + CUDA
###############################################################################

############################### 1. Targets ####################################
CPU_APPS  := time_smh build_sketch selection
CUDA_APPS := selection_cuda time_smh_cuda
APPS      := $(CPU_APPS) $(CUDA_APPS)

############################ 2. Archivos fuente ###############################
SRC_CPU  := \
    experiments/src/time_smh.cpp \
    src/build_sketch.cpp          \
    src/selection.cpp

SRC_CUDA := \
    src/selection_kernels.cu      \
    src/selection_cuda.cpp        \
    experiments/src/time_smh_cuda.cpp

############################ 3. Herramientas ##################################
CXX  ?= c++
NVCC ?= nvcc

############################ 4. Flags comunes #################################
CUDA_ARCH ?= sm_86                # cámbialo si tu GPU es distinta

INC_DIRS  := -I. -Isketch -Isketch/include -Isketch/include/blaze \
             -Iseqan-library-2.4.0/include -Iinclude              \
             -I/usr/local/cuda/include

CXXFLAGS  := -O3 -std=c++17 -Wall -Wextra -march=native \
             -fopenmp -DSEQAN_HAS_ZLIB=1 -DNDEBUG $(INC_DIRS)

CUDAFLAGS := -O3 --std=c++17 -arch=$(CUDA_ARCH) \
             -Xcompiler "-Wall -Wextra -fopenmp -march=native" \
             -lineinfo -DNDEBUG $(INC_DIRS)

LDFLAGS_CPU  := -lz -pthread -fopenmp
LDFLAGS_CUDA := -lcudart -lz -Xcompiler -fopenmp

############################ 5. Carpetas build ################################
BUILD_DIR := build
OBJ_DIR   := $(BUILD_DIR)/objects
BIN_DIR   := $(BUILD_DIR)

########################### 6. Objetos derivados ##############################
CPP_OBJS := $(addprefix $(OBJ_DIR)/,$(SRC_CPU:.cpp=.o)) \
            $(addprefix $(OBJ_DIR)/,$(filter %.cpp,$(SRC_CUDA:.cpp=.o)))
CU_OBJS  := $(addprefix $(OBJ_DIR)/,$(SRC_CUDA:.cu=.o))
ALL_OBJS := $(CPP_OBJS) $(CU_OBJS)

########################## 7. Reglas de alto nivel ############################
.PHONY: all clean
all: $(addprefix $(BIN_DIR)/,$(APPS))

clean:
	@echo "Cleaning build directory…"
	@rm -rf $(BUILD_DIR)

########################### 8. Reglas patrón ##################################
# Compilación .cpp → .o
$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compilación .cu → .o
$(OBJ_DIR)/%.o: %.cu
	@mkdir -p $(@D)
	$(NVCC) $(CUDAFLAGS) -c $< -o $@

########################### 9. Enlazado final #################################
# --- CPU ---
$(BIN_DIR)/time_smh: $(addprefix $(OBJ_DIR)/,experiments/src/time_smh.o)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS_CPU) -o $@

$(BIN_DIR)/build_sketch: $(addprefix $(OBJ_DIR)/,src/build_sketch.o)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS_CPU) -o $@

$(BIN_DIR)/selection: $(addprefix $(OBJ_DIR)/,src/selection.o)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS_CPU) -o $@

# --- GPU ---
# selection_cuda  ← host (selection_cuda.cpp) + kernels (.cu)
$(BIN_DIR)/selection_cuda: \
	$(addprefix $(OBJ_DIR)/,src/selection_cuda.o src/selection_kernels.o)
	@mkdir -p $(@D)
	$(NVCC) $(CUDAFLAGS) $^ $(LDFLAGS_CUDA) -o $@

# time_smh_cuda   ← host (time_smh_cuda.cpp) + kernels (.cu)
$(BIN_DIR)/time_smh_cuda: \
	$(addprefix $(OBJ_DIR)/,experiments/src/time_smh_cuda.o src/selection_kernels.o)
	@mkdir -p $(@D)
	$(NVCC) $(CUDAFLAGS) $^ $(LDFLAGS_CUDA) -o $@
