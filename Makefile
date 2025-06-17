# Executables to build
TARGETS := metrics_hll metrics_smh time_hll time_smh build_sketch selection

# Source directories
SRC_DIRS := experiments/src src

# Compiler and flags
CXX      := c++
CXXFLAGS := -Wall -Wextra -std=c++14 -fopenmp -O3 -march=native -DSEQAN_HAS_ZLIB=1 -DNDEBUG
LDFLAGS  := -lstdc++ -lm -lz -pthread
INCLUDE  := -I. -Isketch/ -Isketch/include -Isketch/include/blaze
INCLUDE  += -I/home/centos/dw/seqan/include

# Directories
BUILD    := build
OBJ_DIR  := $(BUILD)/objects
BIN_DIR  := $(BUILD)

# Map target names to their .cpp files
metrics_hll_src            := experiments/src/metrics_hll.cpp
metrics_smh_src            := experiments/src/metrics_smh.cpp
time_hll_src               := experiments/src/time_hll.cpp
time_smh_src               := experiments/src/time_smh.cpp
build_sketch_src           := src/build_sketch.cpp
selection_src              := src/selection.cpp

# Derive all source files and object files
SRCS    := $(foreach t, $(TARGETS), $($(t)_src))
OBJECTS := $(SRCS:%.cpp=$(OBJ_DIR)/%.o)
BINARIES := $(addprefix $(BIN_DIR)/, $(TARGETS))

# Default target: build all
all: build $(BINARIES)

# Compile each source file to object file
$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

# Link each object to its binary
# One rule for each binary using eval
$(foreach t, $(TARGETS), \
  $(eval $(BIN_DIR)/$(t): $(OBJ_DIR)/$($(t)_src:.cpp=.o)) \
)

$(BIN_DIR)/%:
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ $(LDFLAGS) -o $@

# Create build directories
build:
	@mkdir -p $(OBJ_DIR)/experiments/src
	@mkdir -p $(OBJ_DIR)/src

# Clean rule: delete build directory 
clean:
	@echo "Removing all build directory"
	@rm -rf $(BUILD)

.PHONY: all clean build

