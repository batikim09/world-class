# ===========================================
#
#  Lazy Makefile
#
#  Copyright 2018, Yukara Ikemiya
#
# ===========================================

PROGRAM_NAME := test
SHAREDLIB_NAME := libworld.so

# compiler
CXX := g++ -std=c++11
CXXFLAGS := -Wall -Wextra -O3 -mavx -fopenmp -pipe
CXX_SHARED_FLAGS := -shared -fPIC -Wall

# external libraries
LIBS := -lfftw3 -lsndfile
EXT_INCLUDES := -I/usr/local/include/eigen3/Eigen

# source directories
ROOT_DIR := .
SOURCE_DIR := $(ROOT_DIR)
HEADER_DIR := $(ROOT_DIR)

# get all source files
SOURCES := $(shell find . -name "*.cpp" -or -name "*.c")
HEADERS := $(shell find . -name "*.hpp" -or -name "*.h")

# get all sub-directories (except a .git directory)
INCLUDES := $(addprefix -I,$(shell find $(ROOT_DIR) -path $(ROOT_DIR)/.git -prune -o -type d -print)) $(EXT_INCLUDES)

# output directories
OUT_DIR := build
PROGRAM_DIR := $(OUT_DIR)/bin
SHARED_DIR := $(OUT_DIR)/lib
OBJ_DIR := $(OUT_DIR)/obj
DEPEND_DIR := $(OUT_DIR)/depend

# output files
PROGRAM := $(PROGRAM_DIR)/$(PROGRAM_NAME)
SHARED := $(SHARED_DIR)/$(SHAREDLIB_NAME)
SOURCE_NAMES = $(SOURCES) # $(notdir $(SOURCES))
OBJS := $(addprefix $(OBJ_DIR)/,$(SOURCE_NAMES:.cpp=.o))
DEPENDS := $(addprefix $(DEPEND_DIR)/,$(SOURCE_NAMES:.cpp=.depend))

# rules
.PHONY: all
all: $(DEPENDS) $(PROGRAM)
$(PROGRAM): $(SHARED)
	@echo "\ngenerating $@"
	@mkdir -p $(PROGRAM_DIR)
	cp $(SHARED) $(PROGRAM_DIR)/$(SHAREDLIB_NAME)
	cp $(SHARED) $(ROOT_DIR)/$(SHAREDLIB_NAME) #this is just to play a demo. You should install the shared lib into system.
	$(CXX) $(CXXFLAGS) -I$(INCLUDES) $(SOURCE_DIR)/test/test.cpp -o $(PROGRAM) -lworld -L$(SHARED_DIR)

.PHONY: lib
lib: $(DEPENDS) $(SHARED)
$(SHARED): $(OBJS)
	@echo "\ngenerating $@"
	@mkdir -p $(SHARED_DIR)
	$(CXX) $(CXX_SHARED_FLAGS) -o $(SHARED) $(LIBS) $(OBJS) -lc

$(DEPEND_DIR)/%.depend: %.cpp $(HEADERS)
	@echo "\ngenerating $@"
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	@$(CXX) $(CXXFLAGS) $(INCLUDES) -MM $< > $@ $(LIBS)

$(OBJ_DIR)/%.o: %.cpp
	@echo "\ngenerating $@"
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CXX) $(CXXFLAGS) -fPIC $(INCLUDES) -c $^ -o $@ $(LIBS)

ifneq "$(MAKECMDGOALS)" "clean"
-include $(DEPENDS)
endif

.PHONY : clean
clean:
	$(RM) -r $(OUT_DIR)
