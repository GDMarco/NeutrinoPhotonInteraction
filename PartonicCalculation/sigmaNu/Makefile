# Only the Cuba library is required for the numerical integration
# GSL path
GSL_PATH = $(shell gsl-config --prefix)
# LHAPDF path
LHAPDF_PATH = $(shell lhapdf-config --prefix)
# CUBA path
CUBA_PATH = $(shell cuba-config --prefix)
# your CUBA 4.0 path (testing: /usr/local/Cuba-4.0)
# CUBA_PATH = /Users/rhorry/Research/BUILDS/CUBA/

# compiler flags: -Wall -O2
CXXFLAGS = -Wall -O2 -std=c++11

# compiler name: g++ is default
FF = gfortran
#
INCLUDE_PATH =
LIBS_PATH = 

ifneq "$(GSL_PATH)" ""
INCLUDE_PATH += $(shell gsl-config --cflags)
LIBS_PATH += $(shell gsl-config --libs)
endif

ifneq "$(LHAPDF_PATH)" "" 
INCLUDE_PATH +=  -I$(shell lhapdf-config --incdir)
LIBS_PATH += $(shell lhapdf-config --libs)
endif

# Note cuba-config not released with installation, i created my own
ifneq "$(CUBA_PATH)" ""
INCLUDE_PATH += -I$(CUBA_PATH)/include
LIBS_PATH += -L$(CUBA_PATH)/lib -lcuba
endif

# Define code compilation structure
SRCDIR = ./src
BINDIR = ./bin
OBJDIR = ./obj
INCDIR = ./inc

# Define various source files (i.e. the core files)
SOURCES = integration.cc cuts.cc tools.cc var.cc ME2_Analytic.cc dsigma.cc 
# Append to this those for the various processes
# SIGMA_SOURCES = test_dsigma.cc

# Combine these
SOURCES += $(SIGMA_SOURCES)

# Object list obtained from the source list
_OBJ = $(SOURCES:.cc=.o)
OBJ = $(patsubst %,$(OBJDIR)/%,$(_OBJ))

# Add local source/header files 
INCLUDE_PATH += -I$(SRCDIR) -I$(INCDIR)

CXXFLAGS += $(INCLUDE_PATH)
CXXFLAGS += "-Wl,-no_warn_duplicate_libraries"

# Build all binaries
all: $(BINDIR)/Sigma.exe

# Program for QQbar production in e+e- collisions
$(BINDIR)/Sigma.exe: $(SRCDIR)/test_dsigma.cc $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@ $< $(LIBS_PATH)

$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	$(CXX) -c -o $@ $< $(CXXFLAGS)

clean:
	rm -f $(BINDIR)/*.exe
	rm -vf $(OBJDIR)/*.o

.PHONY: clean
