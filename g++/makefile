# ONON
# Set main options (can be at command line)
ifndef $(SYSTEM) 
SYSTEM=mac1
endif

# Possible modes are debug and run
ifndef $(MODE)
MODE=debug
endif

# Set main compile option
CXX = g++

# Set local paths
SRCLIBCPP = ../src/libcpps
SRCLIBHEADER = ../src/libheaders
SRCLIBMAIN = ../src/projects/ebola
NRCHEADDIR = ../nrcpptmp/other
NRCSRCDIR = ../nrcpptmp/recipes

# Set header paths and library paths for gsl if required
# Note windows options are already on the default compiler paths
# Linux likely to be the same as windows
ifeq ($(SYSTEM),mac)
GSLHEADER = /opt/local/include
GSLLIB = /opt/local/lib
else ifeq ($(SYSTEM),windows)
GSLHEADER = /usr/include
GSLLIB = /lib
else ifeq ($(SYSTEM),mac1)
GSLHEADER = /usr/local/include
GSLLIB = /usr/local/lib
endif

INCLUDEDIRS = -I$(NRCHEADDIR) -I$(SRCLIBHEADER) -I$(GSLHEADER)
LDFLAGS =   -L$(GSLLIB) -lgsl -lgslcblas

# System variable to tell the compiler where to look for .cpps and .hs
VPATH = $(SRCLIBCPP) $(SRCLIBHEADER) $(SRCLIBMAIN) $(NRCHEADDIR) $(NRCSRCDIR)

# Set standatd compile variables
ifeq ($(MODE),run)
CXXFLAGS = -O3 $(INCLUDEDIRS) 
else ifeq ($(MODE),debug)
 
# Possible flags -Wall for all warnings
CXXFLAGS = -Wall -ggdb $(INCLUDEDIRS) 
endif

# Define the object files
OBJFILES = 	ran1.o ran2.o gamdev.o gammln.o poidev.o \
			SR_BytePackedInt.o SR_CachedLookups.o SR_DensityField.o SR_EventMatrix.o \
			SR_GatherRunInformation.o SR_GridHex.o SR_InitialConditions.o SR_Kernels.o \
			SR_Label.o SR_PagesForThings.o SR_Parameter.o SR_Stats.o SR_StripComments.o \
			SR_Utility.o SR_Workplaces.o ebola_common.o

all : ebola_build.exe ebola_run.exe

ebola_run.exe : $(OBJFILES) ebola_run.o  
	$(CXX) $(OBJFILES) ebola_run.o -o $@ $(LDFLAGS) 

ebola_build.exe : $(OBJFILES) ebola_build.o  
	$(CXX) $(OBJFILES) ebola_build.o -o $@ $(LDFLAGS)

clean:
	rm *.o *.exe
