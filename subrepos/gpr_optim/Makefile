# Just a helper
RM := rm -rf

# List of libraries and additional header files
LIBS_PATHS := -L./googletest/googletest/build/lib
INCS_PATHS := -I./googletest/googletest/include
LIBS := -lgtest -lgtest_main -lpthread
INCL :=

# Compiler options
CXXFLAGS := -O3 -std=c++11 -Wall -DNDEBUG

# Additional options
EXTRA_OPTS :=

# Name of the compiled executable
EXE_NAME := gpr.out

# List of source files
SRCS := \
	tests/data_types/FieldTest.cpp \
	tests/gpr/auxiliary/DistanceTest.cpp \
	tests/gpr/auxiliary/ProblemSetUpTest.cpp \
	tests/gpr/auxiliary/GradientTest.cpp \
	tests/gpr/covariance_functions/SexpatCFTest.cpp \
	tests/gpr/dimer/DimerTest.cpp \
	tests/gpr/observation_models/LikGaussianTest.cpp \
	tests/gpr/ml/GaussianProcessRegressionTest.cpp \
	tests/gpr/ml/SCGTest.cpp \
	tests/gpr/prior/PriorGaussianTest.cpp \
	tests/gpr/prior/PriorLogUnifTest.cpp \
	tests/gpr/prior/PriorSqrttTest.cpp \
	tests/gpr/prior/PriorTTest.cpp \
	tests/gpr/potentials/EAMPotentialTest.cpp \
	tests/gpr/AtomicDimerTest.cpp \
	tests/managers/io/FileManagerTest.cpp \
	backend/DistributionFunctions.cpp \
	gpr/optimization/LBFGS.cpp \
	gpr/auxiliary/Distance.cpp \
	gpr/auxiliary/ProblemSetUp.cpp \
	gpr/auxiliary/Gradient.cpp \
	gpr/covariance_functions/ConstantCF.cpp \
	gpr/covariance_functions/SexpatCF.cpp \
	gpr/dimer/Dimer.cpp \
	gpr/ml/GaussianProcessRegression.cpp \
	gpr/ml/SCG.cpp \
	gpr/observation_models/LikGaussian.cpp \
	gpr/prior/PriorGaussian.cpp \
	gpr/prior/PriorLogUnif.cpp \
	gpr/prior/PriorSqrtt.cpp \
	gpr/prior/PriorT.cpp \
	gpr/prior/PriorUnif.cpp \
	gpr/potentials/EAMPotential.cpp \
	gpr/AtomicDimer.cpp \
	gpr/AtomicDimerInit.cpp \
	managers/io/FileManager.cpp \
	main.cpp

# Filter of source files
FILTER :=

# Directory for object files
OBJDIR = ./obj

# Object files
OBJECTS = $(patsubst %.cpp, $(OBJDIR)/%.o,$(SRCS))

# All Target
all: $(OBJECTS) make_exec

$(OBJECTS): | obj

# Create all necessary directories for the object files
obj:
	@mkdir -p $(OBJDIR)
	@mkdir -p $(OBJDIR)/tests/backend
	@mkdir -p $(OBJDIR)/tests/data_types
	@mkdir -p $(OBJDIR)/tests/gpr
	@mkdir -p $(OBJDIR)/tests/gpr/auxiliary
	@mkdir -p $(OBJDIR)/tests/gpr/covariance_functions
	@mkdir -p $(OBJDIR)/tests/gpr/dimer
	@mkdir -p $(OBJDIR)/tests/gpr/ml
	@mkdir -p $(OBJDIR)/tests/gpr/observation_models
	@mkdir -p $(OBJDIR)/tests/gpr/prior
	@mkdir -p $(OBJDIR)/tests/gpr/potentials
	@mkdir -p $(OBJDIR)/tests/managers
	@mkdir -p $(OBJDIR)/tests/managers/io
	@mkdir -p $(OBJDIR)/backend
	@mkdir -p $(OBJDIR)/gpr
	@mkdir -p $(OBJDIR)/gpr/optimization
	@mkdir -p $(OBJDIR)/gpr/auxiliary
	@mkdir -p $(OBJDIR)/gpr/covariance_functions
	@mkdir -p $(OBJDIR)/gpr/dimer
	@mkdir -p $(OBJDIR)/gpr/ml
	@mkdir -p $(OBJDIR)/gpr/observation_models
	@mkdir -p $(OBJDIR)/gpr/prior
	@mkdir -p $(OBJDIR)/gpr/potentials
	@mkdir -p $(OBJDIR)/managers
	@mkdir -p $(OBJDIR)/managers/io
	@echo ' '
	@echo 'Building options:'
	@echo -e '   Compiler:        $(CXX)'
	@echo -e '   Flags:           $(CXXFLAGS)'
	@echo -e '   Extra options:   $(EXTRA_OPTS)'
	@echo -e '   Executable name: $(EXE_NAME)'
	@echo ' '

# Build object files and print out their names
$(OBJDIR)/%.o: %.cpp
	@echo 'Building file: $<'
	@$(CXX) $(filter $@.cpp,$(FILTER)) $(CXXFLAGS) $(INCS_PATHS) $(INCL) $(EXTRA_OPTS) -c $< -o $@

# Create lib file
make_exec: $(OBJECTS)
	$(CXX) -o $(EXE_NAME) $(OBJECTS) $(LIBS_PATHS) $(LIBS)

# Other targets
clean:
	-$(RM) $(OBJDIR) $(EXE_NAME)

.PHONY: all clean
.SECONDARY:
