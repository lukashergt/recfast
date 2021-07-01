# Makefile for recfast
# Lukas Hergt, 01.06.2021
#
# Usage:
# make        # compile all binaries
# make test   # execute a basic run test
# make clean  # remove ALL binaries and objects


# Whether to compile in debugging mode (default: false)
DEBUG=1


###############################################################################
############################# Prepare directories #############################
###############################################################################
MAKE_DIR = $(PWD)
BUILD_DIR = $(MAKE_DIR)/build
TEST_DIR = $(MAKE_DIR)/test/example_data
.base:
	if ! [ -e $(BUILD_DIR) ]; then mkdir $(BUILD_DIR) ; fi;
	touch build/.base
vpath %.o build
vpath .base build


###############################################################################
############################# Running with intel ##############################
###############################################################################
ifeq "$(shell which ifort >/dev/null 2>&1; echo $$?)" "0" 
FC = ifort

# default flags
# -------------
# fpp                    : perform preprocessing
# fpic                   : shared object libraries
FFLAGS += -fpp -fpic -heap-arrays
MODFLAG = -module $(BUILD_DIR)

ifeq ($(DEBUG),1)
# debugging mode
# --------------
# g              : enable gnu debugger compatibility
# O0             : no optimisation
# traceback      : create a backtrace on system failure
# check all      : all checks (whilst compiling)
# warn all       : all warnings (whilst running)
# ftrapuv        : Traps uninitialized variables by setting them to very large values
# debug all      : additional debugging information
# gen-interfaces : generate an interface block for each routine
# warn-interfaces: warn on these interface blocks
FFLAGS += -g -O0 -traceback -check all,noarg_temp_created -warn all -ftrapuv -debug all -gen-interfaces -warn-interfaces
else
# optimised mode
# --------------
#   ipo          : interprocedural optimization (optimize entire program)
#   O3           : maximum optimisation
#   no-prec-div  : slightly less precise floating point divides, but speeds up
#   static       : link intel libraries statically
#   xHost        : maximise architecture usage
#   w            : turn off all warnings
#   vec-report0  : disables printing of vectorizer diagnostic information
#   opt-report0  : disables printing of optimization reports
IPO = -ipo
FFLAGS += $(IPO) -O3 -no-prec-div $(HOST) -w -vec-report0 -qopt-report0
endif


###############################################################################
############################# Running with gnu ################################
###############################################################################
else ifeq "$(shell which gfortran >/dev/null 2>&1; echo $$?)" "0"
FC = gfortran

# default flags
# --------------
# free-line-length-none : turn of line length limitation (why is this not a default??)
# cpp                   : perform preprocessing
# fPIC                  : for compiling a shared object library
FFLAGS += -ffree-line-length-none -cpp -fPIC -fno-stack-arrays 
MODFLAG = -J $(BUILD_DIR)

ifeq ($(DEBUG),1)
# debugging mode
# --------------
# g             : enable gnu debugger compatibility
# O0            : no optimisation
# Wall          : all warnings
# Wextra        : even more warnings
# pedantic      : check for language features not part of f95 standard
# implicit-none : specify no implicit typing 
# backtrace     : produce backtrace of error
# fpe-trap      : search for floating point exceptions (dividing by zero etc)
FFLAGS += -g -O0 -Wall -Wextra -pedantic -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=zero,overflow 
else
# optimised mode
# --------------
# Ofast : maximum optimisation
FFLAGS += -Ofast
endif


endif


###############################################################################
################################ Make targets #################################
###############################################################################

SRCS := $(wildcard *.f03)
BINS := $(SRCS:%.f03=%)

all: $(BINS) test

# Build step for executable
%: $(BUILD_DIR)/%.o
	$(FC) $(FFLAGS) $(MODFLAG) $< -o $@

# Build step for fortran source
$(BUILD_DIR)/%.o: %.f03 .base
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

# Run a basic test with input from example.ini
test: $(BINS)
	@echo
	@echo "Test with test/example_data/example.ini"
	@echo "======================================="
	./recfast < $(TEST_DIR)/example.ini
	@echo

clean:
	@echo "Cleaning recfast"
	@echo "================"
	rm -rvf $(BUILD_DIR)
	rm -vf $(MAKE_DIR)/$(BINS)
	rm -vf $(MAKE_DIR)/test.out

