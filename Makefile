# Makefile for recfast
# Lukas Hergt, 01.06.2021
#
# Usage:
# make        # compile all binaries and install the python wrapper
# make test   # execute a basic run test
# make clean  # remove ALL objects
# make purge  # remove ALL objects and binaries


# Whether to compile in debugging mode (default: false)
DEBUG=1


###############################################################################
############################# Prepare directories #############################
###############################################################################
MAKE_DIR = $(PWD)
SOURCE_DIR = $(MAKE_DIR)/src
BUILD_DIR = $(MAKE_DIR)/build
TEST_DIR = $(MAKE_DIR)/test/example_data
.base:
	if ! [ -e $(BUILD_DIR) ]; then mkdir $(BUILD_DIR) ; fi;
	touch build/.base
vpath %.o build
vpath .base build


###############################################################################
############################## Compilation flags ##############################
###############################################################################
FFLAGS = -fPIC
#CFLAGS = -fPIC -I/usr/include/python3.9 -I/usr/lib/python3.9/site-packages/numpy/core/include/
#CFLAGS = -fPIC $(shell pkg-config --dont-define-prefix --cflags python3) $(shell python -c "import numpy; print('-I' + numpy.get_include())")
CFLAGS = -fPIC $(shell python3-config --includes) $(shell python -c "import numpy; print('-I' + numpy.get_include())")


###############################################################################
############################# Running with intel ##############################
###############################################################################
ifeq "$(shell which ifort >/dev/null 2>&1; echo $$?)" "0" 
FC = ifort
CC = icc

# default flags
# -------------
# fpp                    : perform preprocessing
# fpic                   : shared object libraries
FFLAGS += -fpp -heap-arrays
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
CFLAGS += -Wall
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
CC = gcc

# default flags
# --------------
# free-line-length-none : turn of line length limitation (why is this not a default??)
# cpp                   : perform preprocessing
# fPIC                  : for compiling a shared object library
FFLAGS += -ffree-line-length-none -cpp -fno-stack-arrays 
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
CFLAGS += -Wall
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

SRCS := precision constants ode_solver recfast
OBJS := $(patsubst %,$(BUILD_DIR)/%.o,$(SRCS))

all: .base $(OBJS) recfast pyrecfast.so test

recfast: $(OBJS)
	$(FC) $(FFLAGS) $(MODFLAG) $^ -o $@

pyrecfast: .base pyrecfast.so

pyrecfast.so: $(OBJS) $(BUILD_DIR)/python_wrapper.o $(BUILD_DIR)/pyrecfast.o
	$(FC) $(FFLAGS) $(MODFLAG) -shared $^ -o $@

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.f08 .base
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.c .base
	$(CC) $(CFLAGS) -c $< -o $@

$(SOURCE_DIR)/%.c: $(SOURCE_DIR)/%.pyx $(BUILD_DIR)/python_wrapper.o
	python setup.py sdist bdist_wheel egg_info --egg-base $(BUILD_DIR)


# Run a basic test with input from example.ini
test: recfast
	@echo
	@echo
	@echo
	@echo "Test with test/example_data/example.ini"
	@echo "======================================="
	./recfast < $(TEST_DIR)/example.ini
	@echo
	tail -n 5 $(TEST_DIR)/example_new_CODATA_AME_2photon_ionisationlevel.out
	tail -n 5 $(MAKE_DIR)/test.out
	@echo

clean:
	@echo "Cleaning recfast build"
	@echo "======================"
	rm -rvf $(BUILD_DIR)
	rm -vf $(SOURCE_DIR)/pyrecfast.c
	@echo

purge: clean
	@echo "Purging recfast folder from executables and any remaining .mod, .o, .c etc."
	@echo "==========================================================================="
	rm -vf $(MAKE_DIR)/recfast
	rm -vf $(MAKE_DIR)/*.mod
	rm -vf $(MAKE_DIR)/*.o
	rm -vf $(MAKE_DIR)/*.c
	rm -vf $(MAKE_DIR)/*.so
	rm -rvf $(MAKE_DIR)/dist
	rm -rvf $(MAKE_DIR)/pyrecfast.egg-info
	rm -vf $(MAKE_DIR)/test.out
	pip uninstall pyrecfast
	@echo

