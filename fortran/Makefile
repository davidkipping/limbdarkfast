FC = gfortran
FFLAGS ?= -std=f2008 -O2 -Wall -Wextra -pedantic -fimplicit-none

all: build/test_limbdarkfast build/example

build/limbdarkfast.o: limbdarkfast.f90
	@mkdir -p build
	$(FC) $(FFLAGS) -Jbuild -c $< -o $@

build/%: %.f90 build/limbdarkfast.o
	$(FC) $(FFLAGS) -Ibuild -o $@ $< build/limbdarkfast.o

test: build/test_limbdarkfast
	./build/test_limbdarkfast ../tests/vectors

clean:
	rm -rf build

.PHONY: all test clean
