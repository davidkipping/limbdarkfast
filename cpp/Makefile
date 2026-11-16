CXX ?= c++
CXXFLAGS ?= -std=c++11 -O2 -Wall -Wextra -pedantic

all: build/test_limbdarkfast build/example

build/%: %.cpp limbdarkfast.hpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -o $@ $<

test: build/test_limbdarkfast
	./build/test_limbdarkfast ../tests/vectors

clean:
	rm -rf build

.PHONY: all test clean
