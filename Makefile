cpu_executables := sd1 os1
cpu_executables += t_instrument t_interstellar t_nsx
cpu_executables += RingEq
cpu_headers := $(wildcard cpu/*)

ifneq (command line,$(origin CXX))
  CXX := clang++
endif

CXXFLAGS := -std=c++20 -O3 -march=native -ffast-math

$(cpu_executables): % : cpu/%.cpp $(cpu_headers) build
	$(CXX) $(CXXFLAGS) $< -o build/$@

build:
	mkdir -p build
