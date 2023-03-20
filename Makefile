bin = ./bin/
src = ./src/
obj = ./obj/

CXX = g++
CXXFLAGS = -O3 -g -c

HDRS = $(src)ShallowWater.h
OBJS = $(obj)main.o $(obj)TimeIntegrate.o $(obj)derivatives.o $(obj)SetInitialConditions.o $(obj)output.o $(obj)SetParameters.o
LIBS = -lblas -lboost_program_options -fopenmp

default: $(bin)main.out

$(bin)main.out: $(OBJS)
	$(CXX) -o $@ $^ $(LIBS)

$(obj)main.o: $(src)main.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(obj)TimeIntegrate.o: $(src)TimeIntegrate.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(obj)derivatives.o: $(src)derivatives.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(obj)SetInitialConditions.o: $(src)SetInitialConditions.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(obj)output.o: $(src)output.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(obj)SetParameters.o: $(src)SetParameters.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<


.PHONY: test1
test1: main.out
	./main.out --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 1 --type L

.PHONY: test2
test2: main.out
	./main.out --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 2 --type L

.PHONY: test3
test3: main.out
	./main.out --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 3 --type L

.PHONY: test4
test4: main.out
	./main.out --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 4 --type L

.PHONY: clean:
    rm -rf obj bin