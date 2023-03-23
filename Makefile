bin = ./bin
src = ./src
obj = ./obj

CXX = g++
CXXFLAGS = -O3 -Wall -g -c

HDRS = $(src)/ShallowWater.h
OBJS = $(obj)/main.o $(obj)/TimeIntegrate.o $(obj)/derivatives.o $(obj)/SetInitialConditions.o $(obj)/output.o $(obj)/ConstructorDestructor.o
LIBS = -lblas -lboost_program_options -fopenmp

.PHONY = test1, test2, test3, test4, clean, doc

default: $(bin)/main.out

$(bin)/main.out: $(OBJS)
	@mkdir -p $(@D)
	$(CXX) -o $@ $^ $(LIBS)

$(obj)/%.o: $(src)/%.cpp $(HDRS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $@ -c $< -fopenmp

test1: $(bin)/main.out
	./bin/main.out --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 1 --type L

test2: $(bin)/main.out
	./bin/main.out --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 2 --type L

test3: $(bin)/main.out
	./bin/main.out --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 3 --type L

test4: $(bin)/main.out
	./bin/main.out --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 4 --type L

clean:
	rm -rf obj bin

doc:
	doxygen Doxyfile