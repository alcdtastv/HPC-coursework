default: main.out

.PHONY: runblas
runblas: main.out
	./main.out --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 3 --type B

.PHONY: runloop
runloop: main.out
	./main.out --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 3 --type L
	
main.out: main.cpp
	g++ -O3 -o main.out main.cpp -lblas -lboost_program_options