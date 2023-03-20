default: main.out

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
	
main.out: main.cpp
	g++ -O3 -o main.out main.cpp -lblas -lboost_program_options -fopenmp
