default: main.out

main.out: main.cpp
	g++ -O3 -o main.out main.cpp -lboost_program_options