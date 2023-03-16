default: main.out

main.out: main.cpp
	g++ -o main.out main.cpp -lboost_program_options