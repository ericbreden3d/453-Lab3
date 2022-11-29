lab3: run

run: comp
	./lab3

comp: lab3.cpp Matrix.cpp Matrix.h
	mpicxx lab3.cpp Matrix.cpp -o lab3

