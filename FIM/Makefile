main: FIM_Ising.cpp
	g++ -fopenmp -Wall -std=c++11  -O3  FIM_Ising.cpp -o main -I/home/ccheng32/local/include -L/home/ccheng32/local/lib64 -larmadillo
clean:
	rm -rf main *.o
format:
	clang-format -i --style=google *.cpp *.h *.cu
