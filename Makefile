HOST=$(shell hostname)

include Makefile.$(HOST)

targets= fastn2mex.mexa64 # brent bisection newton secant steffenson 

all: $(targets)

brent: brent.cpp
	g++ -std=c++11 -fopenmp -L/home/honglin/local/lib/ -I/home/honglin/local/include/ ./brent.cpp -lgsl -lgslcblas -o $@

bisection: bisection.cpp
	g++ -std=c++11 -fopenmp -L/home/honglin/local/lib/ -I/home/honglin/local/include/ ./bisection.cpp -lgsl -lgslcblas -o $@

newton: newton.cpp
	g++ -std=c++11 -fopenmp -L/home/honglin/local/lib/ -I/home/honglin/local/include/ ./newton.cpp -lgsl -lgslcblas -o $@

secant: secant.cpp
	g++ -std=c++11 -fopenmp -L/home/honglin/local/lib/ -I/home/honglin/local/include/ ./$@.cpp -lgsl -lgslcblas -o $@

steffenson: steffenson.cpp
	g++ -std=c++11 -fopenmp -L/home/honglin/local/lib/ -I/home/honglin/local/include/ ./$@.cpp -lgsl -lgslcblas -o $@

%.o: %.f
	g++ -fPIC -O3 -c $^

fastn2mex.o: fastn2mex.cpp fzero.cpp
	g++ -std=c++11 -Wall -c -fPIC -O3 -I$(MATLAB_HOME)/extern/include/ fastn2mex.cpp

fastn2mex.mexa64: fastn2mex.o $(LAPACK) $(BLAS) dpofa.o dposl.o
	echo $(HOST)
	$(MATLAB_HOME)/bin/mex $^ -lgfortran

clean:
	rm *.o
	rm $(targets)
