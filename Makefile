IFlags = -I$(HOME)/sw/include
LFlags = -L$(HOME)/sw/lib

all:
	g++ ${IFlags} RGSW.cpp -O3 -o RGSW_output ${LFlags} -lntl -lgmp -lstdc++
	
	./RGSW_output