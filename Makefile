galsim: galsim.c 
	gcc -O3 -fopt-info-vec -ffast-math -Wall -fopenmp -march=native -o galsim galsim.c
clean:
	rm -f galsim
