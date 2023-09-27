build:
	mpicxx -g -o sequential sequential.c
	mpicxx -fopenmp -g -c parallel.c
	nvcc -gencode arch=compute_61,code=sm_61 -g -c cuda.cu -o cuda.o
	mpicxx -fopenmp -g -o parallel parallel.o cuda.o -L/usr/local/cuda/lib64 -lcudart

run:
	mpirun -np 4 ./parallel 

clean:
	rm *.o sequential parallel