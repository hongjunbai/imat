all: matrix_demo vector_demo

#CC=icpc
#flags=-fast
CC = g++ 
flags= -O3 #-g

# !!! After GSL installed, check and modify this accordingly
#gsl= -DGSL -lgsl -lgslcblas -lm #-I/opt/local/include -L/opt/local/lib

vector_demo: demos/vector_demo.cpp imat/vector.h imat/vector_inl.h imat/base.h imat/base_inl.h
	$(CC) $(flags) demos/vector_demo.cpp -o vector_demo
matrix_demo: demos/matrix_demo.cpp imat/matrix.h imat/matrix_inl.h imat/base.h imat/base_inl.h
	$(CC) $(flags) demos/matrix_demo.cpp -o matrix_demo $(gsl)

clean:
	/bin/rm -f vector_demo matrix_demo.dat matrix_demo
