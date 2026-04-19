#
CC=gcc
COPTS= -O2 -Wall

main_mesh: mesh.o main_mesh.o eigen.o libmesh6.o
	$(CC) $(COPTS) -o main_mesh mesh.o main_mesh.o eigen.o libmesh6.o -lm

main_quality: mesh.o main_quality.o eigen.o libmesh6.o
	$(CC) $(COPTS) -o main_quality mesh.o main_quality.o eigen.o libmesh6.o -lm

cavity: mesh.o main_cavity_test.o eigen.o libmesh6.o
	$(CC) $(COPTS) -o cavity mesh.o main_cavity_test.o eigen.o libmesh6.o -lm

triangulation_naive: mesh.o main_triangulation_naive.o eigen.o libmesh6.o
	$(CC) $(COPTS) -o triangulation_naive mesh.o main_triangulation_naive.o eigen.o libmesh6.o -lm

triangulation_delaunay: mesh.o main_triangulation_delaunay.o eigen.o libmesh6.o
	$(CC) $(COPTS) -o triangulation_delaunay mesh.o main_triangulation_delaunay.o eigen.o libmesh6.o -lm

swap_test: mesh.o main_swap_test.o eigen.o libmesh6.o
	$(CC) $(COPTS) -o swap_test mesh.o main_swap_test.o eigen.o libmesh6.o -lm

libmesh6.o :	libmesh6.c libmesh6.h 
	$(CC) -c $(COPTS)  -I. libmesh6.c  

mesh.o :	mesh.c mesh.h libmesh6.h 
	$(CC) -c $(COPTS)  -I. mesh.c

eigen.o :	eigen.c mesh.h 
	$(CC) -c $(COPTS)  -I. eigen.c

main_mesh.o :	main_mesh.c mesh.h
	$(CC) -c $(COPTS)  -I. main_mesh.c 

main_quality.o :	main_quality.c mesh.h
	$(CC) -c $(COPTS)  -I. main_quality.c 

main_triangulation_naive.o: main_triangulation_naive.c mesh.h
	$(CC) -c $(COPTS)  -I. main_triangulation_naive.c

main_triangulation_delaunay.o: main_triangulation_delaunay.c mesh.h
	$(CC) -c $(COPTS)  -I. main_triangulation_delaunay.c

main_swap_test.o: main_swap_test.c mesh.h
	$(CC) -c $(COPTS)  -I. main_swap_test.c

main_cavity_test.o: main_cavity_test.c mesh.h
	$(CC) -c $(COPTS)  -I. main_cavity_test.c

clean :
	-rm main_mesh main_quality cavity triangulation_naive triangulation_delaunay swap_test mesh.o main_mesh.o main_quality.o main_triangulation_naive.o main_triangulation_delaunay.o eigen.o libmesh6.o main_swap_test.o main_cavity_test.o
