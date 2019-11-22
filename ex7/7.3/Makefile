CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.x : MolDyn_NVE.o #random.o
	$(CC) MolDyn_NVE.o -o main.x
MolDyn_NVE.o : main.cpp
	$(CC) -c main.cpp -o MolDyn_NVE.o $(CFLAGS)
#random.o : rng/random.c rng/random.h
#	$(CC) -c rng/random.c -o random.o $(CFLAGS)
clean :
	rm *.o *.x
	rm -rf *epot*
	rm -rf *ekin*
	rm -rf *etot*
	rm -rf *temp*
	rm -rf frames/*.xyz
