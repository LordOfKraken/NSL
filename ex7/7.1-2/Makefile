CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.x : Monte_Carlo_NVT.o random.o
	$(CC) random.o Monte_Carlo_NVT.o -o main.x
Monte_Carlo_NVT.o : Monte_Carlo_NVT.cpp
	$(CC) -c Monte_Carlo_NVT.cpp -o Monte_Carlo_NVT.o $(CFLAGS)
random.o : rng/random.c rng/random.h
	$(CC) -c rng/random.c -o random.o $(CFLAGS)
clean :
	rm *.o *.x
	rm -rf *gofr*
	rm -rf *gave*
	rm -rf *epot*
	rm -rf *pres*
	rm -rf frames/*.xyz
	rm -rf seed.out
	rm -rf config.final
