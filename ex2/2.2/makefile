CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.x : main.o random.o
	$(CC) random.o main.o -o main.x
main.o : main.c
	$(CC) -c main.c -o main.o $(CFLAGS)
random.o : rng/random.c rng/random.h
	$(CC) -c rng/random.c -o random.o $(CFLAGS)
clean :
	rm *.o main.x seed.out
