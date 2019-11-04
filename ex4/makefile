CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.x : main.o
	$(CC) main.o -o main.x
main.o : main.c
	$(CC) -c main.c -o main.o $(CFLAGS)
clean :
	rm *.o main.x
	rm -rf *epot*
	rm -rf *ekin*
	rm -rf *etot*
	rm -rf *temp*
	rm -rf *pres*
	rm -rf frames/*.xyz
