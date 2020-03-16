CC = gcc
CFLAGS = -Wall -pedantic -lm -O3

all: traveling

traveling: travelingMain.c travelingGreedy.c travelingMC.c
	$(CC) -o travelingSalesman travelingMain.c $(CFLAGS)
	
clean:
	rm -f travelingSalesman
