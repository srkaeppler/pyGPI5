CFLAGS = -Wall -fPIC -ansi -D_GNU_SOURCE -D_DARWIN
CC = gcc
OBJS = nrlmsise-00.o nrlmsise-00_data.o

all: clean nrlmsise-00.o nrlmsise-00_test.o nrlmsise-00_data.o msislib nrlmsise-test

nrlmsise-test: nrlmsise-00.o nrlmsise-00_test.o nrlmsise-00_data.o
	$(CC) $(CFLAGS) nrlmsise-00.o -lm nrlmsise-00_test.o nrlmsise-00_data.o -o nrlmsise-test

nrlmsise-00.o: nrlmsise-00.c nrlmsise-00.h
	$(CC) $(CFLAGS) -c nrlmsise-00.c -o nrlmsise-00.o

nrlmsise-00_test.o:	nrlmsise-00_test.c nrlmsise-00.h
	$(CC) $(CFLAGS) -c nrlmsise-00_test.c -o nrlmsise-00_test.o

nrlmsise-00_data.o:	nrlmsise-00_data.c nrlmsise-00.h
	$(CC) $(CFLAGS) -c nrlmsise-00_data.c -o nrlmsise-00_data.o

msislib: nrlmsise-00.o nrlmsise-00_data.o
	$(CC) -shared -o libnrlmsise-00.so nrlmsise-00.o nrlmsise-00_data.o

clean:
	rm -rf  nrlmsise-test nrlmsise-00.o nrlmsise-00_test.o nrlmsise-00_data.o libnrlmsise-00.so
