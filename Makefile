CC = gcc
CFLAGS = -O3
all:	hac
hac:	test-main.o hac.o hac.h
	$(CC) $(CFLAGS) -o hac test-main.o hac.o
test-main.o:	test-main.c hac.h
	$(CC) $(CFLAGS) -c test-main.c
hac.o:	hac.c
	$(CC) $(CFLAGS) -c hac.c
clean:
	rm -f hac test-main.o hac.o
