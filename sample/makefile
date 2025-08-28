
#LINUX
COMMON=-O2 -I../include -L../lib -pthread -Wl,-no-as-needed -Wl,-rpath,'$$ORIGIN'/../lib
LIBS = -lmujoco -lglfw
CC = g++


ROOT = basic

all:
	$(CC) $(COMMON) basic.cc $(LIBS) -o $(ROOT)

main.o:
	$(CC) $(COMMON) -c basic.cc

clean:
	rm *.o $(ROOT)
