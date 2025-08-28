
# #MAC
# COMMON=-O2 -I../../include -L../../bin -mavx -pthread
# LIBS = -w -lmujoco200 -lglfw.3
# CC = gcc

#LINUX
# COMMON=-O2 -I../../include -L../../bin -mavx -pthread -Wl,-rpath,'$$ORIGIN'
# LIBS = -lmujoco200 -lGL -lm -lglew ../../bin/libglfw.so.3
# CC = gcc

COMMON=-O2 -I../../include -L../../lib -mavx -pthread -Wl,-no-as-needed -Wl,-rpath,'$$ORIGIN'/../../lib
LIBS = -lmujoco -lglfw -lm
CC = gcc

#WINDOWS
#COMMON=/O2 /MT /EHsc /arch:AVX /I../../include /Fe../../bin/
#LIBS = ../../bin/glfw3.lib  ../../bin/mujoco200.lib
#CC = cl

ROOT = cart_pend

all:
	$(CC) $(COMMON) main.c $(LIBS) -o $(ROOT)
# 	$(CC) $(COMMON) main.c $(LIBS) -o ../../bin/$(ROOT)
	

main.o:
	$(CC) $(COMMON) -c main.c

clean:
	rm *.o $(ROOT)
