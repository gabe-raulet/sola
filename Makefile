INCLUDES=-I/opt/homebrew/Cellar/portaudio/19.7.0/include -I/opt/homebrew/Cellar/portmidi/2.0.6/include -I/opt/homebrew/Cellar/libsndfile/1.2.2_1/include
LIBRARY=-L/opt/homebrew/Cellar/portaudio/19.7.0/lib -L/opt/homebrew/Cellar/portmidi/2.0.6/lib -L/opt/homebrew/Cellar/libsndfile/1.2.2_1/lib

D?=0
FLAGS=-lsndfile

ifeq ($(D),1)
FLAGS+=-O0 -g -fsanitize=address -fno-omit-frame-pointer
else
FLAGS+=-O2
endif

all: sola

sola: sola.c
	clang $(INCLUDES) $(LIBRARY) $(FLAGS) -o sola sola.c

clean:
	rm -rf *.dSYM *.o *.out sola
