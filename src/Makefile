# Makefile for compiling NSGA-III source code
CC = g++
LD = g++
RM = rm -f
OBJS:=$(patsubst %.cpp,%.o,$(wildcard *.cpp))
MAIN = nsga3

all: $(MAIN)

$(MAIN):main_nsga3.cpp
	$(LD) main_nsga3.cpp -o $(MAIN) -lm

clean:
	$(RM) $(MAIN) $(OBJS)
