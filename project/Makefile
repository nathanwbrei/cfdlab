# Include files
SOURCES=initLB.c visualLB.c boundary.c collision.c streaming.c computeCellValues.c flag.c main.c helper.c checks.c

# Compiler
# --------
CC=gcc

CFLAGS=-Wpedantic -Wall -std=c99 -O3 -g -lpthread -fopenmp -msse2 -ftree-vectorize -ftree-vectorizer-verbose=5 -fopt-info-vec-optimized

# Linker flags
# ------------
LDFLAGS = -lm -lpthread -fopenmp -lrt -ldl

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=lbsim

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) 

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)


$(OBJECTS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@
