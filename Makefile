
LIBS = -lgomp -lm
FLAGS = -fopenmp -DNUM_THREADS=4

all: capsule.o paralela.o mymath.o
	gcc -pg $? $(FLAGS) $(LIBS) -o ep

%.o: %.c
	gcc -pg $(FLAGS) -c $^ -o $@

clean:
	rm *.o; rm ep
