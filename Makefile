
LIBS = -lgomp -lm
ARCH = native
NUM_CORES = `cat /proc/cpuinfo | grep processor | wc -l`
DEBUG = OFF

ifeq ($(DEBUG), OFF)
	FLAGS = -mtune=$(ARCH) -O3 -fopenmp -DNUM_THREADS=$(NUM_CORES) # -pg -DDEBUG -floop-block -floop-interchange
else # DEBUG == ON
	FLAGS = -O0 -fopenmp -DNUM_THREADS=1 -pg -DDEBUG
endif

all: capsule.o paralela.o mymath.o
	gcc $? $(FLAGS) $(LIBS) -o ep

%.o: %.c
	gcc $(FLAGS) -c $^ -o $@

clean:
	rm *.o; rm ep
