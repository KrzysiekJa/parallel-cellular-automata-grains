CC = g++
CFLAGS  = -std=c++17 -O3 -Wall -pedantic
MPFLAGS = -Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include/
LIB = -L/usr/local/opt/llvm/lib -lomp
HDS = headers
CLS = classes
UTLS = utils
OBJECTS = main.o grainsgrowth.o utils_functions.o cell.o


all: simula clean run

variables.o: $(HDS)/variables.hpp

cell.o: $(CLS)/cell.cpp $(HDS)/cell.hpp $(HDS)/variables.hpp
	$(CC) $(CFLAGS) -c $(CLS)/cell.cpp

grainsgrowth.o: $(CLS)/grainsgrowth.cpp $(HDS)/grainsgrowth.hpp $(HDS)/variables.hpp $(HDS)/cell.hpp
	$(CC) $(CFLAGS) -c $(CLS)/grainsgrowth.cpp $(MPFLAGS)

utils_functions.o: $(UTLS)/utils_functions.cpp $(UTLS)/utils_functions.hpp $(HDS)/variables.hpp $(HDS)/cell.hpp
	$(CC) $(CFLAGS) -c $(UTLS)/utils_functions.cpp

main.o: main.cpp $(UTLS)/utils_functions.cpp $(HDS)/grainsgrowth.hpp $(HDS)/variables.hpp $(HDS)/cell.hpp
	$(CC) $(CFLAGS) -c main.cpp

simula: $(OBJECTS)
	$(CC) $(OBJECTS) -o simula $(LIB)

clean:
	rm -f *.o
	echo Clean done

run:
	./simula
	echo Simulation done
