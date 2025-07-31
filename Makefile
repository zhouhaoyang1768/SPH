CC         = g++
EXECUTABLE = sph
CFLAGS     = -c -Wall -O3 -fopenmp

SOURCES    = simulator.cpp SPH2.cpp
OBJECTS    = $(SOURCES:.cpp=.o)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -fopenmp -lGL -lGLU -lglut -lstdc++ -o $(EXECUTABLE)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o $(EXECUTABLE)
