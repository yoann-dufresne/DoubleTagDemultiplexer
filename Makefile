CC=g++
CFLAGS=-c -Wall -std=c++14 -O3
LDFLAGS=-lstdc++fs
SOURCES=$(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=dtd

all: $(SOURCES) $(EXECUTABLE)
    
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	@rm -rf *.o

mrproper: clean
	@rm -rf $(EXEC)
	