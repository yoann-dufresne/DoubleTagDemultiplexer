CC=g++
CFLAGS=-c -Wall -std=c++14 -O3 -lboost_filesystem
LDFLAGS=-lboost_filesystem -lboost_system -I /usr/include/boost
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
	