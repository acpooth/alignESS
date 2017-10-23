CC=g++
CFLAGS=-c -Wall  
LDFLAGS=/usr/lib/libboost_regex.a
SOURCES=cpp/main.cpp cpp/AG.cpp cpp/Alineamiento.cpp cpp/Poblacion.cpp cpp/Secuencia.cpp cpp/Individuo.cpp cpp/TFitness.cpp cpp/FitnessEntropy.cpp cpp/FitnessInfo.cpp cpp/FitnessSumaPares.cpp cpp/Selector.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=Ejecutable/AlineaMultiple

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
clean:
	rm -rf cpp/*o Ejecutable/AlineaMultiple
