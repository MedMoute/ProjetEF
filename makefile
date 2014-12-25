DEBUG=yes

ifeq ($(DEBUG),yes)
	CFLAGS=-g -c
else
	CFLAGS=-c
endif

EXEC=main
CXX=mpicxx
LDFLAGS= 
SRC=src/parallel.cpp src/nonParallel.cpp src/probleme.cpp src/maillage.cpp src/main.cpp
OBJ=$(SRC:.cpp=.o)

all: ${EXEC}

ifeq ($(DEBUG),yes)
	@echo "generation en mode debug"
else
	@echo "generation en mode release"
endif

main: $(OBJ)
	@$(CXX) -o $@ $^ $(LDFLAGS)

main.o: ./include/param.h

%.o: %.cpp
	@$(CXX) $(CFLAGS) $< -o $@  -I ./eigen  



clean:
	@find -name '*~' -exec rm {} \;
	@find -name '*.o' -exec rm {} \;

