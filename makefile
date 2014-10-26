DEBUG=no

ifeq ($(DEBUG),yes)
	CFLAGS=-g -c
else
	CFLAGS=-c
endif

EXEC=main
CXX=g++
LDFLAGS= 
SRC=src/probleme.cpp src/maillage.cpp src/main.cpp
OBJ=$(SRC:.cpp=.o)

all: ${EXEC}

ifeq ($(DEBUG),yes)
	@echo "generation en mode debug"
else
	@echo "generation en mode release"
endif

main: $(OBJ)
	@$(CXX) -o $@ $^ $(LDFLAGS)

main.o: ./include/probleme.h

%.o: %.cpp
	@$(CXX) $(CFLAGS) $< -o $@  -I ./eigen  



clean:
	@find -name '*~' -exec rm {} \;
	@find -name '*.o' -exec rm {} \;

