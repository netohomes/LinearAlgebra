# VARIABLES DE COMPILACIÓN
# Compilador
CC = g++
# Opciones de compilación
CFLAGS = -Wall -g
# Ejecutable
EXE = graficador

all:
	$(CC) $(CFLAGS) main.cpp graficador.cpp pincel.cpp ../extra/matrix.cpp fparser/fparser.cc `pkg-config --cflags --libs gtk+-3.0` -o $(EXE)

clean:
	@echo Cleanup
	rm $(EXE)


