# VARIABLES DE COMPILACIÓN
# Compilador
CC = g++
# Opciones de compilación
FLAGS = -Wall -g
# Nombre del programa
EXE = output

all:
	@echo Compilando...
	$(CC) $(FLAGS) -I./usr/local/Cellar/cairo/1.16.0_5/include/cairo main.cpp eigen_val_vec.cpp derivatives_and_integrals.cpp linear_solvers.cpp metInterpolacion.cpp graficador/graficador.cpp graficador/pincel.cpp datos.cpp graficador/fparser/fparser.cc graficador/grafFunc.cpp `pkg-config --cflags --libs gtk+-3.0` -o $(EXE)
clean:
	@echo Borrando ejecutable...
	rm $(EXE)