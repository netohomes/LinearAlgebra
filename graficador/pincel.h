#ifndef PINCEL_H
#define PINCEL_H
#include <cairo/cairo.h>
#include <iostream>
#include <gtk/gtk.h>
#include <gdk/gdk.h>
//============================================================================
class pincel{
	cairo_t *cr;
public:
	// Propiedades para las líneas
	double r;    // Rojo
	double g;    // Verde
	double b;    // Azul
	double a;    // Transparencia
	double grosor;
	// Propiedades para el texto
	double tamText;
	double rT, gT, bT, aT;    // Colores
	std::string tipoLetra;	
	pincel( );
	pincel( GdkWindow* );
	pincel( GdkWindow*, double, double, double, double );
	pincel( GdkWindow*, double, double, double, double, double );
	pincel( GdkWindow*, double, double, double, double, double,
            double, double, double, double, double, std::string );
	pincel( const pincel& );
	void reinicia( GdkWindow* );
	void dibujaLinea( double, double, double, double );
	void stroke( void );
	void dibujaTexto( double, double, std::string );
	void dibujaRectangulo( double, double, double, double );
	void dibujaPunto( double x, double y, double r );
	~pincel( ); 
};
//============================================================================
#endif
