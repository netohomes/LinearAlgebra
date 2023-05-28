#include "pincel.h"
using namespace std;														
//============================================================================
// PINCEL PARA PINTAR CON CAIRO
//............................................................................
pincel::pincel( ){
	// Fija el color a negro
	r = 0.0;
	g = 0.0;
	b = 0.0;
	a = 1.0;
	grosor = 10.0; 
	rT = gT = bT = 0.0;
	aT = 1.0;
	cr = NULL;
}
//............................................................................
void pincel::reinicia( GdkWindow *sup ){
	// Crea el pincel para cuando se emplea el constructor por default
	cr = gdk_cairo_create( sup );
}
//............................................................................
pincel::pincel( GdkWindow *sup ){
  	// Fija el color a negro
	r = 0.0;
	g = 0.0;
	b = 0.0;
	a = 1.0;
	grosor = 10.0; 
	rT = gT = bT = 0.0;
	aT = 1.0;
	// Crea el pincel
	cr = gdk_cairo_create( sup );
}
//............................................................................
pincel::pincel( GdkWindow *sup,  double rr, double gg, double bb, 
                double ggrosor ){
	r = rr;
	g = gg;
	b = bb;
	grosor = ggrosor;
	rT = gT = bT = 0.0;
	aT = 1.0;
	tamText = 1.0;
	tipoLetra = "serif";
	// Crea el pincel
	cr = gdk_cairo_create( sup );
}
//............................................................................
pincel::pincel( GdkWindow *sup,  double rr, double gg, double bb, 
                double aa, double ggrosor ){
	r = rr;
	g = gg;
	b = bb;
	a = aa;
	grosor = ggrosor;
	rT = gT = bT = 0.0;
	aT = 1.0;
	tamText = 1.0;
	tipoLetra = "serif";
	// Crea el pincel
	cr = gdk_cairo_create( sup );
}
//............................................................................
pincel::pincel( GdkWindow *sup, double rr, double gg, double bb,
                double aa, double ggrosor, double rrT, double ggT, double bbT,
                double aaT, double hText, string letra ){
	r = rr;
	g = gg;
	b = bb;
	a = aa;
	grosor = ggrosor;
	rT = rrT;
	gT = ggT;
	bT = bbT;
	aT = aaT;
	tamText   = hText;
	tipoLetra = letra;
	// Crea el pincel
	cr = gdk_cairo_create( sup );
}
//............................................................................
// Constructor por copia
pincel::pincel( const pincel &P ){
	cr = P.cr;
	r = P.r;
	g = P.g;
	b = P.b;
	a = P.a;
	grosor  = P.grosor;
	tamText = P.tamText;
	rT = P.rT;
	gT = P.gT;
	bT = P.bT;
	aT = P.aT;
	tipoLetra = P.tipoLetra;
}
//............................................................................
pincel::~pincel( ){
	if( cr != NULL ) cairo_destroy( cr );
}
//............................................................................
void pincel::dibujaLinea(  double x1, double y1,  double x2, double y2 ){
	cairo_set_source_rgba( cr, r, g, b, a );
	cairo_set_line_width( cr, grosor );
	cairo_move_to( cr, x1, y1 );
	cairo_line_to( cr, x2, y2 );
}
//............................................................................
void pincel::stroke( void ){
	// Pinta las línes que se han trazado
	cairo_stroke( cr );
}
//............................................................................
void pincel::dibujaTexto( double x, double y, string texto ){
	cairo_select_font_face( cr, tipoLetra.c_str( ), CAIRO_FONT_SLANT_NORMAL,
	                            CAIRO_FONT_WEIGHT_BOLD );
	cairo_set_font_size( cr, tamText );
	cairo_set_source_rgba( cr, rT, gT, bT, aT );
	cairo_move_to( cr, x, y );
	cairo_show_text( cr, texto.c_str( ) );
}
//............................................................................
void pincel::dibujaRectangulo( double x, double y, double largo, double ancho ){
	cairo_set_source_rgba( cr, r, g, b, a );
	cairo_move_to( cr, x, y );
	cairo_rectangle( cr, x, y, largo, ancho );
	cairo_fill( cr );
}	
//............................................................................
void pincel::dibujaPunto( double x, double y, double r ){	
	double pi = 3.1416;
	cairo_set_source_rgba( cr, r, g, b, a );
	cairo_move_to( cr, x, y );
    cairo_arc( cr, x, y, r, 0.0, 2.0*pi);
    cairo_fill(cr);
}
//............................................................................
//............................................................................
//============================================================================

