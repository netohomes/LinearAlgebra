#ifndef GRAFICADOR_H
#define GRAFICADOR_H
#include "../matrix.h"
#include "pincel.h"
#include "fparser/fparser.hh"
//============================================================================
class graficador{
	std::string var;
	matriz<double> coord;
	GdkWindow *sup;         // Superficie sobre la cual se pintará
	FunctionParser fparser; // parser para función
	double largo, alto;
	matriz<double> areaVista, areaGraf;
	int divEje, nPunt;
	double yMin, yMax, esc, uni;
	int tipoFunc, tipoEsc; // Tipo de escala 0->alto 1->largo
	double (*func) ( double );

 public:
	void grafica( double xMin, double xMax, int nPunt );
	void graficaFija( double xMin, double xMax, int nPunt );
	double escala( double yMin, double yMax );
	double escala( double xMin, double xMax, double yMin, double yMax );
	void generaPuntos( double n, double xMin, double xMax, double &yMin,
	                   double &yMax );
	void generaPuntos( double n, double xMin, double xMax );
	void dibujaGraf( int nPunt, double **coord, double **areaVista,	
                     double **areaGraf, double esc );
	void dibujaFunc( int n, double **coord, double **areaGraf,
					 double esc, pincel &pinFunc );
	void dibujaEjes( double **areaVista, double esc, pincel &pinEjes );
	void transforma( double x, double y, double xMin, double xMax,
	                 double yMin, double yMax, double esc, double &xp,
	                 double &yp );

	void funcionEval( std::string func, std::string x );
	void funcionEval( double (*func) ( double ) );
	double evaluaFuncion( const double val );

	void trazaCotasXEjes( double x1, double y1, double x2, double y2,
                      double uni, double esc, double **areaVista, double lcP,
	                  double lcS, pincel &cotasP, pincel &cotasS );
	void trazaCotasXBordes( double x1, double y1, double x2, double y2,
                            double uni, double esc, double **areaVista,
                            double lcP, double lcS, pincel &cotasP,
	                        pincel &cotasS, bool arriba );
	void trazaCotasYEjes( double x1, double y1, double x2, double y2,
                          double uni, double esc, double **areaVista,
                          double lcP, double lcS, pincel &cotasP,
	                      pincel &cotasS );
	void trazaCotasYBordes(double x1, double y1, double x2, double y2,
                           double uni, double esc, double **areaVista,
                           double lcP, double lcS, pincel &cotasP,
	                       pincel &cotasS, bool izq );
	double tamUnidad( double yMax, double yMin );
	double redondea( const double cval );

	void desplaza( int dir, double desp );
	void desplazaFija( int dir, double desp );
	void zoom( double f );
	void zoomFija( double f );

	void dibujaPuntos( int n, double *X, double *Y, pincel &pinPto );
	void dibujaPuntos( int n, double *X, double *Y );
	
	graficador( double L, double A, GdkWindow *areaDibujo );
	~graficador( );
};
//============================================================================
class data{
	// Almacena los datos para pasarlos a gtk
 public:
	std::string funcion, var;
	double largo, alto, a, b;
	graficador *func2D;
	double (*func)( double );
	int nPunt;
	double  *X, *Y;
	bool funcion_en_una_cadena;

	data( std::string F, std::string v, double L, double A, double aa,
	      double bb );
	data( double L, double A, double aa, double bb, double (*f)( double ), 
	      int n, double *X, double *Y );	
	~data( );
};
//============================================================================
void grafica( double a, double b, double (*func)( double ), int argc,
              char **argv, int L, int A );
void max_min_X( int n, double *X, double &a, double &b );
//============================================================================
#endif
