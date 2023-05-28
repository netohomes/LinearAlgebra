#include <iostream>
#include "graficador.h"
#include <iomanip>
using namespace std;
#define NUM_DIGITOS 5
//============================================================================
// DATOS
//============================================================================
data::data( string F, string v, double L, double A, double aa, double bb ):
            funcion( F ), var( v ), largo( L ), alto( A ), a( aa ), b( bb ) {
	func2D = NULL;
	func   = NULL;
	funcion_en_una_cadena = true;
	X = NULL;
	Y = NULL;
}
data::data( double L, double A, double aa, double bb, double (*f)( double ),
            int n, double *Xp, double *Yp ): largo( L ), alto( A ),
            a( aa ), b( bb ), nPunt( n ), X( Xp ), Y( Yp ){
	func2D = NULL;
	func   = f;
	funcion_en_una_cadena = false;
}
data::~data( ){ }
//============================================================================
// GRAFICADOR DE FUNCIONES 2D
//============================================================================
//............................................................................
graficador::graficador( double L, double A, GdkWindow *areaDibujo ):
                        coord( ), areaVista( 4, 2 ), areaGraf( 4, 2 ){
	largo = L;
	alto  = A;
	sup = areaDibujo;
	divEje = 15;
	tipoFunc = -1;
}
//............................................................................
void graficador::grafica( double xMin, double xMax, int nPunt ){
	// Dibuja una función en un intervalo definido por las coordenadas "y"
	// máxima y mínima en un área de dibujo definida
	//double esc, yMin, yMax;
	double aux;
	// Genera puntos
	this->nPunt = nPunt;
	// Pide el tamaño para la matriz de coordenadas
	coord.reinicia( nPunt, 2 );
	generaPuntos( nPunt, xMin, xMax, yMin, yMax );
	// Calcula la escala
	esc = escala( yMin, yMax );
	// Límites
	// Del área visible
	aux = 0.5 * ( largo / esc - ( xMax - xMin ) );
	areaVista.val[0][0] = xMin - aux;
	areaVista.val[0][1] = yMin;
	areaVista.val[1][0] = xMax + aux;
	areaVista.val[1][1] = yMin;
	areaVista.val[2][0] = xMax + aux;
	areaVista.val[2][1] = yMax;	
	areaVista.val[3][0] = xMin - aux;
	areaVista.val[3][1] = yMax;
	// Del área con parte de la función (en un incio coinciden)
	areaGraf.val[0][0] = xMin;
	areaGraf.val[0][1] = yMin;
	areaGraf.val[1][0] = xMax;
	areaGraf.val[1][1] = yMin;
	areaGraf.val[2][0] = xMax;
	areaGraf.val[2][1] = yMax;	
	areaGraf.val[3][0] = xMin;
	areaGraf.val[3][1] = yMax;
	// Dibuja la gráfica
	dibujaGraf( nPunt, coord.val, areaVista.val, areaGraf.val, esc );
}
//............................................................................
void graficador::funcionEval( string func, string x ){
	// Para cuando la función es pasada en una cadena
	fparser.Parse( func, x );
	tipoFunc = 0;
}
void graficador::funcionEval( double (*func) ( double ) ){
	// Cuando se le pasa un puntero a una función
	this->func = func;
	tipoFunc = 1;
}
//............................................................................
void graficador::desplaza( int dir, double desp ){
	// Desplaza la gráfica
	//double yMin, yMax, esc;
	// Cambio del área que se puede ver y el área visible de la gráfica
	// Desplazamiento vertical (se mueven las dos áreas igual (no se debería
	// hace por los posibles errores de desbordamiento))
	// Desplazamineto vertical
	if( dir == 0 ){
		// Del área visible
		areaVista.val[0][1] += desp;
		areaVista.val[1][1] += desp;
		areaVista.val[2][1] += desp;	
		areaVista.val[3][1] += desp;
	}
	// Desplazamiento horizontal
	else{
		// Del área visible
		areaVista.val[0][0] += desp;
		areaVista.val[1][0] += desp;
		areaVista.val[2][0] += desp;
		areaVista.val[3][0] += desp;
		// Del área con parte de la función (en un incio coinciden)
		areaGraf.val[0][0] += desp;
		areaGraf.val[1][0] += desp;
		areaGraf.val[2][0] += desp;
		areaGraf.val[3][0] += desp;
	}
	// Nuevos puntos
	generaPuntos( nPunt, areaGraf.val[0][0], areaGraf.val[1][0] );

	// Dibuja la gráfica
	dibujaGraf( nPunt, coord.val, areaVista.val, areaGraf.val, esc );	
}
//............................................................................
void graficador::zoom( double f ){
	double aux1;
	// Calcula la escala
	esc = escala( areaGraf.val[3][1] + f, areaGraf.val[0][1] - f );
	aux1 = f * largo / alto;
	// Límites
	// Del área visible
	areaVista.val[0][0] -= aux1;
	areaVista.val[0][1] -= f;
	areaVista.val[1][0] += aux1;
	areaVista.val[1][1] -= f;
	areaVista.val[2][0] += aux1;
	areaVista.val[2][1] += f;
	areaVista.val[3][0] -= aux1;
	areaVista.val[3][1] += f;
	// Del área con parte de la función (en un incio coinciden)
	areaGraf.val[0][0] -= aux1;
	areaGraf.val[0][1] -= f;
	areaGraf.val[1][0] += aux1;
	areaGraf.val[1][1] -= f;
	areaGraf.val[2][0] += aux1;
	areaGraf.val[2][1] += f;
	areaGraf.val[3][0] -= aux1;
	areaGraf.val[3][1] += f;
	// Genera puntos
	generaPuntos( nPunt, areaGraf.val[0][0], areaGraf.val[1][0] );
	// Dibuja la gráfica
	dibujaGraf( nPunt, coord.val, areaVista.val, areaGraf.val, esc );
}
//............................................................................
double graficador::escala( double yMin, double yMax ){
	// Determina la escala
	tipoEsc = 0; // Respecto al alto del área visible
	return alto / fabs( yMax - yMin  );
}
//............................................................................
void graficador::generaPuntos( double n, double xMin, double xMax,
                               double &yMin, double &yMax ){
	// Genera los puntos para hacer la gráfica
	double dx, x, y;
	// Tamaño de paso
	dx = fabs( xMax - xMin ) / static_cast < double > ( n - 1 );
	// Genera los puntos
	x = xMin;
	yMax = yMin = evaluaFuncion( x );	
	for( int i = 0; i < n; i++ ){
		y = evaluaFuncion( x );
		// Guarda los puntos
		coord.val[i][0] = x;
		coord.val[i][1] = y;
		if( y > yMax ) yMax = y;
		if( y < yMin ) yMin = y;
		// Avanza
		x += dx;
	}
	x = xMax;
	y = evaluaFuncion( x );
	coord.val[nPunt-1][0] = x;
	coord.val[nPunt-1][1] = y;
	if( y > yMax ) yMax = y;
	if( y < yMin ) yMin = y;
}
void graficador::generaPuntos( double n, double xMin, double xMax  ){
	// Genera los puntos para hacer la gráfica
	double dx, x;
	// Tamaño de paso
	dx = fabs( xMax - xMin ) / static_cast < double > ( n - 1 );
	// Pide el tamaño para la matriz de coordenadas
	coord.reinicia( n, 2 );
	// Genera los puntos
	x = xMin;
	for( int i = 0; i < n; i++ ){
		// Guarda los puntos
		coord.val[i][0] = x;
		coord.val[i][1] = evaluaFuncion( x );
		// Avanza
		x += dx;
	}
	coord.val[nPunt-1][0] = xMax;
	coord.val[nPunt-1][1] = evaluaFuncion( xMax );
}
//............................................................................
void graficador::dibujaGraf( int nPun, double **coord, double **areaVista,
                             double **areaGraf, double esc ){
	pincel pinFunc ( sup, 0.13, 0.55, 0.0, 2.5, 2.5 );
	pincel pinFondo( sup, 1.0, 1.0, 1.0, 1.0, 1.0 );
	pincel pinEjes( sup, 0.0, 0.0, 0.0, 3.0, 1.0 );
	pinFondo.dibujaRectangulo( 0.0, 0.0, largo, alto );
	// Dibuja la función
	dibujaFunc( nPunt, coord, areaVista, esc, pinFunc );
	// Dibuja ejes
	dibujaEjes( areaVista, esc, pinEjes );
	// Pinta todo
	pinFunc.stroke( );
	pinEjes.stroke( );
}
//............................................................................
void graficador::dibujaFunc( int nPunt, double **coord, double **areaVista,
                             double esc, pincel &pinFunc ){
	// Traza líneas entre los puntos
	double x1, y1, x2, y2;
	for( int i = 1; i < nPunt; i++ ){
		transforma( coord[i-1][0], coord[i-1][1], areaVista[0][0],
                    areaVista[1][0], areaVista[0][1], areaVista[3][1], esc,
		            x1, y1 );
		transforma( coord[i][0], coord[i][1], areaVista[0][0],
                    areaVista[1][0], areaVista[0][1], areaVista[3][1], esc,
		            x2, y2 );
	 	pinFunc.dibujaLinea( x1, y1, x2, y2 );
	}
}
//............................................................................
void graficador::dibujaEjes( double **areaVista, double esc, pincel &pinEjes ){
	// Traza los ejes con sus cotas
	double x1, y1, x2, y2, lcP, lcS;
	pincel cotasP( sup, 0.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 1.0, 12.0,
	               "serif" );
	pincel cotasS( sup, 0.5, 0.5, 0.5, 1.0, 1.0 );
	uni = tamUnidad( areaVista[0][1], areaVista[3][1] );
	lcP = 4.0;
	lcS = 4.0;
	// Eje "x"
	// Cotas sobre los ejes
	if( areaVista[0][1] * areaVista[3][1] < 0.0 ){
		transforma( areaVista[0][0], 0.0, areaVista[0][0], areaVista[1][0],
		            areaVista[0][1], areaVista[3][1], esc, x1, y1 );
		transforma( areaVista[1][0], 0.0, areaVista[0][0], areaVista[1][0],
		            areaVista[0][1], areaVista[3][1], esc, x2, y2 );
	 	pinEjes.dibujaLinea( x1, y1, x2, y2 );
		// Cotas
		trazaCotasXEjes( x1, y1, x2, y2, uni, esc, areaVista, lcP, lcS,
		                 cotasP, cotasS );
	}
	// Cotas sobre los bordes
	else{
		// Borde superior
		if( areaVista[3][1] <= 0.0 ){
			transforma( areaVista[3][0], areaVista[3][1], areaVista[0][0],
		                areaVista[1][0], areaVista[0][1], areaVista[3][1],
			            esc, x1, y1 );
			transforma( areaVista[2][0], areaVista[2][1], areaVista[0][0],
			            areaVista[1][0], areaVista[0][1], areaVista[3][1],
			            esc, x2, y2 );
			trazaCotasXBordes( x1, y1, x2, y2, uni, esc, areaVista, lcP, lcS,
		                       cotasP, cotasS, true );
		}
		// Borde inferior
		else{
			transforma( areaVista[0][0], areaVista[0][1], areaVista[0][0],
		                areaVista[1][0], areaVista[0][1], areaVista[3][1],
			            esc, x1, y1 );
			transforma( areaVista[1][0], areaVista[1][1], areaVista[0][0],
			            areaVista[1][0], areaVista[0][1], areaVista[3][1],
			            esc, x2, y2 );
			trazaCotasXBordes( x1, y1, x2, y2, uni, esc, areaVista, lcP, lcS,
		                       cotasP, cotasS, false );
		} 	
	}
	// Eje "y"
	if( areaVista[0][0] * areaVista[1][0] < 0.0 ){
		transforma( 0.0, areaVista[0][1], areaVista[0][0], areaVista[1][0],
		            areaVista[0][1], areaVista[3][1], esc, x1, y1 );
		transforma( 0.0, areaVista[3][1], areaVista[0][0], areaVista[1][0],
		            areaVista[0][1], areaVista[3][1], esc, x2, y2 );
	 	pinEjes.dibujaLinea( x1, y1, x2, y2 );
		// Cotas
		trazaCotasYEjes( x1, y1, x2, y2, uni, esc, areaVista, lcP, lcS,
		                 cotasP, cotasS );		
	}
	else{
		// Borde izquierdo
		if( areaVista[0][0] >= 0.0 ){
			transforma( areaVista[0][0], areaVista[0][1], areaVista[0][0],
			            areaVista[1][0], areaVista[0][1], areaVista[3][1],
			            esc, x1, y1 );
			transforma( areaVista[3][0], areaVista[3][1], areaVista[0][0],
			            areaVista[1][0], areaVista[0][1], areaVista[3][1],
			            esc, x2, y2 );
			// Cotas
			trazaCotasYBordes( x1, y1, x2, y2, uni, esc, areaVista, lcP, lcS,
		                       cotasP, cotasS, true );
		}
		// Borde derecho
		else{
			transforma( areaVista[1][0], areaVista[1][1], areaVista[0][0],
			            areaVista[1][0], areaVista[0][1], areaVista[3][1],
			            esc, x1, y1 );
			transforma( areaVista[2][0], areaVista[2][1], areaVista[0][0],
			            areaVista[1][0], areaVista[0][1], areaVista[3][1],
			            esc, x2, y2 );
			// Cotas
			trazaCotasYBordes( x1, y1, x2, y2, uni, esc, areaVista, lcP, lcS,
		                       cotasP, cotasS, false );
		}
	}
	cotasS.stroke( );
	cotasP.stroke( );
}
//............................................................................
void graficador::transforma( double x, double y, double xMin, double xMax,
                             double yMin, double yMax, double esc, double &xp,
                             double &yp ){
	// Transforma la ubicación de un punto a la correcta para graficarlo
	xp =  ( x - xMin ) * esc;
	yp =  ( -( y - yMin ) + fabs( yMax - yMin ) ) * esc;
}
//............................................................................
double graficador::evaluaFuncion( const double val ){
  	// Evalúa la función
  	double Fval, vals[] = { 0.0 };
	switch( tipoFunc ){
		// Con fparser
		case 0:
			vals[0] += val;
			Fval = fparser.Eval( vals );
			break;
		// Con otra función pasada a través de un puntero 
		case 1:
			Fval = func( val );
			break;
		default:
			cout << endl;
			cout << "ERROR: no has definido una funcion" << endl;
	}
	return Fval;
}
//............................................................................
void graficador::trazaCotasXEjes( double x1, double y1, double x2, double y2,
                                  double uni, double esc, double **areaVista,
                                  double lcP, double lcS, pincel &cotasP,
                                  pincel &cotasS ){
	// Todas las variables están en el sistema para dibujar y la unidad
	// también
	// x1, y1, x2, y2 YA ESTÁN en el sistema para dibujar
	double xR, xD, uniDir, res;
	int conta;
	// La unidad con la dirección adecuada
	uniDir = areaVista[0][0] / fabs( areaVista[0][0] ) * uni;
	xR = 0.0;
	conta  = 0;
	// Recorre xR hasta llegar antes del fin del área de vista
	while( fabs( xR ) < fabs( areaVista[0][0] ) ){
		xR += uniDir;
		if( areaVista[0][0] < 0.0 ) conta--;
		else                        conta++;
	}
	if( areaVista[0][0] < 0.0 ){
		xR -= uniDir; // Resta la iteración que se hizo de más
		conta++;
	}
	// Traza las cotas
	res = fabs( xR - areaVista[0][0] );
	xD = x1 + res * esc; // x1 ya está en el sistema de dibujo
	// Sobre el eje
	while( xR <= areaVista[1][0] ){
		// Cotas principales
		if( fabs( conta % 5 ) == 0  ){
			ostringstream s; string num;
			s << setprecision( NUM_DIGITOS ) << xR;
			num = s.str( );
			cotasP.dibujaLinea( xD, y1 + lcP, xD, y1 - lcP );
			cotasP.dibujaTexto( xD, y1 - lcP - 3.0, num ); 
		}
		// Cotas secundarias
		else{
			cotasS.dibujaLinea( xD, y1 + lcS, xD, y1 - lcS );
		}
		xR += uni; // Ojo: aquí es "uni" y no "uniDir"
		xD += uni * esc;
		conta++;
	}
}
//............................................................................
void graficador::trazaCotasXBordes( double x1, double y1, double x2, double y2,
                                    double uni, double esc, double **areaVista,
                                    double lcP, double lcS, pincel &cotasP,
	                                pincel &cotasS, bool arriba ){
	// Todas las variables están en el sistema para dibujar y la unidad
	// también
	// x1, y1, x2, y2 YA ESTÁN en el sistema para dibujar
	double xR, xD, uniDir, res;
	int conta;
	// La unidad con la dirección adecuada
	uniDir = areaVista[0][0] / fabs( areaVista[0][0] ) * uni;
	xR = 0.0;
	conta  = 0;
	// Recorre xR hasta llegar antes del fin del área de vista
	while( fabs( xR ) < fabs( areaVista[0][0] ) ){
		xR += uniDir;
		if( areaVista[0][0] < 0.0 ) conta--;
		else                        conta++;
	}
	if( areaVista[0][0] < 0.0 ){
		xR -= uniDir; // Resta la iteración que se hizo de más
		conta++;
	}
	// Traza las cotas
	res = fabs( xR - areaVista[0][0] );
	xD = x1 + res * esc; // x1 ya está en el sistema de dibujo
	while( xR <= areaVista[1][0] ){
		// Cotas principales
		if( fabs( conta % 5 ) == 0  ){
			ostringstream s; string num;
			s << setprecision( NUM_DIGITOS ) << xR;
			num = s.str( );
			// Borde superior
			if( arriba ){
				cotasP.dibujaLinea( xD, 0.0, xD, 0.0 + lcP );
				cotasP.dibujaTexto( xD, 0.0 + lcP + 10.0, num ); 
			}
			// Borde inferior
			else{
				cotasP.dibujaLinea( xD, alto, xD, alto - lcP );
				cotasP.dibujaTexto( xD, alto - lcP - 5.0, num ); 
			}
		}
		// Cotas secundarias
		else{
			// Borde superior
			if( arriba ){
				cotasS.dibujaLinea( xD, 0.0, xD, 0.0 + lcS );
			}
			// Borde inferior
			else{
				cotasS.dibujaLinea( xD, alto, xD, alto - lcS );
			}
		}
		xR += uni; // Ojo: aquí es "uni" y no "uniDir"
		xD += uni * esc;
		conta++;
	}
}
//............................................................................
void graficador::trazaCotasYEjes( double x1, double y1, double x2, double y2,
                                  double uni, double esc, double **areaVista,
                                  double lcP, double lcS, pincel &cotasP,
	                              pincel &cotasS ){
	// Todas las variables están en el sistema para dibujar y la unidad
	// también
	// x1, y1, x2, y2 YA ESTÁN en el sistema para dibujar
	double yR, yD, uniDir, res;
	int conta;
	// La unidad con la dirección adecuada
	uniDir = areaVista[0][1] / fabs( areaVista[0][1] ) * uni;
	yR = 0.0;
	conta = 0;
	// Recorre yR hasta llegar antes del fin del área de vista
	while( fabs( yR ) < fabs( areaVista[0][1] ) ){
		yR += uniDir;
		if( areaVista[0][1] < 0.0 ) conta--;
		else                        conta++;
	}
	if( areaVista[0][1] < 0.0 ){
		yR -= uniDir; // Resta la iteración que se hizo de más
		conta++;
	}
	// Traza las cotas
	res = fabs( yR - areaVista[0][1] );
	yD = y1 - res * esc; // x1 ya está en el sistema de dibujo
	// Sobre el eje
	while( yR <= areaVista[3][1] ){
		// Cotas principales
		if( fabs( conta % 5 ) == 0  ){
			ostringstream s; string num;
			s << setprecision( NUM_DIGITOS ) << yR;
			num = s.str( );
			cotasP.dibujaLinea( x1 - lcP, yD, x1 + lcP, yD );
			cotasP.dibujaTexto( x1 + lcP + 3.0, yD, num ); 
		}
		// Cotas secundarias
		else{
			cotasS.dibujaLinea( x1 - lcP, yD, x1 + lcP, yD );
		}
		yR += uni; // Ojo: aquí es "uni" y no "uniDir"
		yD -= uni * esc;
		conta++;
	}
}
//............................................................................
void graficador::trazaCotasYBordes(double x1, double y1, double x2, double y2,
                                  double uni, double esc, double **areaVista,
                                  double lcP, double lcS, pincel &cotasP,
								   pincel &cotasS, bool izq ){
	// Todas las variables están en el sistema para dibujar y la unidad
	// también
	// x1, y1, x2, y2 YA ESTÁN en el sistema para dibujar
	double yR, yD, uniDir, res;
	int conta;
	// La unidad con la dirección adecuada
	uniDir = areaVista[0][1] / fabs( areaVista[0][1] ) * uni;
	yR = 0.0;
	conta = 0;
	// Recorre yR hasta llegar antes del fin del área de vista
	while( fabs( yR ) < fabs( areaVista[0][1] ) ){
		yR += uniDir;
		if( areaVista[0][1] < 0.0 ) conta--;
		else                        conta++;
	}
	if( areaVista[0][1] < 0.0 ){
		yR -= uniDir; // Resta la iteración que se hizo de más
		conta++;
	}
	// Traza las cotas
	res = fabs( yR - areaVista[0][1] );
	yD = alto; // x1 ya está en el sistema de dibujo
	while( yR <= areaVista[3][1] ){
		// Cotas principales
		if( fabs( conta % 5 ) == 0  ){
			ostringstream s; string num;
			s << setprecision( NUM_DIGITOS ) << yR;
			// A la izquierda
			num = s.str( );
			if( izq ){
				cotasP.dibujaLinea( 0.0, yD, 0.0 + lcP, yD );
				cotasP.dibujaTexto( 0.0 + lcP + 5.0, yD, num );
			}
			// A la derecha
			else{
				cotasP.dibujaLinea( largo - lcP, yD, largo, yD );
				cotasP.dibujaTexto( largo - lcP - 25.0, yD, num );
			}
		}
		// Cotas secundarias
		else{
			// A la izquierda
			if( izq ) cotasS.dibujaLinea( 0.0, yD, 0.0 + lcP, yD );
			// A la derecha
			else      cotasS.dibujaLinea( largo - lcP, yD, largo, yD );
		}
		yR += uni; // Ojo: aquí es "uni" y no "uniDir"
		yD -= uni * esc;
		conta++;
	}
}
//............................................................................
double graficador::tamUnidad( double yMax, double yMin ){
	// Determina el tamaño adecuado de la unidad para la gráfica
	double uni;
	//uni = ( yMax - yMin ) / static_cast < double > ( divEje );
	// La redondea a la potencia de 10 más cercana
	uni = alto / static_cast < double > ( divEje ) / esc;
	uni = redondea( uni );
	return uni;
}
//............................................................................
double graficador::redondea( const double cval ){
	// Redondea un número a la potencia de 10 más cercana
	int pot10 = 0;
	double aux, val = fabs( cval );
	if( fabs( val ) < 1.0 ){
        while( aux < 1 ){
            val *= 10.0;
            aux  = static_cast < int > ( val );
            pot10--;
        }
        if( val >= 5.0 ) pot10++;
    }
    else{
        while( val > 1.0 ){
            val /= 10.0;
            pot10++;
        }
        val *= 10.0;
        if( val < 5.0 ) pot10--;
		//pot10--;
                              
    }
    return pow( 10.0, pot10 );
}
//...........................................................................
void graficador::dibujaPuntos( int n, double *X, double *Y, pincel &pinDot ){
	// Dibuja una serie de puntos
	double x, y;
	for( int i = 0; i < n; i++ ){
		transforma( X[i], Y[i], areaVista.val[0][0], areaVista.val[1][0],
		            areaVista.val[0][1], areaVista.val[3][1], esc, x, y );
		pinDot.dibujaPunto( x, y, 5.0 );
	}
}
void graficador::dibujaPuntos( int n, double *X, double *Y ){
	// Dibuja una serie de puntos
	double x, y;
	pincel pinDot( sup, 0.9, 0.1, 0.1, 3.0, 1.0 );
	for( int i = 0; i < n; i++ ){
		transforma( X[i], Y[i], areaVista.val[0][0], areaVista.val[1][0],
		            areaVista.val[0][1], areaVista.val[3][1], esc, x, y );
		pinDot.dibujaPunto( x, y, 4.0 );
	}
}
//............................................................................
void graficador::graficaFija( double xMin, double xMax, int nPunt ){
	// Dibuja una función en un intervalo definido por las coordenadas "y"
	// máxima y mínima en un área de dibujo definida
	//double esc, yMin, yMax;
	double aux;
	// Genera puntos
	this->nPunt = nPunt;
	yMax = yMin = 0.0;
	// Pide el tamaño para la matriz de coordenadas
	coord.reinicia( nPunt, 2 );
	generaPuntos( nPunt, xMin, xMax, yMin, yMax );
	// Calcula la escala
	esc = escala( xMin, xMax, yMin, yMax );
	// Límites
	// Del área visible
	aux = 0.5 * ( largo / esc - ( xMax - xMin ) );
	areaVista.val[0][0] = xMin - aux;
	areaVista.val[0][1] = 0.5 * ( yMax + yMin - alto / esc );
	areaVista.val[1][0] = xMax + aux;
	areaVista.val[1][1] = 0.5 * ( yMax + yMin - alto / esc );
	areaVista.val[2][0] = xMax + aux;
	areaVista.val[2][1] = 0.5 * ( yMax + yMin + alto / esc );	
	areaVista.val[3][0] = xMin - aux;
	areaVista.val[3][1] = 0.5 * ( yMax + yMin + alto / esc );
	// Del área con parte de la función (en un incio coinciden)
	areaGraf.val[0][0] = xMin;
	areaGraf.val[0][1] = yMin;
	areaGraf.val[1][0] = xMax;
	areaGraf.val[1][1] = yMin;
	areaGraf.val[2][0] = xMax;
	areaGraf.val[2][1] = yMax;	
	areaGraf.val[3][0] = xMin;
	areaGraf.val[3][1] = yMax;
	// Dibuja la gráfica
	dibujaGraf( nPunt, coord.val, areaVista.val, areaGraf.val, esc );
}
//............................................................................
double graficador::escala( double xMin, double xMax, double yMin,
                           double yMax ){
	// Determina la escala mayor
	double fA, fL, esc;
	fA = alto / fabs( yMax - yMin );
	fL = largo / fabs( xMax - xMin );
	esc = fA;
	tipoEsc = 0;
	if( esc > fL ){
		esc = fL;
		tipoEsc = 1;
	}
	/*cout << "yMax = " << yMax << endl;
	cout << "yMin = " << yMin << endl;
	cout << "xMax = " << xMax << endl;
	cout << "xMin = " << xMin << endl << endl;*/
	return esc;
}
//............................................................................
void graficador::desplazaFija( int dir, double desp ){
	// Desplaza la gráfica
	//double yMin, yMax, esc;
	// Cambio del área que se puede ver y el área visible de la gráfica
	// Desplazamiento vertical (se mueven las dos áreas igual (no se debería
	// hace por los posibles errores de desbordamiento))
	// Desplazamineto vertical
	if( dir == 0 ){
		// Del área visible
		areaVista.val[0][1] += desp;
		areaVista.val[1][1] += desp;
		areaVista.val[2][1] += desp;	
		areaVista.val[3][1] += desp;
	}
	// Desplazamiento horizontal
	else{
		// Del área visible
		areaVista.val[0][0] += desp;
		areaVista.val[1][0] += desp;
		areaVista.val[2][0] += desp;
		areaVista.val[3][0] += desp;
	}
	// Dibuja la gráfica
	dibujaGraf( nPunt, coord.val, areaVista.val, areaGraf.val, esc );	
}
//............................................................................
void graficador::zoomFija( double f ){
	double aux1;
	// Calcula la escala
	esc = escala( areaGraf.val[3][1] + f, areaGraf.val[0][1] - f );
	aux1 = f * largo / alto;
	// Límites
	// Del área visible
	areaVista.val[0][0] -= aux1;
	areaVista.val[0][1] -= f;
	areaVista.val[1][0] += aux1;
	areaVista.val[1][1] -= f;
	areaVista.val[2][0] += aux1;
	areaVista.val[2][1] += f;
	areaVista.val[3][0] -= aux1;
	areaVista.val[3][1] += f;
	// Dibuja la gráfica
	dibujaGraf( nPunt, coord.val, areaVista.val, areaGraf.val, esc );
}
//============================================================================
