#include "derivatives_and_integrals.h"
#include "matrix.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include "linear_solvers.h"
using namespace std;
//==================================================================
//                         INTEGRACIÓN
//==================================================================
int_deriv::int_deriv(  ){ }
int_deriv::~int_deriv( ){ }
//..................................................................
double int_deriv::simpson( double a, double b, int n, double *x ){
	// Integra sobre una variable usando la regla de simpson para
	// una función tipo f(x,y,z)
	double sum1, sum2, I;
	sum1 = 0.0;
	for( int k = 1; k < n - 1; k+=2 ) sum1 += x[k];
	sum2 = 0.0;
	for( int k = 2; k < n - 2; k+=2 ) sum2 += x[k];
	I = ( b - a ) * ( x[0] + 4.0 * sum1 + 2.0 * sum2 + x[n-1] ) /
	    ( 3.0 * static_cast < double > ( n - 1 ) ); // OJO!!! se
	//  tiene que dividir por el número de divisiones del rango,
	// no entre el número de nodos que se necesitan para la división
	return I;
}
//..................................................................
double int_deriv::simpson3( double (*f)( double, double, double ),
                            double x1, double x2, double y1, double y2,
                            double z1, double z2, int n ){
	// Aplica la regla se Simpson 1/3 para integrar una función dada
	double dx, dy, dz, I, x, y, z;
	//int n = 5;
	miVector < double > v( n );
	miVector < double > Ix( n );
	miVector < double > Iy( n );
	// Deltas en cada dirección
	dx = ( x2 - x1 ) / static_cast < double > ( n - 1 );
	dy = ( y2 - y1 ) / static_cast < double > ( n - 1 );
	dz = ( z2 - z1 ) / static_cast < double > ( n - 1 );
	// Hace la integral triple
	z = z1;
	for( int i = 0; i < n; i++ ){           // Para Z
		y = y1;
		for( int j = 0; j < n; j++ ){       // Para Y
			x = x1;                         // Para X
			for( int k = 0; k < n; k++ ){
				v.val[k] = f( x, y, z );
				x += dx;
			}
			Ix.val[j] = simpson( x1, x2, n, v.val );
			y += dy;
		}
		Iy.val[i] = simpson( y1, y2, n, Ix.val );
		z += dz;
	}
	I = simpson( z1, z2, n, Iy.val );
	return I;
}
//..................................................................
double int_deriv::monteCarlo3( double (*f)( double, double, double ),
                               double x1, double x2, double y1,
                               double y2, double z1, double z2, int n ){
	// Aplica el método de Monte Carlo a una integral triple
	double x, y, z, sum, I;
	//int n = 1000;
	sum = 0.0;
	for( int i = 0; i < n; i++ ){
		x = myRandom( x1, x2 );
		y = myRandom( y1, y2 );
		z = myRandom( z1, z2 );
		sum += f( x, y, z );
	}
	I = ( z2 - z1 ) * ( y2 - y1 ) * ( x2 - x1 ) * sum /
	    static_cast < double > ( n );
	return I;
}
//..................................................................
double int_deriv::cuadGauss3( double (*f)( double, double, double ),
                            double x1, double x2, double y1, double y2,
                            double z1, double z2, int n ){
	// Aplica el método de integación por Gauss para integrales triples
	// usando hexaedros (cubos)
	// n -> Puntos de integración
	double sum, x, y, z, t;
	sum = 0.0;
	for( int i = 0; i < n; i++ ){
		t = ptosGauss( n, i+1 );
		z = ( ( z2 - z1 ) * t + z2 + z1 ) / 2.0;
		for( int j = 0; j < n; j++ ){
			t = ptosGauss( n, j+1 );
			y = ( ( y2 - y1 ) * t + y2 + y1 ) / 2.0;
			for( int k = 0; k < n; k++ ){
				t = ptosGauss( n, k+1 );
				x = ( ( x2 - x1 ) * t + x2 + x1 ) / 2.0;
				sum += pesosGauss( n, k+1 ) * pesosGauss( n, j+1 ) *
				       pesosGauss( n, i+1 ) * f( x, y, z );
			}
		}
	}
	sum *= ( x2 - x1 ) * ( y2 - y1 ) * ( z2 - z1 ) / 8.0;
	return sum;
}
//..................................................................
double int_deriv::myRandom( double a, double b ){
	double r;
	r = static_cast < double > ( rand( ) ) /
	    static_cast < double > (RAND_MAX - 1) * ( b - a ) + a;
	return r;
}
//..................................................................
double int_deriv::ptosGauss( const int n, const int i ){
	// Determina los puntos de Gauss que se van a emplear
	// n -> Número total de puntos de Gauss que se emplearán para
	// la integral i -> Punto "i" de integración
	// CASI TODOS LOS PUNTOS IMPARES ESTAN MAL ESCRITOS. HAY QUE
	// CORREGIR LOS PUNTOS Y LOS PESOS
	double x;
	switch( n ){
		case 1:
			x = 0.0;
			break;
		case 2:
			if( i == 1 ) x =  1.0 / sqrt( 3.0 );
			if( i == 2 ) x = -1.0 / sqrt( 3.0 );
			break;
		case 3:
			// Antes estaban mal
			if( i == 1 ) x =   sqrt( 3.0 / 5.0 );
			if( i == 2 ) x =   0.0;
			if( i == 3 ) x =  -sqrt( 3.0 / 5.0 );
			break;
		case 4:
			if( i == 1 ) x =  0.861136311594053;
			if( i == 2 ) x =  0.339981043584856;
			if( i == 3 ) x = -0.339981043584856;
			if( i == 4 ) x = -0.861136311594053;
			break;
		case 5:
			if( i == 1 ) x =  0.9061798459386640;
			if( i == 2 ) x =  0.5384693101056831;
			if( i == 3 ) x =  0.0;
			if( i == 4 ) x = -0.5384693101056831;
			if( i == 5 ) x = -0.9061798459386640;
		case 6:
			if( i == 1 ) x =  0.9324695142031521;
			if( i == 2 ) x =  0.6612093864662645;
			if( i == 3 ) x =  0.2386191860831969;
			if( i == 4 ) x = -0.2386191860831969;
			if( i == 5 ) x = -0.6612093864662645;
			if( i == 6 ) x = -0.9324695142031521;
			break;
		case 7:
			if( i == 1 ) x =  0.9491079123427585;
			if( i == 2 ) x =  0.7415311855993945;
			if( i == 3 ) x =  0.4058451513773972;
			if( i == 4 ) x =  0.0;
			if( i == 5 ) x = -0.4058451513773972;
			if( i == 6 ) x = -0.7415311855993945;
			if( i == 7 ) x = -0.9491079123427585;
			break;
		case 8:
			if( i == 1 ) x =  0.9602898564975363;
			if( i == 2 ) x =  0.7966664774136267;
			if( i == 3 ) x =  0.5255324099163290;
			if( i == 4 ) x =  0.1834346424956498;
			if( i == 5 ) x = -0.1834346424956498;
			if( i == 6 ) x = -0.5255324099163290;
			if( i == 7 ) x = -0.7966664774136267;
			if( i == 8 ) x = -0.9602898564975363;
			break;
        case 9:
			if( i == 1 ) x =  0.9681602395076261;
			if( i == 2 ) x =  0.8360311073266358;
			if( i == 3 ) x =  0.6133714327005904;
			if( i == 4 ) x =  0.3242534234038089;
			if( i == 5 ) x =  0.0;
			if( i == 6 ) x = -0.3242534234038089;
			if( i == 7 ) x = -0.6133714327005904;
			if( i == 8 ) x = -0.8360311073266358;
			if( i == 9 ) x = -0.9681602395076261;
			break;
		default:
			cout << "ERROR: num. de ptos. de int. no valido" << endl;
	}
	return x;
}
//..................................................................
double int_deriv::pesosGauss( const int n, const int i ){
	// Determina los puntos de Gauss que se van a emplear
	// n -> Número total de puntos de Gauss que se emplearán para
	// ñla integral i -> Peso "i" de integración
	double p;
	switch( n ){
		case 1:
			p = 2.0;
			break;
		case 2:
			if( i == 1 ) p = 1.0;
			if( i == 2 ) p = 1.0;
			break;
		case 3:
			if( i == 1 ) p =  5.0 / 9.0;
			if( i == 2 ) p =  8.0 / 9.0;
			if( i == 3 ) p =  5.0 / 9.0;
			break;
		case 4:
			if( i == 1 ) p = 0.3478548451374538;
			if( i == 2 ) p = 0.6521451548625461;
			if( i == 3 ) p = 0.6521451548625461;
			if( i == 4 ) p = 0.3478548451374538;
			break;
		case 5:
			if( i == 1 ) p = 0.2369268850561891;
			if( i == 2 ) p = 0.4786286704993665;
			if( i == 3 ) p = 0.5688888888888889;
			if( i == 4 ) p = 0.4786286704993665;
			if( i == 5 ) p = 0.2369268850561891;
		case 6:
			if( i == 1 ) p = 0.1713244923791704;
			if( i == 2 ) p = 0.3607615730481386;
			if( i == 3 ) p = 0.4679139345726910;
			if( i == 4 ) p = 0.4679139345726910;
			if( i == 5 ) p = 0.3607615730481386;
			if( i == 6 ) p = 0.1713244923791704;
			break;
		case 7:
			if( i == 1 ) p = 0.1294849661688697;
			if( i == 2 ) p = 0.2797053914892766;
			if( i == 3 ) p = 0.3818300505051189;
			if( i == 4 ) p = 0.4179591836734694;
			if( i == 1 ) p = 0.3818300505051189;
			if( i == 2 ) p = 0.2797053914892766;
			if( i == 3 ) p = 0.1294849661688697;
			break;
		case 8:
			if( i == 1 ) p = 0.1012285362903763;
			if( i == 2 ) p = 0.2223810344533745;
			if( i == 3 ) p = 0.3137066458778873;
			if( i == 4 ) p = 0.3626837833783620;
			if( i == 5 ) p = 0.3626837833783620;
			if( i == 6 ) p = 0.3137066458778873;
			if( i == 7 ) p = 0.2223810344533745;
			if( i == 8 ) p = 0.1012285362903763;
			break;
        case 9:
            if( i == 1 ) p = 0.0812743883615744;
            if( i == 2 ) p = 0.1806481606948574;
            if( i == 3 ) p = 0.2606106964029354;
            if( i == 4 ) p = 0.3123470770400029;
            if( i == 5 ) p = 0.3302393550012598;
            if( i == 6 ) p = 0.3123470770400029;
            if( i == 7 ) p = 0.2606106964029354;
            if( i == 8 ) p = 0.1806481606948574;
            if( i == 9 ) p = 0.0812743883615744;
            break;
		default:
			cout << "ERROR: num. de ptos. de int. no valido" << endl;
	}
	return p;
}
//..................................................................
double int_deriv::trapecio(double (*f)( double, double, double ), 
                           double x1, double x2, double y1, double y2, 
                           double z1, double z2 ) {
	// Aplica la regla del trapecio para hacer la integración triple
	// de una función f(x,y,z)
	double dx, dy, dz, I, x, y, z;
	int nDiv = 500;
	// Deltas en cada dirección
	dx = ( x2 - x1 ) / static_cast < double > ( nDiv );
	dy = ( y2 - y1 ) / static_cast < double > ( nDiv );
	dz = ( z2 - z1 ) / static_cast < double > ( nDiv );
	// Hace la integral triple
	I = 0.0;
	z = z1;
	for( int i = 0; i < nDiv; i++ ){           // Para Z
		y = y1;
		for( int j = 0; j < nDiv; j++ ){       // Para Y
			x = x1;
			for( int k = 0; k < nDiv; k++ ){   // Para X
				if( k == 0 || k == nDiv - 1 )
					I += 0.5 * f(x,y,z) * dx * dy * dz;
				else                         
					I += f(x,y,z) * dx * dy * dz;
				x += dx;
			}
			y += dy;
		}
		z += dz;
	}
	return I;
}
//==================================================================
//                         DERIVACIÓN
//==================================================================
double int_deriv::difAdelante( double (*f)( double ), double x, int d,
                            int p, double dx ){
    // Obtiene la n-ésima con diferencias hacia adelante
    // f(x) = función a derivar
    // x    = punto en el cual se desea evaluar la derivada
    // d    = orden de la derivada que se desea calcular
    // p    = orden del error
    int n;
    double sum, z, fact, pot;
    solver res;
    n = d + p;
    matriz <double> A( n, n );
    miVector <double> b( n );
    miVector <double> C( n );
    for( int j = 0; j < n; j++ ) A.val[0][j] = 1.0;
    for( int j = 0; j < n; j++ ) A.val[1][j] = static_cast <double> (j);
    for( int i = 2; i < n; i++ ){
        for( int j = 0; j < n; j++ )
            A.val[i][j] = A.val[1][j] * A.val[i-1][j];
    }
    //A.print( );
    b.ini( 0.0 );
    b.val[ d ] = 1.0;
    //b.print();
    res.factoriza_LU( n, A.val );
    res.resuelve_LU( n, A.val, C.val, b.val );
    fact = 1.0;
    for( int i = 1; i <= d; i++ ) fact *= static_cast <double> (i);
    for( int i = 0; i < n ; i++ ) C.val[i] *= fact;
    //C.print();
    sum = 0.0;
    z = x;
    for( int i = 0; i < n; i++ ){
        sum += C.val[i] * f( z );
        z += dx;
    }
    pot = pow( dx, d );
    sum /= pot;
    return sum;
}
//..................................................................
double int_deriv::difCentrada( double (*f)( double ), double x, int d,
                               int p, double dx ){
    // Obtiene la n-ésima con diferencias hacia adelante
    // f(x) = función a derivar
    // x    = punto en el cual se desea evaluar la derivada
    // d    = orden de la derivada que se desea calcular
    // p    = orden del error
    int n, m;
    double sum, z, fact, pot;
    solver res;
    m = ( ( d + p - 1 ) - ( d + p - 1 )%2 ) / 2;
    n =  2 * m + 1;
    matriz <double> A( n, n );
    miVector <double> b( n );
    miVector <double> C( n );
    for( int j = 0; j < n; j++ ) A.val[0][j] = 1.0;
    for( int j = 0; j < n; j++ ) A.val[1][j] = static_cast <double> (j-m);
    for( int i = 2; i < n; i++ ){
        for( int j = 0; j < n; j++ )
            A.val[i][j] = A.val[1][j] * A.val[i-1][j];
    }
    A.print( );
    b.ini( 0.0 );
    b.val[ d ] = 1.0;
    b.print();
    res.factoriza_LU( n, A.val );
    res.resuelve_LU( n, A.val, C.val, b.val );
    fact = 1.0;
    for( int i = 1; i <= d; i++ ) fact *= static_cast <double> (i);
    C.print();
    for( int i = 0; i < n; i++ ) C.val[i] *= fact;
    C.print();
    sum = 0.0;
    z = x - static_cast <double> (m) * dx;
    for( int i = 0; i < n; i++ ){
        sum += C.val[i] * f( z );
        z += dx;
    }
    pot = pow( dx, d );
    sum /= pot;
    return sum;
}
//==================================================================
