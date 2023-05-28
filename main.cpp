#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "derivatives_and_integrals.h"
#include "metInterpolacion.h"
#include "datos.h"
#include "eigen_val_vec.h"

using namespace std;

//============================================================================
// Algunas funciones
double esfera( double x, double y, double z ){
	return x*x + y*y + z*z;
}
double cono( double x, double y, double z ){
	return x*x + y*y - z*z;
}
double toro( double x, double y, double z ){
	return ( 1.0 - sqrt( x*x + y*y ) ) * ( 1.0 - sqrt( x*x + y*y ) ) + z * z;
}
double length( double x, double y, double z){
	return sqrt(x*x + y*y + z*z);
}
double f1( double x ){
    return sin( x );
}
double f2( double x ){
    return pow(x,9) + pow(x,8) + pow(x,7) + pow(x,6) + pow(x,5) + pow(x,4) +
           pow(x,3) + pow(x,2) + x;
}
double f3( double x ){
    return sin( x ) * cos( x );
}
//============================================================================
// Ejemplos de derivadas
void derivadas() {
    double (*fd)( double ), deriv;
    int i, d, p;
    float x, dx;
    int_deriv D;

    cout << "Num. de funcion = ";
    scanf( "%d", &i );
    cout << "Orden de la derivada (d) = ";
    scanf( "%d", &d );
    cout << "Orden del error (p) = ";
    scanf( "%d", &p);
    cout << "Punto de evaluacion (x) = ";
    scanf( "%f", &x );
    cout << "Epsilon (h) = ";
    scanf( "%f", &dx );

    if( i == 1 ) fd = f1;
    if( i == 2 ) fd = f2;
    if( i == 3 ) fd = f3;
    deriv = D.difAdelante( fd, x, d, p, dx );
    cout << "-------------------" << endl;
    if( i == 1 ) cout << "Funcion 1" << endl;
    if( i == 2 ) cout << "Funcion 2" << endl;
    if( i == 3 ) cout << "Funcion 3" << endl;
    cout << "-------------------" << endl;
    cout << "d = " << d << " p = " << p << endl;
    cout << "x = " << x << endl;
    cout << "f'(x) = " << deriv << endl;
}
//============================================================================
// Ejemplos de integrales
void integrales() {
    double (*f)( double, double, double );
	double Ig, Is, Im;
	// Integración sobre la esfera
	f = esfera;
	int_deriv I;
	cout << "-------------------------------------" << endl;
	cout << "Integracion sobre la esfera:" << endl;
	cout << "-------------------------------------" << endl;
	Ig = I.cuadGauss3(  f, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 2 );
	Is = I.simpson3(    f, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 3 );
	Im = I.monteCarlo3( f, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1000000 );
	cout << "I (Gauss)       = " << Ig << endl;
	cout << "I (Simpson 1/3) = " << Is << endl;
	cout << "I (Monte Carlo) = " << Im << endl;
	// Integración sobre el cono
	f = cono;
	cout << "-------------------------------------" << endl;
	cout << "Integracion sobre el cono:" << endl;
	cout << "-------------------------------------" << endl;
	Ig = I.cuadGauss3(  f, -1.0, 1.0, -1.0, 1.0, 0.0, 1.0, 2 );
	Is = I.simpson3(    f, -1.0, 1.0, -1.0, 1.0, 0.0, 1.0, 3 );
	Im = I.monteCarlo3( f, -1.0, 1.0, -1.0, 1.0, 0.0, 1.0, 1000000 );
	cout << "I (Gauss)       = " << Ig << endl;
	cout << "I (Simpson 1/3) = " << Is << endl;
	cout << "I (Monte Carlo) = " << Im << endl;
	// Integración sobre el cono
	f = toro;
	cout << "-------------------------------------" << endl;
	cout << "Integracion sobre el toro:" << endl;
	cout << "-------------------------------------" << endl;
	Ig = I.cuadGauss3(  f, -1.0, 1.0, -1.0, 1.0, 0.5, 1.0, 8 );
	Is = I.simpson3(    f, -1.0, 1.0, -1.0, 1.0, 0.5, 1.0, 81 );
	Im = I.monteCarlo3( f, -1.0, 1.0, -1.0, 1.0, 0.5, 1.0, 1000000 );
	cout << "I (Gauss)       = " << Ig << endl;
	cout << "I (Simpson 1/3) = " << Is << endl;
	cout << "I (Monte Carlo) = " << Im << endl;
	// Integración de longitud
	f = length;
	cout << "-------------------------------------" << endl;
	cout << "Integracion de longitud:" << endl;
	cout << "-------------------------------------" << endl;
	Ig = I.cuadGauss3(  f, -1.0, 1.0, -1.0, 1.0, 0.5, 1.0, 3 );
	Is = I.simpson3(    f, -1.0, 1.0, -1.0, 1.0, 0.5, 1.0, 50 );
	Im = I.monteCarlo3( f, -1.0, 1.0, -1.0, 1.0, 0.5, 1.0, 100000000 );
	cout << "I (Gauss)       = " << Ig << endl;
	cout << "I (Simpson 1/3) = " << Is << endl;
	cout << "I (Monte Carlo) = " << Im << endl;	
}
//============================================================================
void interpolacionEjemplos(int argc, char **argv) {
	int metodo, calcError;
	std::string coord, datError;
	interpolacion ajustar;
	miVector<double> X, Y, xError;
	datos dat;
	// Método para la interpolación
	metodo = atoi( argv[1] );
		/*
		1 = polinomios
		2 = polinomios de Lagrange
		3 = Splines C0
		4 = Splines C1
		5 = Splines C2
		6 = Polinomios de Newton
        7 = Minimos cuadrados
		*/
	// Archivo con las coordenadas
	coord = argv[2];
	// Se quiere medir el error?
	calcError = atoi( argv[3] );
	if( calcError == 1 ){
		datError = argv[4];
		dat.leeVector( datError, xError );
	}
	// Lectura de los puntos
	dat.leeCoord2D( coord, X, Y );
	// Interpola
	if( calcError == 1 )
		ajustar.interpola( metodo, X.n, X.val, Y.val, argc, argv,
	                       1, xError.n, xError.val );
	else
		ajustar.interpola( metodo, X.n, X.val, Y.val, argc, argv,
	                       0, 0, NULL );
}
//============================================================================
// Eigenvalores y eigenvectores
void eingenvalores(int argc, char **argv) {
    eigen_vec_val P1;
	datos dat;
	matriz<double> A, PHI;
	// MÈtodo de soluciÛn
	int metodo = atoi( argv[1] ); // 1->Potencia, 2->Potencia inversa
	                              // 3->Cociente de Rayleigh, 4->Jacobi
	int nEigen;
	string archivo_A, archivo_V0, archivo_sigma;
	miVector<double> V0, sigma;
	// Cantidad de valores y vectores propios que se desean calcular
	// Lee matriz
	archivo_A += argv[2];
	dat.leeMatriz( archivo_A, A );
	// Corre el mÈtodo elegido
	switch( metodo ){
		// MÈtodo de la potencia
		case 1:
			nEigen = atoi( argv[3] );
			if( nEigen > A.m ) nEigen = A.m; 
			P1.potencia( A.m, A.val, nEigen, 1, false );
			break;
		// MÈtodo de la potencia inversa 
		case 2:
			nEigen = atoi( argv[3] );
			if( nEigen > A.m ) nEigen = A.m;
			P1.potencia_inversa( A.m, A.val, nEigen, 1, false );
			break;
		// MÈtodo Cociente de Rayleigh
		case 3:
			archivo_V0 += argv[3];
			archivo_sigma += argv[4];
			dat.leeVector( archivo_V0, V0 ); 
			dat.leeVector( archivo_sigma, sigma );
			P1.cociente_Rayleigh2( A.m, A.val, sigma.val[0], V0.val );
			break;
		// MÈtodo de Jacobi 
		case 4:
			PHI.reinicia( A.m, A.n ); // Debe ser cuadrada
			P1.jacobi2( A.m, A.val, PHI.val );
			break;
		// MÈtodo de IteraciÛn en Subespacio
		case 5:
			nEigen = atoi( argv[3] );
			if( nEigen > A.m ) nEigen = A.m; 
			if( nEigen < 2 ) nEigen = 2; 
			PHI.reinicia( A.m, nEigen );
			P1.iterSubEspacio2( A.m, nEigen, A.val, PHI.val );
			break;

		default: 
			std::cout << "ERROR: el metodo que elegiste no es valido";
			std::cout << std::endl;
	}
}
//============================================================================
int main(int argc, char **argv) {
    //integrales();
    //derivadas();
    //interpolacionEjemplos(argc, argv);
    eingenvalores(argc, argv);
    return 0;
}