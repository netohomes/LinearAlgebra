#include "metInterpolacion.h"
#include "linear_solvers.h"
//#include "../extra/matrix.h"
#define ANCHO_IMP 10
using namespace std;
pol_Newton Newton; 
splines_C0 C0;
splines_C1 C1;
splines_C2 C2;
polinomios P;
pol_Lagrange Lag;
minCuad MC;
//============================================================================
//                            INTERPOLACIÓN 2D
//============================================================================
interpolacion::interpolacion( void ){}
interpolacion::~interpolacion( void ){ }
//............................................................................
void interpolacion::interpola( const int metodo, const int nCoord, double *X,
                               double *Y, int argc, char **argv,
                               int calcError, int nE, double *datE ){
	// Hace la interpolacion de los puntos con el método escogido
	// Selección del método
	//vector operVec;
	grafFunc GG;
	double p, q, error;
	double (*func)(double);
	double (*Ferror)(double);
	Ferror = F;
	switch( metodo ){
		// Con polinomios
		case 1:
			std::cout << "Interpolando con polinomios" << std::endl;
			P.iniciar( nCoord, X, Y );
			P.calc_coef( );
			func = evalFunc_pol;
			break;
		case 2:
			std::cout << "Interpolando con polinomios de Lagrange";
			std::cout << std::endl;
			Lag.iniciar( nCoord, X, Y );
			Lag.calc_coef( );
			func = evalFunc_polLagrange;
		//operVec.print( nCoord, a, 1, "Coeficientes del polinimio (a_i)" );
			break;
		case 3:
			std::cout << "Interpolando con splines C0";
			std::cout << std::endl;
			C0.iniciar( nCoord, X, Y );
			C0.ordenDat( );
			func = evalFunc_splinesC0;
			break;
		case 4:
			std::cout << "Interpolando con splines C1";
			std::cout << std::endl;
			C1.iniciar( nCoord, X, Y );
			C1.ordenDat( );
			C1.calc_zi( );
			func = evalFunc_splinesC1;
			break;
		case 5:
			std::cout << "Interpolando con splines C2";
			std::cout << std::endl;
			C2.iniciar( nCoord, X, Y );
			C2.ordenDat( );
			C2.calc_coef( );
			func = evalFunc_splinesC2;
			break;
		case 6:
			std::cout << "Interpolando con polinomios de Newton";
			std::cout << std::endl;
			Newton.coef_Newton( nCoord, X, Y );
			func = evalFunc_Newton;
			break;
		case 7:
			std::cout << "Regresion con Minimos Cuadrados";
			std::cout << std::endl;
			MC.iniciar( nCoord, X, Y );
			MC.calc_coef( );
			func = evalFunc_minCuad;
			break;
		default:
			std::cout << "ERROR: metodo no valido" << std::endl;
	}
	// Calcula el error
	if( calcError == 1 ){
		error = mideError( nE, datE, Ferror, func );
		cout << endl << "Error = " << error << endl;
	}
	// Grafica los resultados
	p = operVec.min( nCoord, X ); // mínimo de "x"
	q = operVec.max( nCoord, X ); // máximo de "x"
	GG.grafica( p, q, func, 1000, 800, argc, argv, nCoord,  X, Y );
}
//............................................................................
// Datos para interpolacion
//............................................................................
datInterpolacion::datInterpolacion( ): X( ), Y( ), ordenar( ){ }
datInterpolacion::~datInterpolacion( ){ }
//............................................................................
void datInterpolacion::iniciar( int n, double *Xp, double *Yp ){
	X.reinicia( n );
	Y.reinicia( n );
	operVec.copia_val( n , X.val, Xp );
	operVec.copia_val( n , Y.val, Yp );
} 
//............................................................................
//                             Con POLINOMIOS
//............................................................................
polinomios::polinomios( ): datInterpolacion( ){ }
polinomios::~polinomios( ){ } 
//............................................................................
void polinomios::calc_coef( void ){
	// Interpola con polinomios de grado "n" de forma que éste pasa por todos
	// los puntos. Es de la forma f(x) = a1*x + a2*x^2 + ... + an*x^n
	// nPunt = núm. de puntos
	// dim   = coordenadas de cada punto (dimensión)
	// {X}   = coordenadas en "x" de los puntos
	// {Y}   = coordenadas en "y" de los puntos
	// {a}   = para los coeficientes del polinomio
	int nPunt = X.n;
	solver sol;
	matriz<double> A( nPunt, nPunt ); // matriz<double> de coeficientes
	a.reinicia( nPunt );      // Pide memoria
	double x;
	for( int i = 0; i < nPunt; i++ ){
		x = X.val[i];
		A.val[i][0] = 1.0;
		for( int j = 1; j < nPunt; j++ ){
			A.val[i][j] = x;
			x *= X.val[i];
		}
	}
	sol.factoriza_LU( nPunt, A.val );
	sol.resuelve_LU( nPunt, A.val, a.val, Y.val );
}
//............................................................................
double polinomios::eval_func( double x ){
	double sum, pot;
	sum = 0.0;
	pot = 1.0;
	for( int i = 0; i < X.n; i++ ){
		sum += pot * a.val[i];
		pot *= x;
	}
	return sum;
}
//............................................................................
double evalFunc_pol( double x ){
	return P.eval_func( x );
}
//............................................................................
//                         Con POLINOMIOS DE LAGRANGE
//............................................................................
pol_Lagrange::pol_Lagrange( ): datInterpolacion( ){ }
pol_Lagrange::~pol_Lagrange( ){ }
//............................................................................
void pol_Lagrange::calc_coef( void ){
	// Interpola con polinomios de Lagrange
	int nPunt = X.n;
	solver sol;
	matriz<double> A( nPunt, nPunt ); // matriz<double> de coeficientes
	a.reinicia( nPunt );      // Pide memoria
	double L;
	for( int i = 0; i < nPunt; i++ ){
		for( int j = 0; j < nPunt; j++ ){
			L = 1.0;
			for( int k = 0; k < nPunt; k++ ){
				if( k != j )
					L *= ( X.val[i] - X.val[k] ) / ( X.val[j] - X.val[k] );
			}
			A.val[i][j] = L;
		}
	}
	sol.factoriza_LU( nPunt, A.val );
	sol.resuelve_LU( nPunt, A.val, a.val, Y.val );	
}
//............................................................................
double pol_Lagrange::eval_func( double x ){
	double sum, f;
	sum = 0.0;
	for( int i = 0; i < X.n; i++ ){
		f = 1.0;
		for( int j = 0; j < X.n; j++ ){
			if( i != j ) f *= ( x - X.val[j] ) / ( X.val[i] - X.val[j] );
		}
		sum += a.val[i] * f;
	}
	return sum;
}
//............................................................................
double evalFunc_polLagrange( double x ){
	return Lag.eval_func( x );
}
//............................................................................
//                      Con POLINOMIOS DE NEWTON
//............................................................................
pol_Newton::pol_Newton( void ): F( ), X( ), Y( ){ }
pol_Newton::~pol_Newton( void ){ }
//............................................................................
void pol_Newton::coef_Newton( int nPunt, double *Xp, double *Yp ){
	// Calcula los coeficientes para obtener el polinomio que pasa
	// sobre todos los puntos.
	F.reinicia( nPunt );
	X.reinicia( nPunt );
	Y.reinicia( nPunt );
	operVec.copia_val( nPunt, X.val, Xp );
	operVec.copia_val( nPunt, Y.val, Yp );
	// Calcula las diferencias divididas
	F.val[0] = Yp[0];
	calcF( nPunt - 1, 0, nPunt - 1, F.val, X.val, Y.val );

	for( int i = 0; i < nPunt * ANCHO_IMP; i++ ) std::cout << "-";
	cout << endl;
	operVec.print( nPunt, F.val, 2, "Coeficientes del polinimio de Newton",
	               "----------------------------------------" );
	for( int i = 0; i < nPunt * ANCHO_IMP; i++ ) std::cout << "-";
	cout << endl;
}
//...........................................................................
double pol_Newton::calcF( int n, int ini, int fin, double *F, double *X,
                          double *Y ){
	// Calcula F[X0], F[X1,X0], F[X2,X1,X0], ...
	if( fin - ini == 0 ) return Y[ini];
	F[n] = ( calcF( n - 1, ini + 1, fin, F, X, Y ) -
	         calcF( n - 1, ini , fin - 1, F, X, Y ) ) / ( X[fin] - X[ini] );
}
//............................................................................
double pol_Newton::eval_func( double x ){
	double y, p;
	y = 0.0;
	for( int i = 0; i < F.n ; i++ ){
		p = 1.0;
		for( int j = 0; j < i; j++ ) p *= ( x - X.val[j] );
		y += F.val[i] * p;
	}
	return y;
}
//............................................................................
double evalFunc_Newton( double x ){
	return Newton.eval_func( x );
}
//............................................................................
//                               Con SPLINES
//............................................................................
// Continuidad C0
//............................................................................
splines_C0::splines_C0( ): datInterpolacion( ) { }
splines_C0::~splines_C0( ){ }
//............................................................................
void splines_C0::ordenDat( void  ){
	ordenar.mergeSort( 0, X.n - 1, X.val, true, Y.val );
}
//............................................................................
double splines_C0::eval_func( double x ){
	// Identifica en que lugar está
	int k;
	double y, m, tol = 1e-10;
	for( int i = 1; i < X.n; i++ ){
		if( x >= X.val[i-1] - tol && x <= X.val[i] + tol ){
			k = i;
			break;
		}
		if( x > X.val[X.n-1] + tol ) return 0.0;
		if( x < X.val[0] - tol )     return 0.0;
	}
	// Hace una interpolación con rectas
	m = ( Y.val[k] - Y.val[k-1] ) / ( X.val[k] - X.val[k-1] );
	y = Y.val[k-1] + m * ( x - X.val[k-1] );
	//cout << endl << "x = " << x << endl;
	//cout << "y = " << y << endl;
	return y;
}
//...........................................................................
double evalFunc_splinesC0( double x ){
	return C0.eval_func( x );
}
//............................................................................
// Continuidad C1
//............................................................................
splines_C1::splines_C1( ): datInterpolacion( ) { }
splines_C1::~splines_C1( ){ }
//............................................................................
void splines_C1::ordenDat( void  ){
	ordenar.mergeSort( 0, X.n - 1, X.val, true, Y.val );
}
//............................................................................
void splines_C1::calc_zi( void ){
    // Interpola con splines para un "x" dado
    z.reinicia( X.n );
    // Se escoge z0 arbitrariamente
    z.val[0] = ( Y.val[1] - Y.val[0] ) / ( X.val[1] - X.val[0] );
    // Cálculo de las zi
    for( int i = 0; i < X.n - 1; i++ ){
        z.val[i+1] = -z.val[i] + 2.0 * ( Y.val[i+1] - Y.val[i] ) /
                      ( X.val[i+1] - X.val[i] );
    }
}
//............................................................................
double splines_C1::eval_func( double x ){
    // Interpola para un punto dado
    double Q;
    int k;
    // Identifica entre que puntos se encuentra
    if( x > X.val[X.n-1] ) return Y.val[Y.n-1];
	if( x < X.val[0] )     return Y.val[0];
	for( int i = 1; i < X.n; i++ ){
		if( x >= X.val[i-1] && x <= X.val[i] ){
			k = i;
			break;
		}
	}
	// Hace la interpolacion
	Q = ( z.val[k] - z.val[k-1] ) / ( 2.0 *  ( X.val[k] - X.val[k-1] ) ) *
        ( x - X.val[k-1] ) * ( x - X.val[k-1] ) + z.val[k-1] *
	    ( x - X.val[k-1] ) + Y.val[k-1];
    return Q;
}
//...........................................................................
double evalFunc_splinesC1( double x ){
	return C1.eval_func( x );
}
//...........................................................................
// Continuidad C2
//............................................................................
splines_C2::splines_C2( ): A( ), B( ), S( ), t( ){
	C1 = 0.0;
	C2 = 0.0;
}
splines_C2::~splines_C2( ){ }
//............................................................................
void splines_C2::ordenDat( void  ){
	ordenar.mergeSort( 0, X.n - 1, X.val, true, Y.val );
}
//............................................................................
void splines_C2::calc_coef( void ){
	// Resuelve el sistema tridiagonal para obtener los coeficientes de cada
	// spline. Los puntos ya deben estar ordenados
	solver resuelve;
	double h1, h2, t1, t2, h;
	// Pide la memoria de los vectores
	int n = X.n;
	A.reinicia( n - 2 ); // La diagonal principal
	B.reinicia( n - 2 ); // Digonal secundaria
	S.reinicia( n - 2 ); // Segundas derivadas
	t.reinicia( n - 2 ); // Términos independientes
	// Calcula la matriz<double> y el vector de términos independientes
	h1 = X.val[1] - X.val[0];
	t1 = Y.val[1] - Y.val[0];
	for( int i = 1; i < n - 1; i++ ){
		h2 = X.val[i+1] - X.val[i];
		t2 = Y.val[i+1] - Y.val[i];
		B.val[i-1] = h2;
		A.val[i-1] = ( h1 + h2 ) * 2.0; // Sobre la diagonal
		t.val[i-1] = ( t2 / h2 - t1 / h1 ) * 6.0;
		h1 = h2;
		t1 = t2;
	}
	
	// Impone las condciones de contorno
	h = X.val[1] - X.val[0];	
	t.val[0] -= h * C1;
	h = X.val[n-1] - X.val[n-2];  
	t.val[n-3] -= h * C2;
	// Resuelve el sistema
	resuelve.sistema3diagonalSim( n - 2, A.val, B.val, S.val, t.val ); 
	S.print( );

	/*cout << endl;
	operVec.print( S.n, S.val, 2, "Segundas derivadas (Splines C2)",
	               "-------------------------------------------------" );
		for( int i = 0; i < ( X.n - 2 ) * ANCHO_IMP; i++ ) std::cout << "-";*/
	cout << endl;
}
//............................................................................
double splines_C2::eval_func( double x ){
	// Identifica entre que puntos se encuentra
	int k;
	double a, b, c, d, e, f, h, y, S1, S2;
	if( x > X.val[X.n-1] ) return Y.val[Y.n-1];
	if( x < X.val[0] )     return Y.val[0];
	for( int i = 1; i < X.n; i++ ){
		if( x >= X.val[i-1] && x <= X.val[i] ){
			k = i - 1; // OJO
			break;
		}
	}
	// Evalúa la función
	if( k == 0 ){
		S1 = C1;
		S2 = S.val[0];
	}
	else if( k == X.n - 2 ){
		S1 = S.val[S.n-1];
		S2 = C2;
	}
	else{
		S1 = S.val[k-1];
		S2 = S.val[k];
	}
	h = X.val[k+1] - X.val[k];
	a = S1 * ( X.val[k+1] - x ) * ( X.val[k+1] - x ) *
	    ( X.val[k+1] - x ) / ( 6.0 * h );
	b = S2 * ( x - X.val[k] ) * ( x - X.val[k] ) *
	    ( x - X.val[k] ) / ( 6.0 * h );
	c = Y.val[k+1] * ( x - X.val[k] ) / h;
	d = Y.val[k] * ( X.val[k+1] - x ) / h;
	e = S2 * h * ( x - X.val[k] ) / 6.0;
	f = S1 * h * ( X.val[k+1] - x ) / 6.0;	
	y = a + b + c + d - e - f;
	//cout << "x = " << x << "  y = " << y << endl;
	return y; 
}
//............................................................................
double evalFunc_splinesC2( double x ){
	return C2.eval_func( x );
}
//............................................................................
//                    REGRESIÓN CON MÍNIMOS CUADRADOS
//............................................................................
minCuad::minCuad( ): M( ), C( ), beta( ){ }
minCuad::~minCuad( ){ }
//............................................................................
void minCuad::calc_coef( void ){
	// Calcula los parámetros de la regresión por mínimos cuadrados para
	// una función fija
	solver sol;
	// Calcula la matriz<double> [M] y el vector {C}
	M.reinicia( 5, 5 );
	beta.reinicia( 5 );
	C.reinicia( 5 );
	calc_MC( );
	// Obtiene los parámetros
	sol.factoriza_LU( 5, M.val );
	sol.resuelve_LU( 5, M.val, beta.val, C.val );	
	//beta.print( );

	cout << endl;
	operVec.print( 5, beta.val, 2,
	               "Parametros B para la Regresion con Min. Cuadrados",
	               "-------------------------------------------------" );
	for( int i = 0; i < 5 * ANCHO_IMP; i++ ) std::cout << "-";
	cout << endl;

}
//............................................................................
void minCuad::calc_MC( void ){
	// Calcula la matriz<double> para la regresión
	double x, y;
	C.ini( 0.0 );
	M.iniVal( 0.0 );
	for( int i = 0; i < X.n; i++ ){
		x = X.val[i];
		y = Y.val[i];
		// matriz<double> M 
		M.val[0][0] += sin( x * x ) * sin( x * x );
		M.val[0][1] += cos( x )     * sin( x * x );
		M.val[0][2] += exp( -x )    * sin( x * x );
		M.val[0][3] += log( x )     * sin( x * x );
		M.val[0][4] += x            * sin( x * x );

		M.val[1][0] += sin( x * x ) * cos( x );
		M.val[1][1] += cos( x )     * cos( x );
		M.val[1][2] += exp( -x )    * cos( x );
		M.val[1][3] += log( x )     * cos( x );
		M.val[1][4] += x            * cos( x );

		M.val[2][0] += sin( x * x ) * exp( -x );
		M.val[2][1] += cos( x )     * exp( -x );
		M.val[2][2] += exp( -x )    * exp( -x );
		M.val[2][3] += log( x )     * exp( -x );
		M.val[2][4] += x            * exp( -x );

		M.val[3][0] += sin( x * x ) * log( x );
		M.val[3][1] += cos( x )     * log( x );
		M.val[3][2] += exp( -x )    * log( x );
		M.val[3][3] += log( x )     * log( x );
		M.val[3][4] += x            * log( x ); 

		M.val[4][0] += sin( x * x ) * x;
		M.val[4][1] += cos( x )     * x;
		M.val[4][2] += exp( -x )    * x;
		M.val[4][3] += log( x )     * x; 
		M.val[4][4] += x            * x;
		
		// Vector C
		C.val[0] += y * sin( x * x );
		C.val[1] += y * cos( x );
		C.val[2] += y * exp( -x );
		C.val[3] += y * log( x );
		C.val[4] += y * x;
	}
}
//............................................................................
double minCuad::eval_func( double x ){
	double y;
	y = beta.val[0] * sin( x * x ) + beta.val[1] * cos( x ) +
	    beta.val[2] * exp( -x ) + beta.val[3] * log( x ) + beta.val[4] * x;
	return y;
}
//............................................................................
double evalFunc_minCuad( double x ){
	return MC.eval_func( x );
}
//............................................................................
//............................................................................
double mideError( int n, double *X, double (*F)( double ),
                  double (*I)( double ) ){
	// Obtiene el error entre una función F y una función de interpolación I
	// para un conjunto de x's
	double x, sum = 0.0;
	for( int i = 0; i < n; i++ ){
		x = X[i];
		sum += ( F( x ) - I( x ) ) * ( F( x ) - I( x ) );
	}
	return sqrt( sum );
}
double F( double x ){
	return x + x * sin( 0.5 * x ) / 3.0;
}
//............................................................................
//............................................................................
//============================================================================
