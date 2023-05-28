#include "linear_solvers.h"
#include <fstream>
#include <iomanip>
#include <cmath>
//#include <omp.h> // OpenMP library
#include "matrix.h"
static int ANCHO_IMP = 13;
static int LONG_NUM  = 12;
using namespace std;
//============================================================================
// MÉTODOS PARA RESOLVER SISTEMAS DE ECUACIONES LINEALES
//============================================================================
//............................................................................
//                          Eliminación Gaussiana
//............................................................................
solver::solver( ){ }
void solver::elimGaussiana( int n, double **A, double *x, double *b ){
    // Resuelve un sistema de n x n ecuaciones por eliminación gaussiana
    // Variables de entrada:
   		// n    Tamaño de la matriz A
    	// A    Matriz de tam. n x n
    	// x    Vector que guardará la solución
	// b    Vector de términos independientes
    // Hace la matriz triangular superior
    double factor;
    for( int j = 0; j < n; j++){
        for( int i = j + 1; i < n; i++){
            factor = -A[i][j] / A[j][j];
            b[i]  += factor * b[j];
            for( int k = 0; k < n; k++ ) A[i][k] += factor * A[j][k];
        }
    }
    resuelveSistTriangSUP( n, A, x, b );
	//cout << endl << "ELIMINACION GAUSSIANA:" << endl;
	//imprime( n, A, x, b, false, true, false );
}
//............................................................................
void solver::elimGaussianaTriangInf( int n, double **A, double *x, double *b ){
    // Resuelve un sistema de n x n ecuaciones por eliminación gaussiana
	// haciendo la matriz A triangular inferior
    // Variables de entrada:
   		// n    Tamaño de la matriz A
    	// A    Matriz de tam. n x n
    	// x    Vector que guardará la solución
	// b    Vector de términos independientes
    // Hace la matriz triangular superior
    double factor;
    for( int j = n - 1; j >= 0; j--){
        for( int i =  j - 1; i >= 0; i--){
            factor = -A[i][j] / A[j][j];
            b[i]  += factor * b[j];
            for( int k = n - 1; k >= 0; k-- ) A[i][k] += factor * A[j][k];
        }
    }
    resuelveSistTriangINF( n, A, x, b );
	cout << endl << "ELIMINACION GAUSSIANA (MAT. TRIANG. SUP):" << endl;
	imprime( n, A, x, b, false, true, false );
}
//............................................................................
void solver::elimGaussianaPivoteo( int n, double **A, double *x, double *b ){
    // Resuelve un sistema de n x n ecuaciones por eliminación gaussiana aplicando
    // pivoteo. Como pivote se emplea el valor máximo (en valor absoluto) que queda
    // en la submatriz. Para colocar el pivote en la posición requerida se
    // intercambian renglones y columnas tomando en cuenta que en cada uno de estos
    // intercambios el determinante de la matriz cambia de signo.
    // Variables de entrada:
        // n    Tamaño de la matriz A
        // A    Matriz de tam. n x n
        // x    Vector que guardará la solución
        // b    Vector de términos independientes
    // Variables locales
        // p, q    Posición del pivote
        // factor  Empleado para hacer los valores debajo del pivote 0's
    int p, q, *pos, intercambios = 0;
    double factor;
    pos = new int [n];
    // Posiciones inciales
    for( int i = 0; i < n; i++ ) pos[i] = i;
    // Hace la matriz triangular superior
    for( int j = 0; j < n; j++){
        // Determina la posición del valor máximo absoluto (pivote)
        valMax( j, n - 1, j, n - 1, A, p, q );
        // Mueve el pivote a la posición adecuada
        intercambios += mueveElemento( j, j, p, q, n, n, A, b, pos );
        for( int i = j + 1; i < n; i++ ){
            factor = -A[i][j] / A[j][j];
            b[i]  += factor * b[j];
            for( int k = 0; k < n; k++ ) A[i][k] += factor * A[j][k];
        }
    }
    resuelveSistTriangSUP( n, A, x, b );
	// Reacomoda el vector solución
	reacomodaX( n, x, pos );
    //cout << endl << "ELIMINACION GAUSSIANA CON PIVOTEO:" << endl;
    //imprime( n, A, x, b, false, true, false );
    // Libera memoria
    if( pos != NULL  )delete [] pos;
}
//............................................................................
void solver::reacomodaX( const int n, double *x, int *pos  ){
	// Reacomoda los elementos del vector solución para que estén de acuerdo
	// al orden que tenían originalmente las variables
	double *xAux;
	xAux = new double [n];
	// Hace una copia del vector solución
	for( int i = 0; i < n; i++ ) xAux[i] = x[i];
	for( int i = 0; i < n; i++ ) x[pos[i]] = xAux[i];
	delete xAux;
}
//............................................................................
int solver::mueveElemento( const int ri, const int ci, const int rf, const int cf,
                           const int m, const int n, double **A, double *b,
                           int *pos ){
    // Mueve un elemento de la matriz A de una posición inicial a una final
    // haciendo intercambios de renglones y de columnas. Regresa el número de
    // de cambios que se hicieron.
    // Variables de entrada
        // ri, rf    Renglón inicial y final, donde 0 <= ri,rf <= m-1
        // ci, cf    Columna inicial y final, donde 0 <= ri,rf <= n-1
        // A         Matriz, tam: m x n
        // pos       Guarda los cambios de renglones
    // Variables locales
        // conta     Cuenta el cambio de renglones y columnas realizados
    int conta = 0;
    double aux;
    // Hace el intercambio de renglones
    if( ri != rf ){
        for( int i = 0; i < m; i++ ){
            aux = A[rf][i];
            A[rf][i] = A[ri][i];
            A[ri][i] = aux;
        }
		aux = b[rf];
		b[rf] = b[ri];
		b[ri] = aux;
        conta++;
    }
    // Hace el intercambio de columnas
     if( ci != cf ){
        for( int i = 0; i < n; i++ ){
            aux = A[i][cf];
            A[i][cf] = A[i][ci];
            A[i][ci] = aux;
        }
        conta++;
        // Guarda el cambio de posición
        aux = pos[cf];
        pos[cf] = pos[ci];
        pos[ci] = aux;
    }
    return conta;
}
//............................................................................
void solver::valMax( const int ri, const int rf, const int ci, const int cf,
                     double **A, int &p, int &q ){
    // Obtiene la posición del elemento de la matriz A con el mayor valor
    // absoluto entre los renglones y columnas especificados
    // Variables de entrada
        // ri, rf    Renglón inicial y final, donde 0 <= ri,rf <= m-1
        // ci, cf    Columna inicial y final, donde 0 <= ri,rf <= n-1
        // A         Matriz, tam: m x n
    double max = 0.0;
    for( int i = ri; i <= rf; i++ ){
        for( int j = ci; j <= cf; j++ ){
            if( abs( A[i][j] ) > max ){
                max = abs( A[i][j] );
                p = i; // Renglón
                q = j; // Columna
            }
        }
    }
}
//............................................................................
//                             Factorización LU
//............................................................................
void solver::LU( int n, double **A, double *x, double *b, int solver,
                 string salida ){
	// Resuelve un sistema de ecuaciones lineales de la forma Ax = b, fac-
	// torizando la matriz A = LU, ya sea por el método de Crout o de Doolittle
	// Factoriza A de acuerdo al método especificado
	if( solver == 1 ) LU_Crout_mejorado( n, A ); // Este método debe actua-
	                                                   // lizar a A
	else ;
	// Obtiene y de Ly = b
	miVector <double> y( n );
	resuelveSistTriangINF( n, A, y.val, b );
	// Obtiene x de Ux = y
	for( int i = 0; i < n; i++ ) A[i][i] = 1.0;
	resuelveSistTriangSUP( n, A, x, y.val );
	// Imprime resultados
	//imprime( n, A, x, b, false, true, false );
	escribe_LU( salida, n, A );
}
//............................................................................
void solver::factoriza_LU( int n, double **A ){
    // Factoriza la matriz A y el resultado lo guarda en la misma
    LU_Crout_mejorado( n, A );
}
void solver::resuelve_LU( int n, double **A, double *x, double *b ){
	// Resuelve un sistema de ecuaciones lineales de la forma Ax = b,
	// donde A debe contener su facorización LU
	// Obtiene y de Ly = b
	miVector <double> y( n );
	miVector <double> aux( n );
	resuelveSistTriangINF( n, A, y.val, b );
	// Obtiene x de Ux = y
	for( int i = 0; i < n; i++ ){
        aux.val[i] = A[i][i];
        A[i][i] = 1.0;
    }
	resuelveSistTriangSUP( n, A, x, y.val );
	for( int i = 0; i < n; i++ ) A[i][i] = aux.val[i];
}
//............................................................................
void solver::LU_Crout( int n, double **A, double *x, double *b ){
	// Factoriza la matriz A = LU, donde L es una matriz triangular inferior
	// y U una triangular superior.
	double sum, sumL, sumU;
	matriz <double> L( n, n );
	matriz <double> U( n, n );
	L.val[0][0] = A[0][0];
	U.val[0][0] = 1.0;
	for( int i = 1; i < n; i++ ){
		for( int j = 0; j < i; j++){
			sumU = 0.0;
			sumL = 0.0;
			for( int k = 0; k < j; k++ ){
				sumU += L.val[j][k] * U.val[k][i];
				sumL += L.val[i][k] * U.val[k][j];
			}
			U.val[j][i] = ( A[j][i] - sumU ) / L.val[j][j];
			L.val[i][j] =   A[i][j] - sumL;
		}
		sum = 0.0;
		for( int k = 0; k < i; k++ ) sum += L.val[i][k] * U.val[k][i];
		L.val[i][i] = A[i][i] - sum;
		U.val[i][i] = 1.0;
	}
	//imprime( n, L.val, x, b, true, false, false );
	//imprime( n, U.val, x, b, true, false, false );
}
void solver::LU_Crout_mejorado( int n, double **A ){
	// Factoriza la matriz A = LU, donde L es una matriz triangular inferior
	// y U una triangular superior.
	// NOTA: la mejora está en que todas la factorización la realiza sobre
	// la matriz A
	double sum, sumL, sumU;
	for( int i = 1; i < n; i++ ){
		for( int j = 0; j < i; j++){
			sumU = 0.0;
			sumL = 0.0;
			for( int k = 0; k < j; k++ ){
				sumU += A[j][k] * A[k][i];
				sumL += A[i][k] * A[k][j];
			}
			A[j][i] = ( A[j][i] - sumU ) / A[j][j];
			A[i][j] =   A[i][j] - sumL;
		}
		sum = 0.0;
		for( int k = 0; k < i; k++ ) sum += A[i][k] * A[k][i];
		A[i][i] = A[i][i] - sum;
	}
	//imprime( n, A, x, b, true, false, false );
}
void solver::escribe_LU( string salida, int n, double **LU ){
	// Escribe las matrices L y U al archivo especificado
	string matL( salida );
	matL.insert( matL.length( ), "L.mat" );
	string matU( salida );
	matU.insert( matU.length( ), "U.mat" );
  	ofstream escribeL( matL.c_str( ), std::ofstream::out );
  	ofstream escribeU( matU.c_str( ), std::ofstream::out );
	if( escribeL.is_open( ) && escribeU.is_open( ) ){
		// Dimensiones de las matrices
		escribeL << n << " " << n << endl;
		escribeU << n << " " << n << endl;
		// Escribe las matrices
		for( int i = 0; i < n; i++ ){
			for( int j = 0; j < n; j++ ){
				if( j <= i ) escribeL << setw(10) << setprecision(5) << LU[i][j] << " ";
				else         escribeL << setw(10) << setprecision(5) << 0.0 << " ";
				if( j == i )     escribeU << setw(10) << setprecision(5) << 1.0 << " ";
				else if( j > i ) escribeU << setw(10) << setprecision(5) << LU[i][j] << " ";
				else             escribeU  << setw(10) << setprecision(5) << 0.0 << " ";
			}
			escribeL << endl;
			escribeU << endl;
		}
		escribeL.close( );
		escribeU.close( );
	}
}
//............................................................................
//                           Factorización LDU
//............................................................................
void solver::LDU( int n, double **A, double*x, double *b, string salida ){
	// Resuelve un sistema de ecuaciones lineales de la forma Ax=b.
	// Para ello factoriza la matriz A = LDU.
	// Factoriza a A y el resultado lo guarda en el mismo espacio de A. En A
	// estaría guardado L y D
	cholesky_mejorado( n, A );
	// Crea una matriz auxiliar para guardar el producto LD
	matriz <double> B( n, n );
	// Hace el producto LD
	mult_LD( n, A, B.val );
	// Obtiene y de LDy = b
	miVector <double> y( n );
	resuelveSistTriangINF( n, B.val, y.val, b );
	// Usa la matriz auxiliar para colocar L'
	L_sup( n, A, B.val );
	// Obtiene x de L'x = y
	resuelveSistTriangSUP( n, B.val, x, y.val );
	// Imprime resultados
	//imprime( n, A, x, b, false, true, false );
	// Escribe los resultados
	escribe_LD( salida, n, A );
}
//............................................................................
void solver::L_sup( int n, double **A, double **L_Sup ){
	// Coloca la matriz L' en L_sup. L está almacenada en A, pero sin los 1's
	// de la diagonal
	for( int i = 0; i < n; i++ ){
		for( int j = 0; j <= i; j++ ){
			if( i == j ) L_Sup[i][i] = 1.0;
			else         L_Sup[j][i] = A[i][j];
		}
	}
}
//............................................................................
void solver::mult_LD( int n, double **A, double **LD ){
	// Hace la multiplicación de LxD, considerando que L y D están almacenadas
	// en A. De este producto resulta una matriz triangular inferior.
	for( int i = 0; i < n; i++ ){
		for( int j = 0; j <= i; j++ ){
			if( j == i ) LD[i][i] = A[i][i];
			else         LD[i][j] = A[i][j] * A[j][j];
		}
	}
}
//............................................................................
void solver::cholesky( int n, double **A  ){
	// Factoriza la matriz A = LDL', donde L es una matriz triangular inferior.
	// L tienen 1's en su diagonal. D es una matriz diagonal.
	// NOTA: este método sólo funciona para cuando A es simétrica.
	double sum;
	// Esta matriz contendrá a L y a D
	matriz <double> L( n, n );
	matriz <double> D( n, n );
	for( int i = 0; i < n; i++ ){
		for( int j = 0; j <= i; j++ ){
			// Elementos de la diagonal
			if( i == j ){
				sum = 0.0;
				for( int k = 0; k < i; k++ ) sum += L.val[i][k] * L.val[i][k] * D.val[k][k];
				D.val[i][i] = A[i][i] - sum;
			}
			// Elementos fuera de la diagonal
			else{
				sum = 0.0;
				for( int k = 0; k < j; k++ ) sum += L.val[i][k] * L.val[j][k] * D.val[k][k];
				L.val[i][j] = ( A[i][j] - sum ) / D.val[j][j];
			}
		}
	}
	for( int i = 0; i < n; i++ ) L.val[i][i] = 1.0;
	// Copia los valores a la submatriz triangular superior
	for( int i = 0; i < n; i++ )
	 	for( int j = n - 1; j > i; j-- )
	 		L.val[i][j] = L.val[j][i];
}
void solver::cholesky_mejorado( int n, double **A ){
	// Factoriza la matriz A = LDL', donde L es una matriz triangular inferior.
	// L tienen 1's en su diagonal. D es una matriz diagonal.
	// NOTA: este método sólo funciona para cuando A es simétrica.
	// NOTA: guarda la matriz L y D sobre A
	double sum;
	for( int i = 0; i < n; i++ ){
		for( int j = 0; j <= i; j++ ){
			// Elementos de la diagonal
			if( i == j ){
				sum = 0.0;
				for( int k = 0; k < i; k++ ) sum += A[i][k] * A[i][k] * A[k][k];
				A[i][i] = A[i][i] - sum;
			}
			// Elementos fuera de la diagonal
			else{
				sum = 0.0;
				for( int k = 0; k < j; k++ ) sum += A[i][k] * A[j][k] * A[k][k];
				A[i][j] = ( A[i][j] - sum ) / A[j][j];
			}
		}
	}
	// Imprime resultados
	//imprime( n, A, x, b, false, true, false );
}
void solver::escribe_LD( string salida, int n, double **LD ){
	// Escribe las matrices L y D al archivo especificado
	string matL( salida );
	matL.insert( matL.length( ), "L.mat" );
	string matD( salida );
	matD.insert( matD.length( ), "D.mat" );
  	ofstream escribeL( matL.c_str( ), std::ofstream::out );
  	ofstream escribeD( matD.c_str( ), std::ofstream::out );
	if( escribeD.is_open( ) && escribeL.is_open( ) ){
		// Dimensiones de las matrices
		escribeD << n << endl;
		escribeL << n << " " << n << endl;
		// Escribe las matrices
		for( int i = 0; i < n; i++ ){
			escribeD << setw(10) << setprecision(5) << LD[i][i] << endl;
			for( int j = 0; j < n; j++ ){
				if( j == i )     escribeL << setw(10) << setprecision(5) << 1.0 << " ";
				else if( j < i ) escribeL << setw(10) << setprecision(5) << LD[i][j] << " ";
				else         escribeL << setw(10) << setprecision(5) << 0.0 << " ";
			}
			escribeL << endl;
		}
		escribeL.close( );
		escribeD.close( );
	}
}
//............................................................................
//                          Factorización QR
//............................................................................
void solver::QR( int n, double **A, double *x, double *b ){
	matriz <double> R( n, n );
	// Guarda a [Q] en [A]
	factoriza_QR( n, A, A, R.val );
	resuelve_QR ( n, A, R.val, x, b );
}
//............................................................................
void solver::factoriza_QR( int n, double **A, double **Q, double **R ){
	// Factoriza a [A] = [Q][R] donde [Q] es una matriz con vectores-columna
	// ortonormales y [R] es una matriz triangular inferior.
	double sum;
	miVector <double> v( n );
	for( int j = 0; j < n; j++ ){
		// Inicializa {v} = {a_j} - sum( R_ij{q_i} )
		for( int k = 0; k < n; k++ ) v.val[k] = A[k][j];
		// Calcula los R_ij y los R_ii
		for( int i = 0; i <= j - 1; i++ ){
			sum = 0.0;
			for( int k = 0; k < n; k++ ) sum += Q[k][i] * A[k][j];
			R[i][j] = sum;
			for( int k = 0; k < n; k++ ) v.val[k] -= R[i][j] * Q[k][i];
		}
		// Calcula los R_ii
		sum = 0.0;
		for( int k = 0; k < n; k++ ) sum += v.val[k] * v.val[k];
		R[j][j] = sqrt( sum );
		// Calcula los q_i
		for( int k = 0; k < n; k++ ){
			Q[k][j] = v.val[k] / R[j][j];
		}
	}
	imprime_QR( n, Q, R );
}
//............................................................................
void solver::imprime_QR( int n, double **Q, double **R ){
	// Imprime la factorización QR
    cout << "Matriz [Q]" << endl;
    for( int i = 0; i < n; i++ ){
        for( int j = 0; j < n; j++ )
            cout << setw( LONG_NUM ) << setprecision( 5 ) << Q[i][j] << "  ";
            cout << endl;
	}
	cout << "Matriz [R]" << endl;
    for( int i = 0; i < n; i++ ){
        for( int j = 0; j < n; j++ )
            cout << setw( LONG_NUM ) << setprecision( 5 ) << R[i][j] << "  ";
            cout << endl;
	}
}
//............................................................................
void solver::resuelve_QR( int n, double **Q, double **R, double *x, double *b ){
	// Resuelve sistemas de ecs. linealdes de la forma [A]{x}={b} factorizando
	// la matriz [A] = [Q][R] donde [Q] es una matriz con vectores-columna
	// ortonormales y [R] es una matriz triangular inferior.
	// b' =	Q^(-1) * b = Q^T * b
	matriz <double> operMat;
	miVector <double> bb( n );
	operMat.transpuesta( n, Q );
	operMat.mult_mat_vec( n, n, Q, b, bb.val );
	resuelveSistTriangSUP( n, R, x, bb.val );
	// Regresa [Q] a su origen por si se tiene que resolver para otro {b}
	operMat.transpuesta( n, Q );
}
//............................................................................
//                          Métodos iterativos
//............................................................................
void solver::jacobi( int n, double **A, double *x, double *b ){
	// Calcula la solución de un sistema de ecs. lineales de la forma Ax=b
	// mediante el método iterativo Jacobi.
	// NOTA: Para asegurar un buen funcionamiento A tiene que ser diagonalmente
	// dominante y si no x0 debe de estar cerca de la solución.
	double iniVal, tol, sum;
	int conta;
	miVector <double> x0( n );
	// Inicia variables
	tol    = 0.00000000001;
	iniVal = 0.0;
	conta  = -1;
	// Inicializa vectores
	iniciaVec_jacobi( n, x0.val, x, iniVal );
	cout << "Iteracion " << " Vector x' " << endl;
	while( !jacobi_convergencia( n, x0.val, x, tol ) ){
		// Actualiza el vector, excepto para la 1er iteracion
		if( conta != 0 )
			for( int i = 0; i < n; i++ ) x0.val[i] = x[i];
		for( int i = 0; i < n; i++ ){
			sum = 0.0;
			for( int j = 0; j < i; j++ )     sum += A[i][j] * x0.val[j];
			for( int j = i + 1; j < n; j++ ) sum += A[i][j] * x0.val[j];
			// Nuevo vector
			x[i] = ( b[i] - sum ) / A[i][i];
			//imprime( n, A, x, b, false, true, false );
		}
		// Imprime el vector actual en pantalla
		if( conta >= 0 ){
			for( int i = 0; i < n; i++ ) cout << conta << " " <<  setw( 10 ) <<
			 				  	         setprecision( 3 ) <<  x[i];
			cout << endl;
		}
		conta++;
		//if( conta > 1000000 ) break;
	}
	//imprime( n, A, x, b, false, true, false );
}
void solver::iniciaVec_jacobi( int n, double *x0, double *x, const double val ){
	// Inicia los vectores para jacobi
	for( int i = 0; i < n; i++ ){
		x0[i] = val;
		x[i]  = 2.0 * val + 100.0;
	}
}
bool solver::jacobi_convergencia( int n, double *x0, double *x, double tol ){
	// Checa si el método de jacobi ya convergió
	double sum1, sum2;
	bool converge;
	sum1 = 0.0;
	sum2 = 0.0;
	for( int i = 0; i < n; i++ ){
		sum1 += ( x[i] - x0[i] ) * ( x[i] - x0[i] );
		sum2 += x[i] * x[i];
	}
	if( sum1 / sum2 <= tol ) converge = true;
	else                     converge = false;
	return converge;
}
//............................................................................
void solver::gauss_seidel( int n, double **A, double *x, double *b ){
	// Calcula la solución de un sistema de ecs. lineales de la forma Ax=b
	// mediante el método iterativo Gauss-Seidel.
	// NOTA: este método requiere que la matriz sea diagolamente
	// dominante para funcionar correctamente.
	double iniVal, tol, n1, n2, sum1, sum2;
	int conta;
	// Inicia variables
	tol    = 0.00000000001;
	iniVal = 0.0;
	conta  = 0;
	// Inicializa vectores
	iniciaVec_gauss_seidel( n, x, iniVal );
	// Cantidades para la condición de para
	n2 = norma_gauss_seidel( n, x );
	n1 = n2 * 2.0 + 5.0;
	for( int i = 0; i < n * ANCHO_IMP; i++ ) cout << "-"; cout << endl;
	cout << "Iteracion " << " Vector x' " << endl;
	while( abs( n1 - n2 ) / n2 > tol ){
		for( int i = 0; i < n; i++ ){
			sum1 = 0.0;
			sum2 = 0.0;
			for( int j = 0; j < i; j++ )     sum1 += A[i][j] * x[j];
			for( int j = i + 1; j < n; j++ ) sum2 += A[i][j] * x[j];
			// Nuevo entrada del vector
			x[i] = ( b[i] - sum1 - sum2 ) / A[i][i];
		}
		// Actualiza las n's
		n1 = n2;
		n2 = norma_gauss_seidel( n, x );
		// Imprime el vector actual en pantalla
		imprime_iter( conta, n, x );
		conta++;
		//if( conta > 1000000 ) break;
	}
}
void solver::imprime_iter( int conta, int n, double *x ){
	// Imprime los valores del vector {x} de la iteración i
	cout << conta;
	for( int i = 0; i < n; i++ ) cout << setw( LONG_NUM  ) <<
					  	         setprecision( 5 ) <<  x[i];
	cout << endl;
}
void solver::iniciaVec_gauss_seidel( int n, double *x, const double val ){
	// Inicia los vectores para Gauss-Seidel
	for( int i = 0; i < n; i++ ) x[i] = val;
}
double solver::norma_gauss_seidel( int n, double *x ){
	// Calcula la el producto punto del vector consigo mismo
	double sum;
	sum = 0.0;
	for( int i = 0; i < n; i++ ) sum += x[i] * x[i];
	return sum;
}
//............................................................................
void solver::gradienteConjugadoMejorado( int n, double **A, double *x, double *b ){
	// Resuelve sistemas de la fora [A]{x} = {b}. [A] tiene que ser simétrica
	// y definida positiva. Esta es la versión buena.
	miVector <double> r( n ); // Residuo
	miVector <double> p( n ); // Vector ortogonal a {x_i}
	miVector <double> w( n );  // {w} = [A]{p0}
	int k;          // Contador, en teoría debería k <= n
	double a;       // Factor para generar nueva sol: {x1} = a{p0}
	double B;       // Factor para generar nueva dir: {p1} = -r1 + B{p0}
	double sum, numeA, denoA, numeB;
	double tol = 1e-15, fMax = 1.5;
	int maxIter = static_cast < int > ( fMax * n );
	// NOTA: no se usaron los métodos de la clase matriz y vector con la posi-
	// bilidad de paralelizar el algoritmo aquí mismo.
	// INICIALIZA
	// {x0} = {b}
	for( int i = 0; i < n; i++ ) x[i] = b[i];
	// {r0} = [A]{x0} - {b}
	// {p0} = -{r0}
	numeA = 0.0;
	for( int i = 0; i < n; i++ ){
		sum = 0.0;
		for( int j = 0; j < n; j++ ) sum += A[i][j] * x[j];
		r.val[i] = sum - b[i];
		p.val[i] = -r.val[i];
		numeA += r.val[i] * r.val[i];
	}
	// ITERA
	k = 0;
	/*for( int i = 0; i < n * ANCHO_IMP; i++ ) cout << "-"; cout << endl;
	cout << "Iteracion " << " Vector x' " << endl;
	imprime_iter( k, n, x );*/
	while( k < maxIter ){
		if( abs( numeA ) < tol ) break;
		// {w} = [A]{p0}
		// a = ( {r0}^T {r0} ) / ( {p0}^T {w} )
		denoA = 0.0;
	    for( int i = 0; i < n; i++ ){
			sum = 0.0;
			for( int j = 0; j < n; j++ ) sum += A[i][j] * p.val[j];
			w.val[i] = sum;
			denoA += p.val[i] * w.val[i];
	    }
		a = numeA / denoA;
		// {x1} = {x0} + a {p0}
		// {r1} = {r0} + a {w}
		// B = ( {r1}^T {r1} ) / ( {r0}^T {r0} )
		numeB = 0.0;
		for( int i = 0; i < n; i++ ){
			x[i]      += a * p.val[i];
			r.val[i]  += a * w.val[i];
			numeB     += r.val[i] * r.val[i];
		}
		B = numeB / numeA;
		// {p1} = -{r1} + B{p0}
		for( int i = 0; i < n; i++ ) p.val[i] = -r.val[i] + B * p.val[i];
		// Actualiza variables
		numeA = numeB;
		k++;
		// Imprime iteración
		//imprime_iter( k, n, x );
	}
}
void solver::gradienteConjugadoParalel( int n, double **A, double *x, double *b,
                                        int hilos ){
	// Resuelve sistemas de la fora [A]{x} = {b}. [A] tiene que ser simétrica
	// y definida positiva. Esta es la versión buena.
	miVector <double> r( n ); // Residuo
	miVector <double> p( n ); // Vector ortogonal a {x_i}
	miVector <double> w( n );  // {w} = [A]{p0}
	int k;          // Contador, en teoría debería k <= n
	double a;       // Factor para generar nueva sol: {x1} = a{p0}
	double B;       // Factor para generar nueva dir: {p1} = -r1 + B{p0}
	double sum, numeA, denoA, numeB;
	double tol = 1e-15, fMax = 1.5;
	int maxIter = static_cast < int > ( fMax * n );
	// NOTA: no se usaron los métodos de la clase matriz y vector con la posi-
	// bilidad de paralelizar el algoritmo aquí mismo.
	// INICIALIZA
	// {x0} = {b}
	//unsigned int threads;

	//threads = static_cast < int > ( hilos );
	//omp_set_num_threads( threads );

	for( int i = 0; i < n; i++ ) x[i] = b[i];
	// {r0} = [A]{x0} - {b}
	// {p0} = -{r0}

	numeA = 0.0;

	for( int i = 0; i < n; i++ ){
		sum = 0.0;
		for( int j = 0; j < n; j++ ) sum += A[i][j] * x[j];
		r.val[i] = sum - b[i];
		p.val[i] = -r.val[i];
		numeA += r.val[i] * r.val[i];
	}

	// ITERA
	k = 0;
	/*for( int i = 0; i < n * ANCHO_IMP; i++ ) cout << "-"; cout << endl;
	cout << "Iteracion " << " Vector x' " << endl;
	imprime_iter( k, n, x );*/
	while( k < maxIter ){
		if( abs( numeA ) < tol ) break;
		// {w} = [A]{p0}
		// a = ( {r0}^T {r0} ) / ( {p0}^T {w} )
		denoA = 0.0;
		//#pragma omp parallel for default(shared) private(sum) reduction(+:denoA)
	    for( int i = 0; i < n; ++i ){
			sum = 0.0;
			for( int j = 0; j < n; j++ ) sum += A[i][j] * p.val[j];
			w.val[i] = sum;
			denoA += p.val[i] * w.val[i];
	    }
		a = numeA / denoA;

		// {x1} = {x0} + a {p0}
		// {r1} = {r0} + a {w}
		// B = ( {r1}^T {r1} ) / ( {r0}^T {r0} )
		numeB = 0.0;
		//#pragma omp parallel for default(shared) private(sum) reduction(+:numeB)
		for( int i = 0; i < n; i++ ){
			x[i]      += a * p.val[i];
			r.val[i]  += a * w.val[i];
			numeB     += r.val[i] * r.val[i];
		}
		B = numeB / numeA;

		// {p1} = -{r1} + B{p0}
		//#pragma omp parallel for default(shared)
		for( int i = 0; i < n; i++ ) p.val[i] = -r.val[i] + B * p.val[i];
		// Actualiza variables
		numeA = numeB;
		k++;
		// Imprime iteración
		//imprime_iter( k, n, x );
	}
}
//............................................................................
//                    Sistemas diagonales
//............................................................................
void solver::sistema3diagonalSim( int n, double *V, double *B, double *x,
                                  double *b ){
	// Resuelve un sistema tridiagonal simétrico.
	// En este caso la matriz se almacena en vectores, pero la matriz de coe-
	// ficientes se factoriza en la forma LDL'. D es na matriz diagonal y L
	// es una matriz triangular inferior con 1's en su diagonal.
	// {V} contiene los datos de la diagonal
	miVector <double> D( n );
	miVector <double> L( n );
	// Hace la factorización  LDL'
	D.val[0] = V[0];
	L.val[0] = B[0] / D.val[0];
	for( int i = 1; i < n; i++ ){
		D.val[i]   = V[i] - D.val[i-1] * L.val[i-1] * L.val[i-1];
		L.val[i] = B[i] / D.val[i];
	}
	// Resuelve para Ly = b
	// En este caso y_i = ( b_i - L_(i-1) * D_(i-1) ) / D_(i)
	miVector <double> y( n );
	y.val[0] = b[0] / D.val[0];
	for( int i = 1; i < n; i++ )
		y.val[i] = ( b[i] - L.val[i-1] * D.val[i-1] * y.val[i-1] ) / D.val[i];
	// Resuelve para L'x = y
	// En este caso x_i = y_i - L_i * x_(i+1)
	x[n-1] = y.val[n-1];
	for( int i = n - 2; i >= 0; i-- )
		x[i] = y.val[i] - L.val[i] * x[i+1];
	//cout << "Vector x" << endl;
	//for( int i = 0; i < n; i++ ) cout << x[i] << endl;
}
//............................................................................
//                    Solución de sistemas triangulares
//............................................................................
void solver::resuelveSistTriangSUP( int n, double **A, double *x,
				                    double *b ){
    // Realiza sustitución hacia atrás para resolver un sistema
    // triangular superior
    double sum;
    for( int i = n - 1; i >= 0; i-- ){
        sum = 0.0;
        for( int j = i + 1; j < n; j++ ) sum += A[i][j] * x[j];
        x[i] = ( b[i] - sum ) / A[i][i];
    }
}
//............................................................................
void solver::resuelveSistTriangINF( int n, double **A, double *x,
				                    double *b ){
    // Realiza sustitución hacia atrás para resolver un sistema
    // triangular inferior
    double sum;
    for( int i = 0; i < n; i++ ){
        sum = 0.0;
        for( int j = 0; j < i; j++ ) sum += A[i][j] * x[j];
        x[i] = ( b[i] - sum ) / A[i][i];
    }
}
//............................................................................
// Métodos para impresión
//............................................................................
void solver::imprime( int n, double **A, double *x, double *b,
                      const bool impA, const bool impX, const  bool impB  ){
    // Imprime los resultados de la matriz A y el vector b
	if( impA ){
	    cout << "Matriz [A]" << endl;
	    for( int i = 0; i < n; i++ ){
	        for( int j = 0; j < n; j++ )
	            cout << setw( 13 ) << setprecision( 3 ) << A[i][j] << "  ";
	            cout << endl;
	    }
	}
	if( impX ){
	    cout << endl << "Vector solucion {x}" << endl;
	    for( int i = 0; i < n; i++ )
	    	cout << setw( 15 ) << setprecision( 5 ) << x[i] << endl;
	}
	if( impB ){
	    cout << endl << "Vector de terminos independiente {b}" << endl;
	    for( int i = 0; i < n; i++ )
			cout << setw( 15 ) << setprecision( 3 ) << b[i] << endl;
	}
}
//............................................................................
//============================================================================
