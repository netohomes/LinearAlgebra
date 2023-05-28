#include "eigen_val_vec.h"
#include <cmath>
#include "matrix.h"
#include "linear_solvers.h" // Contiene a iostream
#include <cstdlib>
#include <iomanip>
static int ANCHO_IMP = 22;
static int LONG_NUM  = 15;
using namespace std;
//============================================================================
// MÉTODOS PARA OBTENER LOS EIGENVALORES Y EIGENVECTORES
//============================================================================
eigen_vec_val:: eigen_vec_val( ){}
eigen_vec_val::~eigen_vec_val( ){}
//...........................................................................
//                          Método de la POTENCIA
//...........................................................................
// Con este método se pueden obtener los "p" eigenvalores y eigenvectores
// más grandes (en valor absoluto) de una matriz<double>
//...........................................................................
void eigen_vec_val::potencia( int n, double **A, int p, int tipo, bool imprime ){
	miVector<double> v0( n ), v1( n );
	miVector<double> eigenVal( n );        // Guarda los eigenvalores
	matriz<double> eigenVec( p, n + 1 ); // Guarda los eigenvectores en renglones
	if( tipo == 1 ){
		// Se resta la contribución de los otros eigenvectores restando
		// en cada iteración al miVector<double> resultante
		for( int i = 0; i < p; i++ ){
			for( int j = 0; j < n; j++ ) v0.val[j] = 0.0;
			v0.val[i] = 1.0;
			potencia_i( n, A, v0.val, v1.val, eigenVal.val, eigenVec.val, i,
			            imprime );
		}
	}
	else{
		// Se resta la contribución de los otros eigenvectores haciendo
		// una deflación de la matriz<double> [A] (haciendo cero el eigenvalor
		// que ya ha sido encontrado)
		for( int i = 0; i < p; i++ ){
			for( int j = 0; j < n; j++ ) v0.val[j] = 0.0;
			v0.val[i] = 1.0;
			potencia_i2( n, A, v0.val, v1.val, eigenVal.val, eigenVec.val, i,
			             imprime );
		}
	}
}
//...........................................................................
void eigen_vec_val::restaEigen( int n, int p, double **eigenVec, double *v0 ){
	// Resta la constribución de los "p-1" eigenvectores y eigenvalores
	// anteriores a v0. Los eigenvectores están guardados en renglones.
	double a = 0.0;
	for( int i = 0; i < p; i++ ){
		// Calcula el a_i
		a = 0.0;
        for( int j = 0; j < n; j++ ) a += eigenVec[i][j] * v0[j];
		// Le resta la contribución a v0
		for( int j = 0; j < n; j++ ) v0[j] -= a * eigenVec[i][j];
    }
}
//...........................................................................
void eigen_vec_val::restaEigen2( int n, double **A, double *eigenVal, 
                                 double **eigenVec, int t ){
	// Modifica la matriz<double> A para hacer 0 el eigenvalor ya encontrado
	// Los eigenvectores deben estar normalizados
	for( int i = 0; i < n; i++ ){
		for( int j = 0; j < n; j++ ){
		  A[i][j] = A[i][j] - eigenVal[t] * eigenVec[t][i] * eigenVec[t][j];
		}
	}
}
//...........................................................................
bool eigen_vec_val::potencia_i( int n, double **A, double *V0, double *V1,
                                double *eigenVal, double **eigenVec,  int p,
                                bool imprime ){
	double tol, nume, deno, lambda, lambda_old;
	int conta;
	tol = 0.0000000001;
	matriz<double> operMat;
	miVector<double> operVec;
	// Primer lambda_old y  V0
	restaEigen( n, p, eigenVec, V0 );
	normaliza( n, V0 );
	operMat.mult_mat_vec( n, n, A, V0, V1 );
 	nume = operVec.produc_punto( n, V1, V1 );
	deno = operVec.produc_punto( n, V1, V0 );
	operVec.copia_val( n, V0, V1 );
	lambda_old = nume / deno;
	// Primer lambda y V1
	restaEigen( n, p, eigenVec, V0 );
	normaliza( n, V0 );
	operMat.mult_mat_vec( n, n, A, V0, V1 );
 	nume = operVec.produc_punto( n, V1, V1 );
	deno = operVec.produc_punto( n, V1, V0 );
	lambda = nume / deno;
	conta = 1;
	// Impresión de resultados
	if( imprime ){
		for( int i = 0; i < n * ANCHO_IMP; i++ ) cout << "="; cout << endl;
		cout << "Iter  Eigenvalor  Eigenvector (transpuesto)" << endl;
		imprimeIter( lambda_old, n, V0, 0 );
		imprimeIter( lambda, n, V1, conta );
	}
	// Más aproximaciones
	while( abs( lambda - lambda_old ) > tol ){
		// Actualiza variables
		lambda_old = lambda;
		operVec.copia_val( n, V0, V1 );
		// Quita contribución de eigenvectores anteriores
		restaEigen( n, p, eigenVec, V0 );
		normaliza( n, V0 );
		// Se calcula una nueva aproximación
		operMat.mult_mat_vec( n, n, A, V0, V1 );
	 	nume = operVec.produc_punto( n, V1, V1 );
		deno = operVec.produc_punto( n, V1, V0 );
		lambda = nume / deno;
		conta++;
		normaliza( n, V1 );
		if( imprime ) imprimeIter( lambda, n, V1, conta );
	}
	// Guarda el eigenvector y el eigenvalor
	eigenVal[p] = lambda_old;
	for( int i = 0; i < n; i++ ) eigenVec[p][i] = V0[i];
	imprimeRes( lambda_old, n, V0, conta );
	return true;
}
//...........................................................................
bool eigen_vec_val::potencia_i2( int n, double **A, double *V0, double *V1,
                                 double *eigenVal, double **eigenVec, int p,
                                 bool imprime ){
	// Este método en lugar de restar a V0 la contribución de los eigenvectores
	// ya calculados, modifica la matriz<double> [A] de forma que hace 0 los eigenvalores
	// que ya han sido calculados
	double tol, nume, deno, lambda, lambda_old;
	int conta;
	tol = 0.0000000001;
	matriz<double> operMat;
	miVector<double> operVec;
	// Primer lambda_old y  V0
	if( p > 0 ) restaEigen2( n, A, eigenVal, eigenVec, p - 1 );
	normaliza( n, V0 );
	operMat.mult_mat_vec( n, n, A, V0, V1 );
 	nume = operVec.produc_punto( n, V1, V1 );
	deno = operVec.produc_punto( n, V1, V0 );
	lambda_old = nume / deno;
	normaliza( n, V1 );
	operVec.copia_val( n, V0, V1 );
	// Primer lambda y V1
	operMat.mult_mat_vec( n, n, A, V0, V1 );
 	nume = operVec.produc_punto( n, V1, V1 );
	deno = operVec.produc_punto( n, V1, V0 );
	lambda = nume / deno;
	normaliza( n, V1 );
	conta = 1;
	// Impresión de resultados
	for( int i = 0; i < n * ANCHO_IMP; i++ ) cout << "="; cout << endl;
	if( imprime ){
		cout << "Iter  Eigenvalor  Eigenvector (transpuesto)" << endl;
		imprimeIter( lambda_old, n, V0, 0 );
		imprimeIter( lambda, n, V1, conta );
	}
	// Más aproximaciones
	while( abs( lambda - lambda_old ) > tol ){
		// Actualiza variables
		lambda_old = lambda;
		operVec.copia_val( n, V0, V1 );
		// Se calcula una nueva aproximación
		operMat.mult_mat_vec( n, n, A, V0, V1 );
	 	nume = operVec.produc_punto( n, V1, V1 );
		deno = operVec.produc_punto( n, V1, V0 );
		lambda = nume / deno;
		conta++;
		normaliza( n, V1 );
		if( imprime ) imprimeIter( lambda, n, V1, conta );
	}
	// Guarda el eigenvector y el eigenvalor
	eigenVal[p] = lambda;
	for( int i = 0; i < n; i++ ) eigenVec[p][i] = V1[i];
	imprimeRes( lambda, n, V1, conta );
	return true;
}
//...........................................................................
bool eigen_vec_val::normaliza( int n, double *V ){
	// Normaliza un miVector<double> y regresa si se pudo realizar la normalización
	bool normalizado;
	double sum = 0.0;
	for( int i = 0; i < n; i++ ) sum += V[i] * V[i];
	sum = sqrt( sum );
 	if( sum > 0.000000000000000000001 ){
		normalizado = true;
		for( int i = 0; i < n; i++ ) V[i] /= sum;
	}
	else{
		//for( int i = 0; i < n * ANCHO_IMP; i++ ) cout << "-"; cout << endl;
		cout << endl << "ERROR: el miVector<double> no se pudo normalizar. ";
		cout << "TERMINANDO EL PROGRAMA..." << endl;
		//for( int i = 0; i < n * ANCHO_IMP; i++ ) cout << "-";
		cout << endl;
		normalizado = false;
		exit( EXIT_FAILURE );
	}
	return normalizado;
}
//...........................................................................
void eigen_vec_val::imprimeRes( double lambda, int n, double *v, int iter ){
	// Imprime el eigenvalor y su correspondiente eigenvector
	for( int i = 0; i < n * ANCHO_IMP; i++ ) cout << "-"; cout << endl;
	cout << "Eigenvalor    = " << lambda << endl;
	cout << "Num. de iter. = " << iter << endl;
	//for( int i = 0; i < n * ANCHO_IMP; i++ ) cout << "-";
	//cout << endl;
	cout << "Eigenvector:" << endl;
	for( int i = 0; i < n; i++ ) cout << setprecision(3) << v[i] << "    ";
	cout << endl;
}
//...........................................................................
void eigen_vec_val::imprimeIter( double lambda, int n, double *v, int iter ){
	// Imprime los valores de la iteración
	//for( int i = 0; i < n * ANCHO_IMP; i++ ) cout << "-"; cout << endl;
	cout << iter << setw( 12 ) << lambda;
	for( int i = 0; i < n; i++ )
		cout << "  " << setw( LONG_NUM ) << setprecision( 5 ) << v[i];
	cout << endl;
	//for( int i = 0; i < n * ANCHO_IMP; i++ ) cout << "-"; cout << endl;
}
//...........................................................................
//                       Método de la POTENCIA INVERSA
//...........................................................................
// Con este método se pueden obtener los "p" eigenvalores y eigenvectores
// más pequeños (en valor absoluto) de una matriz<double>
//...........................................................................
void eigen_vec_val::potencia_inversa( int n, double **A, int p, int tipo,
                                      bool imprime ){
	miVector<double> v0( n ), v1( n );
	solver resuelve;
	miVector<double> eigenVal( n );        // Guarda los eigenvalores
	matriz<double> eigenVec( p, n + 1 ); // Guarda los eigenvectores en renglones
    // Factoriza la matriz<double>
    resuelve.factoriza_LU( n, A );
	if( tipo == 1 ){
		for( int i = 0; i < p; i++ ){
			for( int j = 0; j < n; j++ ) v0.val[j] = 0.0;
			v0.val[i] = 1.0;
			potencia_inv_i( n, A, v0.val, v1.val, eigenVal.val, eigenVec.val, i,
			                imprime );
		}
	}
	else{
		for( int i = 0; i < p; i++ ){
			for( int j = 0; j < n; j++ ) v0.val[j] = 0.0;
			v0.val[i] = 1.0;
			potencia_inv_i2( n, A, v0.val, v1.val, eigenVal.val, eigenVec.val, i,
			                imprime );
		}

	}
}
//...........................................................................
bool eigen_vec_val::potencia_inv_i( int n, double **A, double *V0, double *V1,
                                    double *eigenVal, double **eigenVec, int p,
                                    bool imprime ){
    double tol, nume, deno, lambda, lambda_old;
	int conta;
	tol = 0.0000000001;
	miVector<double> operVec;
	solver resuelve;
	// Primer lambda_old y V0
	restaEigen( n, p, eigenVec, V0 );
	normaliza( n, V0 );
    resuelve.resuelve_LU( n, A, V1, V0 );
 	nume = operVec.produc_punto( n, V1, V0 );
	deno = operVec.produc_punto( n, V1, V1 );
	operVec.copia_val( n, V0, V1 );
	lambda_old = nume / deno;
	// Primer lambda y V1
	restaEigen( n, p, eigenVec, V0 );
	normaliza( n, V0 );
    resuelve.resuelve_LU( n, A, V1, V0 );
 	nume = operVec.produc_punto( n, V1, V0 );
	deno = operVec.produc_punto( n, V1, V1 );
	lambda = nume / deno;
	conta = 1;
	// Impresión de resultados
	if( imprime ){
		for( int i = 0; i < n * ANCHO_IMP; i++ ) cout << "="; cout << endl;
		cout << "Iter  Eigenvalor  Eigenvector (transpuesto)" << endl;
		imprimeIter( lambda_old, n, V0, 0 );
	}
	// Más aproximaciones
	while( abs( lambda - lambda_old ) > tol ){
		// Actualiza variables
		lambda_old = lambda;
		operVec.copia_val( n, V0, V1 );
		// Quita contribución de eigenvectores anteriores
		restaEigen( n, p, eigenVec, V0 );
		normaliza( n, V0 );
		// Se calcula una nueva aproximación
		resuelve.resuelve_LU( n, A, V1, V0 );
	 	nume = operVec.produc_punto( n, V1, V0 );
		deno = operVec.produc_punto( n, V1, V1 );
		lambda = nume / deno;
		conta++;
		if( imprime ) imprimeIter( lambda_old, n, V0, conta );
	}
	// Guarda el eigenvector y el eigenvalor
	eigenVal[p] = lambda_old;
	for( int i = 0; i < n; i++ ) eigenVec[p][i] = V0[i];
	imprimeRes( lambda_old, n, V0, conta );
	return true;
}
//...........................................................................
bool eigen_vec_val::potencia_inv_i2( int n, double **A, double *V0, double *V1,
                                    double *eigenVal, double **eigenVec, int p,
                                    bool imprime ){
	// Este método en lugar de restar a V0 la contribución de los eigenvectores
	// ya calculados, modifica la matriz<double> [A] de forma que hace 0 los eigenvalores
	// que ya han sido calculados.
    double tol, nume, deno, lambda, lambda_old;
	int conta;
	tol = 0.0000000001;
	miVector<double> operVec;
	solver resuelve;
	// Primer lambda_old y V0
	if( p > 0 ) restaEigen2( n, A, eigenVal, eigenVec, p - 1 );
	normaliza( n, V0 );
    resuelve.resuelve_LU( n, A, V1, V0 );
 	nume = operVec.produc_punto( n, V1, V0 );
	deno = operVec.produc_punto( n, V1, V1 );
	lambda_old = nume / deno;
	normaliza( n, V1 );
	// Primer lambda y V1
	operVec.copia_val( n, V0, V1 );
    resuelve.resuelve_LU( n, A, V1, V0 );
 	nume = operVec.produc_punto( n, V1, V0 );
	deno = operVec.produc_punto( n, V1, V1 );
	lambda = nume / deno;
	normaliza( n, V1 );
	conta = 1;
	// Impresión de resultados
	if( imprime ){
		for( int i = 0; i < n * ANCHO_IMP; i++ ) cout << "="; cout << endl;
		cout << "Iter  Eigenvalor  Eigenvector (transpuesto)" << endl;
		imprimeIter( lambda_old, n, V0, 0 );
	}
	// Más aproximaciones
	while( abs( lambda - lambda_old ) > tol ){
		// Actualiza variables
		lambda_old = lambda;
		operVec.copia_val( n, V0, V1 );
		// Se calcula una nueva aproximación
		resuelve.resuelve_LU( n, A, V1, V0 );
	 	nume = operVec.produc_punto( n, V1, V0 );
		deno = operVec.produc_punto( n, V1, V1 );
		lambda = nume / deno;
		conta++;
		normaliza( n, V1 );
		if( imprime ) imprimeIter( lambda, n, V1, conta );
	}
	// Guarda el eigenvector y el eigenvalor
	eigenVal[p] = lambda;
	for( int i = 0; i < n; i++ ) eigenVec[p][i] = V1[i];
	imprimeRes( lambda, n, V1, conta );
	return true;
}
//...........................................................................
//                      Método de COCIENTE DE RAYLEIGH
//...........................................................................
// Con este método se pueden obtener el valor propio más cercano a un valor
// dado
//...........................................................................
void eigen_vec_val::cociente_Rayleigh( int n, double **A, double h0, 
                                       double *V0 ){
	double h1, a, b;
	int conta = 0;
	solver resuelve;
	miVector<double> operVec, V1( n );
	double tol = 0.0000000001;
	// Factoriza la matriz<double>
	resuelve.factoriza_LU( n, A );
	// Eigenvalor por eigenvector
	operVec.mult_escalar( h0, n, V0 );
	// Resuelve para el nuevo eigenvector
	resuelve.resuelve_LU( n, A, V1.val, V0 );
	// Calcular el nuevo eigenvalor (cociente de Rayleigh)
	a = operVec.produc_punto( n, V1.val, V0 ) / h0;
	b = operVec.produc_punto( n, V1.val, V1.val );
	h1 = a / b * h0;
	cout << "Iter  Eigenvalor  Eigenvector (transpuesto)" << endl;
	imprimeIter( h1, n, V1.val, conta );
	while( abs( abs( a / b ) - 1.0 ) > tol ){
			h0 = h1;
			operVec.copia_val( n, V0, V1.val );
			operVec.mult_escalar( h0, n, V0 );
			resuelve.resuelve_LU( n, A, V1.val, V0 );
			a = operVec.produc_punto( n, V1.val, V0 ) / h0;
			b = operVec.produc_punto( n, V1.val, V1.val );
			h1 = a / b * h0;
			conta++;
			imprimeIter( h1, n, V1.val, conta );
	}
	imprimeRes( h1, n, V1.val, conta );
}
void eigen_vec_val::cociente_Rayleigh2( int n, double **A, double h0, 
                                       double *V0 ){
	matriz<double> operMat, B( n, n );
	miVector<double> operVec, c( n ), V1( n );
	solver resuelve;
	int conta = 0;
	double h1, tol = 0.0000000001;
	// Una primera iteración
	normaliza( n, V0 );
	operMat.copiaVal( n, n, B.val, A );
	for( int i = 0; i < n; i++ ) B.val[i][i] -= h0;
	resuelve.factoriza_LU( n, B.val );
	resuelve.resuelve_LU( n, B.val, V1.val, V0 );
	normaliza( n, V1.val );
	operMat.mult_mat_vec( n, n, A, V1.val, c.val );
	h1 = operVec.produc_punto( n, V1.val, c.val );
	// Impresión
	for( int i = 0; i < n * ANCHO_IMP; i++ ) cout << "-"; cout << endl;
	cout << "Iter  Eigenvalor  Eigenvector (transpuesto)" << endl;
	for( int i = 0; i < n * ANCHO_IMP; i++ ) cout << "-"; cout << endl;
	imprimeIter( h1, n, V1.val, conta );
	while( abs( h1 - h0 ) > tol  ){
		h0 = h1;
		operVec.copia_val( n, V0, V1.val );
		// [B] = [A] - lambda_0 * [I]
		operMat.copiaVal( n, n, B.val, A );
		for( int i = 0; i < n; i++ ) B.val[i][i] -= h0;
		// Factoriza y resuelve para {V1} de [B]{V1}={V0}
		resuelve.factoriza_LU( n, B.val );
		resuelve.resuelve_LU( n, B.val, V1.val, V0 );
		normaliza( n, V1.val );
		// Obtiene el cociente de Rayleigh para aproximar lambda_1
		// lambda_1 = {V1}^T * [A] * {V1}
		operMat.mult_mat_vec( n, n, A, V1.val, c.val );
		h1 = operVec.produc_punto( n, V1.val, c.val );			
		conta++;
		imprimeIter( h1, n, V1.val, conta );
	}
	imprimeRes( h1, n, V1.val, conta );
}
//...........................................................................
//                      Método de JACOBI
//...........................................................................
// Con este método se calculan todos los eigenvalores y eigenvectores de una
// matriz<double> simétrica
//...........................................................................
void eigen_vec_val::jacobi( int n, double **A, bool calcEigenVec, int ini,
                            int fin ){
	double max, tol, theta, aux;
	int p, q, conta = 0;
	matriz<double> operMat, AA, eigenVec;
	tol = 0.000000000000001;
	max = 1.0;
	if( n != 1 ){
		// Por si se desean calcular los eigenvectores
		if( calcEigenVec ){
			AA.reinicia( n, n );
			// Copia de la matriz<double> A
			for( int i = 0; i < n; i++ )
				for( int j = 0; j < n; j++ ) AA.val[i][j] = A[i][j];
			eigenVec.reinicia( abs( fin - ini ) + 1, n );
		}
		while( abs( max ) > tol ){
			max = maxJacobi( n, A, p, q );
			theta = 0.5 * atan2( 2.0 * A[p][q], A[p][p] - A[q][q] );
			// Rota la matriz<double> A
			// B = A * T
			for( int i = 0; i < n; i++ ){
				aux = A[i][p];
				A[i][p] =  A[i][p] * cos( theta ) + A[i][q] * sin( theta );
				A[i][q] = -aux     * sin( theta ) + A[i][q] * cos( theta );
			}
			// C = T^t * B
			for( int i = 0; i < n; i++ ){
				aux = A[p][i];
				A[p][i] =  A[p][i] * cos( theta ) + A[q][i] * sin( theta );
				A[q][i] = -aux     * sin( theta ) + A[q][i] * cos( theta );
			}
			iterJacobi( n, A, conta, theta, max );
			conta++; 
		}
	}
	// Por si se desean calcular los eigenvectores (LOS EIGENVECTORES
	// SE GUARDAN EN RENGLONES)
	if( calcEigenVec )
		eigenVecJacobi( n, A, AA.val, ini - 1, fin - 1, eigenVec.val );
	resultJacobi( n, A, calcEigenVec, ini - 1, fin - 1, eigenVec.val );
}
void eigen_vec_val::jacobi2( int n, double **A, double **PHI ){
	// Encuentra todos los eigenvalores y eigenvectores de [A]. En
	// este caso los eigenvectores están guardados en COLUMNAS.
	double max, tol, theta, aux;
	int p, q, conta = 0;
	matriz<double> operMat;
	tol = 0.000000000000001;
	max = 1.0;
	if( n != 1 ){
		// Inicializa con la identidad la matriz<double> donde se guardarán los
		// eigenvectores
		operMat.diag( n, PHI, 1.0 );
		while( abs( max ) > tol ){
			max = maxJacobi( n, A, p, q );
			theta = 0.5 * atan2( 2.0 * A[p][q], A[p][p] - A[q][q] );
			// Rota la matriz<double> A
			// B = A * T
			for( int i = 0; i < n; i++ ){
				aux = A[i][p];
				A[i][p] =  A[i][p] * cos( theta ) + A[i][q] * sin( theta );
				A[i][q] = -aux     * sin( theta ) + A[i][q] * cos( theta );
				aux = PHI[i][p];
				PHI[i][p] = PHI[i][p] * cos( theta ) + PHI[i][q] * sin( theta );
				PHI[i][q] = -aux      * sin( theta ) + PHI[i][q] * cos( theta );
			}
			// C = T^t * B
			for( int i = 0; i < n; i++ ){
				aux = A[p][i];
				A[p][i] =  A[p][i] * cos( theta ) + A[q][i] * sin( theta );
				A[q][i] = -aux     * sin( theta ) + A[q][i] * cos( theta );
			}
			//iterJacobi( n, A, conta, theta, max );
			conta++; 
		}
	} 
	//resultJacobi2( n, A, PHI );
}
//...........................................................................
void eigen_vec_val::resultJacobi( int n, double **L, bool calcEigenVec,
                                  int ini, int fin, double **eigenVec ){
	for( int i = 0; i < n * ANCHO_IMP; i++ ) cout << "-"; cout << endl;	
	cout << "Eigenvalores Jacobi:" << endl;
	for( int i = 0; i < n; i++ ) cout << L[i][i] << endl;
	if( calcEigenVec ){
		for( int i = 0; i < abs( fin - ini ) + 1; i++ ){
			for( int j = 0; j < n * ANCHO_IMP; j++ ) cout << "-"; cout << endl;	
			cout << "Eigenvalor = " << L[ini+i][ini+i] << endl;
			for( int j = 0; j < n; j++ ) cout << eigenVec[i][j] << "   ";
			cout << endl;
			for( int j = 0; j < n * ANCHO_IMP; j++ ) cout << "-"; cout << endl;
		}
	}
} 
void eigen_vec_val::resultJacobi2( int n, double **L, double **PHI ){
	for( int i = 0; i < n * ANCHO_IMP; i++ ) cout << "-"; cout << endl;	
	cout << "Eigenvalores Jacobi:" << endl;
	for( int i = 0; i < n; i++ ) cout << setprecision(5) << L[i][i] << endl;
	for( int j = 0; j < n; j++ ){
		for( int i = 0; i < n * ANCHO_IMP; i++ ) cout << "-"; cout << endl;	
		cout << "Eigenvalor = " << setprecision(5) << L[j][j] << endl;
		for( int i = 0; i < n; i++ )
			cout << setprecision(3) << PHI[i][j] << "   ";
		cout << endl;
		for( int i = 0; i < n * ANCHO_IMP; i++ ) cout << "-"; cout << endl;
	}
}
//...........................................................................
void eigen_vec_val::eigenVecJacobi( int n, double **L, double **A, int ini, 
                                    int fin, double **eigenVec ){
	// Determina el eigenvector para un eigenvalor dado fijando una
	// entrada del eigenvector
	// [L]    Matriz con los eigenvalores en la diagonal
	// [A]    Matriz original
	// ini-fin Rango de los eigenvalores de los cuales se desea determinar sus
	//         eigenvectores. Va desde 0 hasta n-1
	// [eigenVec] Matriz donde se guardarán los eigenvectores pedidos por
	//            renglón
	double lambda;
	int conta = 0;
	matriz<double> B( n - 1, n - 1 );
	miVector<double> b( n - 1 );
	solver resuelve;
	// Resuelve para los eigenvalores pedidos dentro de un rango
	for( int i = ini; i <= fin; i++ ){
		// Siempre se fija la primer entrada del eigenvector a 1
		for( int j = 1; j < n; j++ ){
			for( int k = 1; k < n; k++ ){
				B.val[j-1][k-1] = A[j][k];
			}
		}	
		lambda = L[i][i];
		for( int j = 0; j < n - 1; j++ ) B.val[j][j] -= lambda;
		// Vector de términos independientes
		for( int j = 0; j < n - 1; j++ ) b.val[j] = -A[j+1][0] * 1.0;
		// Factoriza la matriz<double>
		resuelve.factoriza_LU( n - 1, B.val );
		resuelve.resuelve_LU( n - 1, B.val, eigenVec[conta], b.val );
		for( int j = n - 1; j >= 1; j-- )
			eigenVec[conta][j] = eigenVec[conta][j-1];
		eigenVec[conta][0] = 1.0;
		// Normaliza el miVector<double>
		normaliza( n, eigenVec[conta] );
		conta++;
	}
}
//...........................................................................
double eigen_vec_val::maxJacobi( int n, double **A, int &p, int &q ){
	// Encuentra la posición del máximo de la matriz<double> A que es simétrica 
	// dado i != j. Por lo tanto se busca sólo en una matriz<double> triangular por
	// abajo o encima de la diagonal.
	double max = A[1][0];
	p = 1; q = 0;
	for( int i = 0; i < n; i++ ){
		for( int j = 0; j < i; j++ ){
			if( abs( A[i][j] ) > abs( max ) ){
				max = A[i][j];
				p = i;
				q = j;
			}
		}
	}
	return max;
}
double eigen_vec_val::maxJacobi( int m, int n, double **A, int &p, int &q ){
	// Encuentra la posición del máximo de la matriz<double> A que NO es simétrica
	// que esté fuera de la diagonal definida dentro de una submatriz de n x n
	double max = A[1][0];
	p = 1; q = 0;
	for( int i = 0; i < m; i++ ){
		for( int j = 0; j < n; j++ ){
			if( abs( A[i][j] ) > abs( max ) && i != j ){
				max = A[i][j];
				p = i;
				q = j;
			}
		}
	}
	return max;
}

//...........................................................................
void eigen_vec_val::iterJacobi( int n, double **A, int iter, double theta,
                                double max ){
	// Imprime la matriz<double> A
	for( int i = 0; i < n * ANCHO_IMP; i++ ) cout << "-"; cout << endl;	
	cout << "Iter = " << iter << "   ";
	cout << "Max = " << max << "   ";
	cout << "Theta = " << theta << endl;
	for( int i = 0; i < n; i++ ){
		for( int j = 0; j < n; j++ )
			cout << setw( LONG_NUM ) << setprecision( 4 ) << A[i][j]; 
			//cout << A[i][j] << "    ";
		cout << endl;
	}
}
//...........................................................................
//                     Método de ITERACIÓN EN SUBESPACIO
//...........................................................................
// Con este método se calculan "p" eigenvalores y eigenvectores de una
// matriz<double> (versión clase)
//...........................................................................
void eigen_vec_val::iterSubEspacio( int n, int p, double ** A, double **X ){
	matriz<double> operMat;
	matriz<double> B( p, p ); // Subespacio de [A]
	matriz<double> T( p, p ); // Matriz de rotación que hace a [B] diagonal
	matriz<double> X1( n, p );
	double N1, N2, tol = 1e-15;
	// Inicio De X con vectores independientes
	operMat.iniVal( n, p, X, 0.0 );
	for( int j = 0; j < p; j++ ) X[j][j] = 1.0;
	// Hace una primera iteración
		XT_A_X( p, n, A, X, B.val );
		jacobi2( p, B.val, T.val );
		N2 = normaB( p, B.val );
		potencia_iterSubEsp( p, n, A, X );	
		operMat.mult_mat( n, p, p , X, T.val, X1.val );
	N1 = 2.0 * ( N2 + tol );
	while( abs( N1 - N2 ) > tol ){
		// Actualiza variables
		operMat.copiaVal( n, p, X, X1.val );
		N1 = N2;
		// Realiza la proyección
		XT_A_X( p, n, A, X, B.val );
		// Resuelve el problema de valores y vectores propios m-dimensional
		jacobi2( p, B.val, T.val );
		// Norma de los eigenvalores
		N2 = normaB( p, B.val );
		// Aplica el método de la potencia
		potencia_iterSubEsp( p, n, A, X );
		// Aprox. a la nueva solución
		operMat.mult_mat( n, p, p, X, T.val, X1.val );
	}
	resultIterSubEsp( n, p, B.val, X );
} 
//...........................................................................
// Versión obtenida del archivo prac_autov
void eigen_vec_val::iterSubEspacio2( int n, int p, double **A, double **X ){
	int iter = 10;
	matriz<double> operMat, AUX( n, p ), B( p, p ), T( p, p );
	double N1, N2, tol = 1e-3;
	// Inicia [X]
	operMat.iniVal( n, p, X, 0.0 );
	for( int j = 0; j < p; j++ ) X[j][j] = 1.0;
	// Primer iteración
		// Aplica el método de la potencia
		for( int k = 0; k < iter; k++ ){
			operMat.mult_mat( n, n, p, A, X, AUX.val );
			operMat.copiaVal( n, p, X, AUX.val );
		}
		// Ortonormaliza
		Gram_Schmidt_Modif( n, p, X );
		// Proyección de Rayleigh-Ritz
		XT_A_X( p, n, A, X, B.val );
		// Aplica Jacobi para resolver el problema de valores y vectores
		// propios p-dimensional
		jacobi2( p, B.val, T.val );
		// Nueva solución
		operMat.mult_mat( n, p, p, X, T.val, AUX.val );
		N2 = normaB( p, B.val );
	N1 = 2.0 * ( N2 + tol );
	while( abs( N1 - N2 ) > tol ){
		// Actualizacion de variables
		N1 = N2;
		operMat.copiaVal( n, p, X, AUX.val );
		// Aplica el método de la potencia
		for( int k = 0; k < iter; k++ ){
			operMat.mult_mat( n, n, p, A, X, AUX.val );
			operMat.copiaVal( n, p, X, AUX.val );
		}
		// Ortonormaliza
		Gram_Schmidt_Modif( n, p, X );
		// Proyección de Rayleigh-Ritz
		XT_A_X( p, n, A, X, B.val );
		// Aplica Jacobi para resolver el problema de valores y vectores
		// propios p-dimensional
		jacobi2( p, B.val, T.val );
		// Nueva solución
		operMat.mult_mat( n, p, p, X, T.val, AUX.val );
		N2 = normaB( p, B.val );
	}
	resultIterSubEsp( n, p, B.val, X );
}
//...........................................................................
void eigen_vec_val::Gram_Schmidt_Modif( int n, int p, double **X ){
	// Aplica el método de Gram-Schmiddt modificado para ortonormalizar
	// los vectores-columna de [X]
	miVector<double> PP1( p ), PP2( p );
	double sum, f;
	for( int j = 1; j < p; j++ ){
		// Productos punto {Zj}{Zj} y {Xj}{Zk} con k < j
		sum = 0.0;
		for( int i = 0; i < n; i++ ) sum += X[i][j-1] * X[i][j-1];
		PP1.val[j-1] = sum;
		for( int k = 0; k < j; k++ ){
			sum = 0.0;
			for( int i = 0; i < n; i++ ) sum += X[i][j] * X[i][k]; 
			PP2.val[k] = sum;
		}
		// Resta a {Xj} la proyección de los vectores {Zk} con 0 <= k <= j-1
		for( int k = 0; k < j; k++ ){
			f = PP2.val[k] / PP1.val[k];
			for( int i = 0; i < n; i++ ) X[i][j] -= f * X[i][k];
		}
	}
	// Normaliza
	for( int j = 0; j < p; j++ ){
		sum = 0.0;
		for( int i = 0; i < n; i++ ) sum += X[i][j] * X[i][j];
		sum = sqrt( sum );
		for( int i = 0; i < n; i++ ) X[i][j] /= sum;
	}
}
//...........................................................................
double eigen_vec_val::normaB( int m, double **B ){
	// Obtiene la norma de los eigenvalores de [B]
	double sum = 0.0;
	for( int i = 0; i < m; i++ ) sum += B[i][i] * B[i][i];
	return sqrt( sum );
}
//...........................................................................
void eigen_vec_val::resultIterSubEsp( int n, int m, double **B, double **X ){
	for( int i = 0; i < m * ANCHO_IMP; i++ ) cout << "-"; cout << endl;	
	cout << "Eigenvalores Iteracion Sub-Espacio:" << endl;
	for( int i = 0; i < m; i++ ) cout << setprecision(5) << B[i][i] << endl;
	for( int j = 0; j < m; j++ ){
		for( int i = 0; i < m * ANCHO_IMP; i++ ) cout << "-"; cout << endl;	
		cout << "Eigenvalor = " << setprecision(5) << B[j][j] << endl;
		for( int i = 0; i < n; i++ )
			cout << setprecision(3) << X[i][m-j-1] << "   ";
		cout << endl;
		for( int i = 0; i < m * ANCHO_IMP; i++ ) cout << "-"; cout << endl;
	}
}
//...........................................................................
void eigen_vec_val::XT_A_X( const int m, const int n, double **A, double **X,
                            double **B ){
	// Realiza el producto [X]^T [A] [X] para el método de iteración en
	// subespacios. De aquí se obtiene el subespacio del cual se obtendrán "m"
	// eigenvalores y eigenvectores. [A] tiene que ser simétrica.
	double sum;
	miVector<double> V( n );
	for( int j = 0; j < m; j++ ){
		for( int i = 0; i < m; i++ ){
			for( int p = 0; p < n; p++ ){
				sum = 0.0;
				for( int q = 0; q < n; q++ ){
					sum += A[p][q] * X[q][j];
				}
				V.val[p] = sum;
		    }
			sum = 0.0;
			for( int k = 0; k < n; k++ ) sum += X[k][i] * V.val[k];
			B[i][j] = sum;
		}
	}
}
//...........................................................................
void eigen_vec_val::potencia_iterSubEsp( int m, int n, double **A, double **X ){
	// Aplica una forma del método de la potencia a los eigenvectores [X]
	// de forma que se evite que todos converjan al mismo eigenvector
	// [X] para evitar que todos coincidan al mismo eigenvector.
	// m = número de eigenvectores
	// n = Tamaño de la matriz<double> [A]: n x n
	// [X]: n x m
	int iter = 2;
	matriz<double> operMat, AUX( n, m );
	for( int k = 0; k < iter; k++ ){
		// Se vale que [A][X] = [[A]{X_1},[A]{X_2},...,[A]{X_n}] porque [A]
		// es simétrica
		operMat.mult_mat( n, n, m, A, X, AUX.val );
		operMat.copiaVal( n, m, X, AUX.val );
	}
	// Ortonormaliza
	Gram_Schmidt_Modif( n, m, X );	
}
//...........................................................................
void eigen_vec_val::AXi_iterSubEsp( int n, int p, double **A, double **X,
                                    double *V ){
	// Hace la multiplicación de la matriz<double> [A] (n x n) con el miVector<double>-columna
	// "p" de [X] (n x m). El resultado se guarda en V (n)
	double sum;
    for( int i = 0; i < n; i++ ){
		sum = 0.0;
		for( int j = 0; j < n; j++ ) sum += A[ i ][ j ] * X[j][p];
		V[i] = sum;
	}
}
//...........................................................................
void eigen_vec_val::restaEigen_iterSubEsp( int n, int p, double **X ){
	// Resta la constribución de los "p-1" eigenvectores y eigenvalores
	// anteriores al eigenvector ubicado en "p" en [X]. Los eigenvectores
	// están guardados en columnas.
	// OJO: los eigenvctores de [X] deben estar normalizados.
	double a = 0.0;
	for( int j = 0; j < p; j++ ){
		// Calcula el a_i
		a = 0.0;
        for( int i = 0; i < n; i++ ) a += X[i][j] * X[i][p];
		// Le resta la contribución a v0
		for( int i = 0; i < n; i++ ) X[i][p] -= a * X[i][j];
    }
}
//...........................................................................
void eigen_vec_val::jacobiGral( int n, int m, double **A, double **PHI ){
	// Encuentra todos los eigenvalores y eigenvectores de [A]. En
	// este caso los eigenvectores están guardados en COLUMNAS.
	// FALTA DESCRIBIR m y n PARA LUEGO NO CONFUNDIR!!! Y TAMAÑOS
	// DE LOS ARREGLOS
	double max, tol, theta, aux;
	int p, q, conta = 0;
	matriz<double> operMat;
	tol = 0.000000000000001;
	max = 1.0;
	if( m != 1 ){
		// Inicializa con la identidad la matriz<double> donde se guardarán los
		// eigenvectores
		operMat.iniVal( n, m, PHI, 0.0 );  
		operMat.diag( m, PHI, 1.0 );
		while( abs( max ) > tol ){
			max = maxJacobi( m, A, p, q );
			theta = 0.5 * atan2( 2.0 * A[p][q], A[p][p] - A[q][q] );
			// Rota la matriz<double> A
			// B = A * T
			for( int i = 0; i < m; i++ ){
				aux = A[i][p];
				A[i][p] =  A[i][p] * cos( theta ) + A[i][q] * sin( theta );
				A[i][q] = -aux     * sin( theta ) + A[i][q] * cos( theta );
				aux = PHI[i][p];
				PHI[i][p] = PHI[i][p] * cos( theta ) + PHI[i][q] * sin( theta );
				PHI[i][q] = -aux      * sin( theta ) + PHI[i][q] * cos( theta );
			}
			// C = T^t * B
			for( int i = 0; i < m; i++ ){
				aux = A[p][i];
				A[p][i] =  A[p][i] * cos( theta ) + A[q][i] * sin( theta );
				A[q][i] = -aux     * sin( theta ) + A[q][i] * cos( theta );
			}
			//iterJacobi( n, A, conta, theta, max );
			conta++; 
		}
	} 
	//resultJacobi2( n, A, PHI );
}
