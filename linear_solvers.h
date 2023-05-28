#ifndef SOLVERS_H
#define SOLVERS_H
#include <iostream>
//============================================================================
class solver{
    // Metodos
    void imprime( int, double**, double*, double*, const bool, const bool,
                  const bool );
    void valMax( const int, const int, const int, const int, double**,
                 int&, int& );
    int mueveElemento( const int, const int, const int, const int, const int,
                       const int, double**, double*, int* );
	void reacomodaX( const int, double*, int* );
    public:
    solver( );

    void elimGaussiana( int, double**, double*, double* );
	void elimGaussianaTriangInf( int, double**, double*, double* );
    void elimGaussianaPivoteo( int, double**, double*, double* );

	void LU( int, double**, double*, double*, int, std::string );
	void LU_Crout( int, double**, double*, double* );
	void LU_Crout_mejorado( int, double** );
	void escribe_LU( std::string, int, double** );
	void factoriza_LU( int n, double **A );
	void resuelve_LU( int n, double **A, double *x, double *b );

	void LDU( int, double**, double*, double*, std::string );
	void cholesky( int, double** );
	void cholesky_mejorado( int, double** );
	void mult_LD( int, double**, double** );
	void   L_sup( int, double**, double** );
	void escribe_LD( std::string, int, double** );

	void jacobi( int, double**, double*, double* );
	void iniciaVec_jacobi( int, double*, double*, const double );
	bool jacobi_convergencia( int, double*, double*, double );
	void gauss_seidel( int, double**, double*, double* );
	void iniciaVec_gauss_seidel( int, double*, const double );
	double norma_gauss_seidel( int, double* );
	void imprime_iter( int, int, double* );

	void sistema3diagonalSim( int, double*, double*, double*, double* );
	
	void QR( int n, double **A, double *x, double *b );
	void resuelve_QR( int n, double **Q, double **R, double *x, double *b );
	void factoriza_QR( int n, double **A, double **Q, double **R );
	void imprime_QR( int n, double **Q, double **R );

	void gradienteConjugadoMejorado( int n, double **A, double *x, double *b );
	void gradienteConjugadoParalel( int n, double **A, double *x, double *b,
	                                int hilos );

    void resuelveSistTriangSUP( int, double**, double*, double* );
    void resuelveSistTriangINF( int, double**, double*, double* );
};
//============================================================================
#endif
