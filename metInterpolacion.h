#ifndef INTERPOLACION_H
#define INTERPOLACION_H
#include "graficador/grafFunc.h" // Ya tiene a matrix.h
#include "ordenamiento.h"
//============================================================================
class datInterpolacion{
 protected:
	miVector<double> X, Y, operVec;
	ordena < double > ordenar;
 public:
	datInterpolacion( );
	~datInterpolacion( );
	void iniciar( int n, double *X, double *Y );
};
//============================================================================
class polinomios : public datInterpolacion{
	miVector<double> a;
 public:
	polinomios( );
	~polinomios( );
	void calc_coef( void );
	double eval_func( double x );	
};
//============================================================================
class pol_Lagrange : public datInterpolacion{
	miVector<double> a;
 public:
	pol_Lagrange( );
	~pol_Lagrange( );
	void calc_coef( void );
	double eval_func( double x );	
};
//============================================================================
class pol_Newton{
	miVector<double> F, X, Y, operVec;
	double calcF( int n, int ini, int fin, double *F, double *X, double *Y );
 public:
	pol_Newton( );
	~pol_Newton( );
	void coef_Newton( int nPunt, double *X, double *Y );
	double eval_func( double x );
};
//============================================================================
class splines_C0 : public datInterpolacion{
 public:
	splines_C0( );
	~splines_C0( );
	void ordenDat( void );
	double eval_func( double x );
};
class splines_C1 : public datInterpolacion{
	miVector<double> z;
 public:
	splines_C1( );
	~splines_C1( );
	void ordenDat( void );
	void calc_zi( void );
	double eval_func( double x );
};
class splines_C2 : public datInterpolacion{
	miVector<double>  A, B, S, t;
	double C1, C2;
 public:
	splines_C2( );
	~splines_C2( );
	void ordenDat( void );
	void calc_coef( void );
	double eval_func( double x );
};
//============================================================================
class minCuad : public datInterpolacion{
	matriz<double> M;
	miVector<double> C, beta;
	void calc_MC( void );
 public:
	minCuad( );
	~minCuad( );
	void calc_coef( void );
	double eval_func( double x );	
};
//============================================================================
class interpolacion{
	miVector<double> operVec;
 public:
	interpolacion( );
	~interpolacion( );
	// Interpola con el método solicitado

	void interpola( const int metodo, const int nCoord, double *X, double *Y,
	                int argc, char **argv, int calcError, int nE,
	                double *datE );
};
//============================================================================
double evalFunc_Newton( double x );
double evalFunc_splinesC0( double x );
double evalFunc_splinesC1( double x );
double evalFunc_splinesC2( double x );
double evalFunc_pol( double x );
double evalFunc_polLagrange( double x );
double evalFunc_minCuad( double x );
double mideError( int n, double *X, double (*F)( double ),
                  double (*I)( double ) );
double F( double x );
//============================================================================
#endif
