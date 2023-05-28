#ifndef EIGEN_VEC_VAL
#define EIGEN_VEC_VAL
//============================================================================
class eigen_vec_val{
	// Contiene los métodos para obtener los eigenvalores y eigenvectores
	// de una matriz
 public:
	eigen_vec_val( );
	~eigen_vec_val( );

	void potencia( int n, double **A, int p, int tipo, bool imprime );
	void restaEigen( int n, int p, double **eigenVec, double *v0 );
	void restaEigen2( int n, double **A, double *eigenVal, double **eigenVec, 
	                  int t );

	bool potencia_i( int n, double **A, double *V0, double *V1, double
	                 *eigenVal, double **eigenVec, int p, bool imprime );
	bool potencia_i2( int n, double **A, double *V0, double *V1, double
	                  *eigenVal, double **eigenVec, int p, bool imprime );

    void potencia_inversa( int n, double **A, int p, int tipo, bool imprime );
    bool potencia_inv_i( int n, double **A, double *V0, double *V1, double
                         *eigenVal, double **eigenVec, int p, bool imprime );
    bool potencia_inv_i2( int n, double **A, double *V0, double *V1, double
                         *eigenVal, double **eigenVec, int p, bool imprime );
	bool normaliza( int n, double *V );

	void cociente_Rayleigh( int n, double **A, double h0, double*V0 );
	void cociente_Rayleigh2( int n, double **A, double h0, double*V0 );

	void jacobi( int n, double **A, bool calcEigenVec = false,int ini = 0,
	             int fin = 0 );
	double maxJacobi( int n, double **A, int &p, int &q );
	void iterJacobi( int n, double **A, int iter, double theta, double max );
	void eigenVecJacobi( int n, double **L, double **A, int ini, int fin,
	                     double **eigenVec );
	void resultJacobi( int n, double **L, bool calcEigenVec, int ini, int fin,
	                   double **eigenVec );
	void resultJacobi( int n, double **L, bool calcEigenVec, int ini, int fin );
	void jacobi2( int n, double **A, double **PHI );
	void resultJacobi2( int n, double **L, double **PHI );

	void iterSubEspacio( int n, int m, double ** A, double **X );
	void iterSubEspacio2( int n, int p, double **A, double **X );
	void XT_A_X( const int m, const int n, double **A, double **X, double **B );
	void potencia_iterSubEsp( int m, int n, double **A, double **X );
	void AXi_iterSubEsp( int n, int p, double **A, double **X, double *V );
	void restaEigen_iterSubEsp( int n, int p, double **X );
	void resultIterSubEsp( int n, int m, double **B, double **X );
	double normaB( int m, double **B );
	void jacobiGral( int m, int n, double **A, double **PHI );
	double maxJacobi( int m, int n, double **A, int &p, int &q );
	void Gram_Schmidt_Modif( int n, int p, double **X );

	void imprimeRes( double lambda, int n, double *v, int iter );
	void imprimeIter( double lambda, int n, double *v, int iter );
};
//============================================================================
#endif
