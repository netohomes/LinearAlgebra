#ifndef INTEGRACION_H
#define INTEGRACION_H
//==================================================================
class int_deriv{
    double simpson( double a, double b, int n, double *x );
    double ptosGauss( const int n, const int i );
    double pesosGauss( const int n, const int i );
    double myRandom( double a, double b );
public:
    int_deriv( );
    ~int_deriv( );

    // Integrales triples
    double simpson3(double (*f)( double, double, double ), double x1,
                    double x2, double y1, double y2, double z1,
                    double z2, int n );
    double monteCarlo3( double (*f)( double, double, double ), double x1,
                       double x2, double y1, double y2, double z1,
                       double z2, int n );
    double cuadGauss3( double (*f)( double, double, double ), double x1,
                      double x2, double y1, double y2, double z1,
                      double z2, int n );
    double trapecio(double (*f)( double, double, double ), double x1, 
                    double x2, double y1, double y2, double z1, 
                    double z2 );
    // Derivación
    double difAdelante( double (*f)( double ), double x, int d,
                        int p, double dx );
    double difCentrada( double (*f)( double ), double x, int d,
                        int p, double dx );

};
//==================================================================
#endif
