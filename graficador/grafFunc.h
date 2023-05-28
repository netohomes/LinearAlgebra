#ifndef GRAF_FUNC_H
#define GRAF_FUNC_H
#include "graficador.h"
//============================================================================
class grafFunc{
	//data *par;
 public:
	grafFunc( );
	grafFunc( double a, double b, double (*func)( double ), int argc,
              char **argv, int L, int A );
	void grafica( double a, double b, double (*func)( double ), int L,
                  int A, int argc, char **argv, int n, double *X, double *Y );
};
//============================================================================
gboolean desplaza( GtkWidget *widget, GdkEventKey *event, gpointer aux );
#endif
