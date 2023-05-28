#include "grafFunc.h"
data *par;
//============================================================================
static gboolean draw_cb( GtkWidget *widget, cairo_t *cr, gpointer datos ){
	par = static_cast < data* > ( datos );
	par->func2D = new graficador( par->largo, par->alto,
	                              gtk_widget_get_window ( widget ) );
	// Escoge el modo de evaluación de la función
	if( par->funcion_en_una_cadena )
		// Cuando se pasa a través de una cadena
		par->func2D->funcionEval( par->funcion, par->var );
	else
		// Cuando se pasa a través de un puntero
		par->func2D->funcionEval( par->func );
	// Hace la gráfica
	if( par->funcion_en_una_cadena ) 
		par->func2D->grafica( par->a, par->b, 100 );
	else{
		par->func2D->graficaFija( par->a, par->b, 500 );
		par->func2D->dibujaPuntos( par->nPunt, par->X, par->Y );
	}
	return FALSE;
}
//............................................................................
gboolean desplazaFija( GtkWidget *widget, GdkEventKey *event, gpointer aux ){
	// Desplaza y hace zoom sobre la vista de la gráfica
	double desp = 0.1;    // Pixeles
	double p;
	p = 1.0;
	switch( event->keyval ){
    	case 's':
			desp *= -1.0;
			par->func2D->desplazaFija( 0, desp );
			par->func2D->dibujaPuntos( par->nPunt, par->X, par->Y );
      		break;
		case 'w':
			par->func2D->desplazaFija( 0, desp );
			par->func2D->dibujaPuntos( par->nPunt, par->X, par->Y );
			break;
		case 'd':
			par->func2D->desplazaFija( 1, desp );
			par->func2D->dibujaPuntos( par->nPunt, par->X, par->Y );
			break;
		case 'a':
			desp *= -1.0;
			par->func2D->desplazaFija( 1, desp );
			par->func2D->dibujaPuntos( par->nPunt, par->X, par->Y );
			break;
		case '+':
			p *= 1.05;
			par->func2D->zoomFija( p );
			par->func2D->dibujaPuntos( par->nPunt, par->X, par->Y );
			break;
		case '-':
			p /= 10.0;
			par->func2D->zoomFija( -p );
			par->func2D->dibujaPuntos( par->nPunt, par->X, par->Y );
			break;
  	}
	return FALSE;
}
//============================================================================
//                  MUESTRA EN UNA VENTA UNA GRÁFICA DE CAIRO
//============================================================================
grafFunc::grafFunc( ){ /*par = NULL;*/ }
//............................................................................
void grafFunc::grafica( double a, double b, double (*func)( double ), int L,
                        int A, int argc, char **argv, int n, double *X,
                        double *Y ){

	// Datos para graficar
	par = new data ( L, A, a, b, func, n, X, Y );

	// Inicia gtk
	gtk_init( &argc, &argv );
	// Ventana y área de dibujo
	GtkWidget *window;
	GtkWidget *areaDib;

	// Tìtulo y demensiones de la ventana
	window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
	gtk_window_set_resizable(GTK_WINDOW(window), FALSE);
	gtk_window_set_title(GTK_WINDOW(window), "GRAFICADOR 1.0");
	
	// Para cerrar
	g_signal_connect (window, "destroy", G_CALLBACK (gtk_main_quit), NULL);
	
	// Área de dibujo
	areaDib = gtk_drawing_area_new( );
	gtk_widget_set_size_request ( areaDib, L, A );
	// Llama a la función para dibujar la gráfica
	g_signal_connect( areaDib, "draw", G_CALLBACK(draw_cb), par );

	// EVENTOS IMPLEMENTADOS
	// Para desplazar la vista y hacer zoom
	g_signal_connect( G_OBJECT( window ), "key_press_event",
	                  G_CALLBACK (desplazaFija), NULL);


	// Borde de la ventana
	gtk_container_set_border_width(GTK_CONTAINER(window), 15);
	// Agrega a la ventana lo que se dinujó
	gtk_container_add (GTK_CONTAINER (window), areaDib);
	// Muestra el área de dibujo
	gtk_widget_show( areaDib );
	// Muestra la ventana
	gtk_widget_show( window );

	gtk_main( );

	delete par;
}
//............................................................................
//============================================================================
