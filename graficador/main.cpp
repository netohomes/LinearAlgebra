#include "graficador.h"
#include <gdk/gdkkeysyms.h>
#include <cstdlib>
data* par;
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
		par->func2D->grafica( par->a, par->b, 100 );
		par->func2D->dibujaPuntos( par->nPunt, par->X, par->Y );
	}
	return FALSE;
}
//............................................................................
gboolean desplaza( GtkWidget *widget, GdkEventKey *event, gpointer aux ){
	// Desplaza y hace zoom sobre la vista de la gráfica
	double desp = 0.1;    // Pixeles
	double p, n;
	p = 1.0;
	n = 1.0;
	switch( event->keyval ){
    	case 's':
			desp *= -1.0;
			par->func2D->desplaza( 0, desp );
      		break;
		case 'w':
			par->func2D->desplaza( 0, desp );
			break;
		case 'd':
			par->func2D->desplaza( 1, desp );
			break;
		case 'a':
			desp *= -1.0;
			par->func2D->desplaza( 1, desp );
			break;
		case '+':
			p *= 1.05;
			par->func2D->zoom( p );
			break;
		case '-':
			p /= 10.0;
			par->func2D->zoom( -p );
			break;
  	}
	return FALSE;
}
//============================================================================
int main( int argc, char *argv[] ){
	int a, b, L, A;
	std::string funcion, var;
	a = atof( argv[1] );
	b = atof( argv[2] );
	funcion = argv[3];
	var = "x";
	// Tamaño de la venta de visulazición
	L = 800; // Pixeles
	A = 600; // Pixeles
	// Guarda los data de la gráfica
	par = new data ( funcion, var, L, A, a, b );
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
	                  G_CALLBACK (desplaza), NULL);

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
	return 0;
}
//============================================================================
