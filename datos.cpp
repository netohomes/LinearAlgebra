#include "datos.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
using namespace std;
// Funciones solo necesarias para este .cpp
//............................................................................
void archivoAbierto( ifstream& datos, string archivo  ){
	// Verifica si el archivo fue abierto
	if( !datos.is_open( ) ){
		cerr << "No se pudo abrir el archivo " << archivo << endl;
    	cerr << "PROGRAMA TERMINADO..." << endl;
    	exit( 1 );
	}
}
void archivoAbierto( ofstream& datos, string archivo  ){
	// Verifica si el archivo fue abierto
	if( !datos.is_open( ) ){
		cerr << "No se pudo abrir el archivo " << archivo << endl;
    	cerr << "PROGRAMA TERMINADO..." << endl;
    	exit( 1 );
	}
}
//............................................................................
//============================================================================
// M�TODOS PARA TRABAJAR CON FICHEROS
//============================================================================
//............................................................................
//                         LECTURA DE ARCHIVOS
//............................................................................
datos::datos( void ){ }
//............................................................................
void datos::existeArchivo( string ruta ){
    // Verifica que el archivo de datos exista y si no
    // termina el programa
    ifstream datos( ruta.c_str( ) );
    if( !datos.is_open( ) ){
        cerr << "El archivo " << ruta;
        cerr << " no se pudo abrir" << endl;
        cerr << "PROGRAMA TERMINADO..." << endl;
        exit( 1 );
    }
    else{
        datos.close( );
    }
}
//............................................................................
void datos::leeMatriz( string archivo, matriz<double> &A ){
    // Lee una matriz de un archivo
    string cad;
	int m, n;
	bool fin = false;
    // Verifica que exista el archivo
    existeArchivo( archivo );
    // Abre el flujo al
    ifstream datos( archivo.c_str( ) );
	// Verifica que se haya abierto el archivo
	archivoAbierto( datos, archivo );
    // Lee las dimensiones de la matriz
    datos >> m;
    datos >> n;
    if( m != n ){
        cerr << "Matriz no cuadrada" << endl;
        cerr << "PROGRAMA TERMINADO" << endl;
        exit( 1 );
    }
    A.reinicia( m, n );
	// Lee la matriz
	for( int i = 0; i < m; i++ ){
		for( int j = 0; j < n; j++ ){
			if( datos.eof( ) ){
				fin = true;
				break;
			}
			else{
				datos >> A.val[i][j];
			}
		}
		if( fin ) break;
	}
	//A.print( );
	datos.close( );
}
//............................................................................

void datos::leeMatrizSparse( string archivo, matriz<double> &A ){
    // Lee una matriz de un archivo
    string cad;
	int m, n, i, j;
    // Verifica que exista el archivo
    existeArchivo( archivo );
    // Abre el flujo al 
    ifstream datos( archivo.c_str( ) );
	// Verifica que se haya abierto el archivo
	archivoAbierto( datos, archivo );
    // Lee las dimensiones de la matriz
    datos >> m;
    datos >> n;
    if( m != n ){
        cerr << "Matriz no cuadrada" << endl;
        cerr << "PROGRAMA TERMINADO" << endl;
        exit( 1 );
    }
    A.reinicia( m, n );
	// Lee la matriz
	while( !datos.eof( ) ){
		datos >> i;
		datos >> j;
		datos >> A.val[i][j];
	}
	//A.print( ); 
	datos.close( );
}
//............................................................................

void datos::leeVector( string archivo, miVector<double> &b ){
	// Lee el miVector de t�rminos independientes de un archivo
	int n;
	// Verifica que exista el archivo
    existeArchivo( archivo );
    // Abre el flujo
    ifstream datos( archivo.c_str( ) );
	// Verifica que se haya abierto el archivo
	archivoAbierto( datos, archivo );
    // Lee la dimensi�n del miVector
    datos >> n;
    b.reinicia( n );
	// Lee el miVector
	for( int i = 0; i < n; i++ ){
		if( datos.eof( ) ) break;
		else               datos >> b.val[i];
	}
	//b.print( );
	datos.close( );
}
//............................................................................

void datos::leeVal_Rayleigh( string archivo, int &nCasos, int &n, matriz<double> &V0,
                             miVector<double> &sigma ){
	// Lee un archivo especial para definir los valores iniciales para
	// usar el m�todo Iteraci�n de Cociente de Rayleigh. Estos son las
	// aprozimaciones al eigenvector (phi_0) y al eigenvalor (sigma_0)
	// Tiene extensi�n .ray
	string aux;
	// Verifica que exista el archivo
    existeArchivo( archivo );
    // Abre el flujo
    ifstream datos( archivo.c_str( ) );
	// Verifica que se haya abierto el archivo
	archivoAbierto( datos, archivo );
	// Lee la cantidad de casos a analizar
	datos >> aux;
	datos >> nCasos;
	// Lee el tama�o de los vectores
	datos >> aux;
	datos >> n;
	// Memoria para guardar todos los vectores y los sigma iniciales
	V0.reinicia( nCasos, n );
	sigma.reinicia( nCasos );
	// Lee los vectores y los sigmas
	for( int i = 0; i < nCasos; i++ ){
		if( datos.eof( ) ) break;
		else{
			datos >> aux;
			datos >> sigma.val[i];
			for( int j = 0; j < n; j++ ) datos >> V0.val[i][j];
		}
	}
	datos.close();
}

void datos::leeCoord2D( string archivo, miVector<double> &X, miVector<double> &Y ){
	// Lee una serie de puntos en 2D
	int n;
	// Verifica que exista el archivo
    existeArchivo( archivo );
    // Abre el flujo
    ifstream datos( archivo.c_str( ) );
	// Verifica que se haya abierto el archivo
	archivoAbierto( datos, archivo );
    // Lee la dimensi�n del miVector
    datos >> n;
    X.reinicia( n );
    Y.reinicia( n );
	// Lee las coordenadas
	for( int i = 0; i < n; i++ ){
		if( datos.eof( ) ) break;
		else{
			datos >> X.val[i];
			datos >> Y.val[i];
		}
	}
	//X.print( );
	//Y.print( );
	datos.close( );
}
//............................................................................
//                         ESCRITURA DE ARCHIVOS
//............................................................................

void datos::escribeMatriz( matriz<double> &A, string nombre, string ruta, bool siRuta ){
	// Escribe una matriz en la ruta especificada y con el nombre pasado
	// Prepara el nombre completo del archivo (con su ruta)
	string archivo( ruta );
	if( siRuta == true ) archivo.insert( archivo.length( ), "/" );
	archivo.insert( archivo.length( ), nombre );
	// Abre el flujo al archivo
	ofstream escribe( archivo.c_str( ), std::ofstream::out );
	archivoAbierto( escribe, nombre );
	// Escribe los valores de la matriz
	escribe << A.m << " " << A.n << endl;
	for( int i = 0; i < A.m; i++ ){
	  	for( int j = 0; j < A.n; j++ ){
	  		escribe << A.val[i][j] << " ";
	  	}
	  	escribe << endl;
	}
	// Cierra el flujo
	escribe.close( );
}
//............................................................................

void datos::escribeVector( miVector<double> &V, string nombre, string ruta, bool siRuta ){
	// Escribe una matriz en la ruta especificada y con el nombre pasado
	// Prepara el nombre completo del archivo (con su ruta)
	string archivo( ruta );
	if( siRuta == true ) archivo.insert( archivo.length( ), "/" );
	archivo.insert( archivo.length( ), nombre );
	// Abre el flujo al archivo
	ofstream escribe( archivo.c_str( ), std::ofstream::out );
	archivoAbierto( escribe, nombre );
	// Escribe los valores de la matriz
	escribe << V.n << endl;
	for( int i = 0; i < V.n; i++ ) escribe << V.val[i] << endl;
	// Cierra el flujo
	escribe.close( );
}
//............................................................................

void datos::escribeCoord2D( miVector<double> &X, miVector<double> &Y, string nombre, string ruta,
                            bool siRuta ){
	// Prepara el nombre completo del archivo (con su ruta)
	string archivo( ruta );
	if( siRuta == true ) archivo.insert( archivo.length( ), "/" );
	archivo.insert( archivo.length( ), nombre );
	// Abre el flujo al archivo
	ofstream escribe( archivo.c_str( ), std::ofstream::out );
	archivoAbierto( escribe, nombre );
	// Escribe los valores
	for( int i = 0; i < X.n; i++ ){
	    escribe << X.val[i] << " " << Y.val[i] <<  endl;
	}
	// Cierra el flujo
	escribe.close( );
}
//............................................................................
// void datos::escribeCad( ofstream salida, string cad ){
// 	// Escribe una cadena en el archivo especificado
//   	archivoAbierto( salida );
// 	salida << cad;
// }
//............................................................................
