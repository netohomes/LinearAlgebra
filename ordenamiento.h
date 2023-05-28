#ifndef ORDENA_H
#define ORDENA_H
//#include <cstdlib>
//==================================================================
template < class T >
class ordena{
	void une( int ini_1, int fin_1, int ini_2, int fin_2, T *V,
	          bool otro, T *X );
public:
	ordena( );
	~ordena( );
	void mergeSort( int ini, int fin, T *V, bool otro = false,
	                T *B = NULL );
	void burbuja( int n, T *V,  bool otro = false, T *B = NULL );
};
//==================================================================
//..................................................................
template < class T > ordena<T>::ordena( ){}
template < class T > ordena<T>::~ordena( ){}
//..................................................................
template < class T >
void ordena<T>::mergeSort( int ini, int fin, T *V, bool otro,
                            T *B ){
	// Emplea el método de Merge Sort para ordenar una serie de
	// datos
	// ini < fin
	// Ordena con respecto a V
	T aux;
	int ini_1, fin_1, ini_2, fin_2;
	if( fin - ini < 1 ) return;
	if( fin - ini == 1 ){
		if( V[ini] > V[fin] ){
			aux = V[ini];
			V[ini] = V[fin];
			V[fin] = aux;
			if( otro ){
				aux = B[ini];
				B[ini] = B[fin];
				B[fin] = aux;
			}
		}
		return;
	}	
	// Determina el inicio y fin de los 2 subarreglos
	ini_1 = ini;
	fin_1 = ( ini + fin - ( fin - ini ) % 2 ) / 2;
	ini_2 = fin_1 + 1;
	fin_2 = fin;
	// Vuelve a llamar la función para cada división
	mergeSort( ini_1, fin_1, V, otro, B );
	mergeSort( ini_2, fin_2, V, otro, B );
	une( ini_1, fin_1, ini_2, fin_2, V, otro, B );
}
//..................................................................
template < class T >
void ordena<T>::une( int ini_1, int fin_1, int ini_2, int fin_2,
                     T *V, bool otro, T *X ){
	// Hace la unión ordenada de las divisiones hechas por
	// merge sort
	int j, k, i;
	/*if( fin_1 != ini_2 - 1 ){
		std::cout << endl << "Error en el ordenamiendo" << endl;
		exit( 1 );
	}*/
	// Vector auxiliares
	T *Y;
	T *B = new T [ fin_1 - ini_1 + 1 ];
	if( otro ) Y = new T [ fin_1 - ini_1 + 1 ];
	for( int i = 0; i <= fin_1 - ini_1; i++ ) B[i] = V[i+ini_1];
	if( otro )
		for( int i = 0; i <= fin_1 - ini_1; i++ ) Y[i] = X[i+ini_1];
	i = j = k = 0;
	while( ( j <= fin_1 - ini_1 ) && ( k <= fin_2 - ini_2 ) ){
		if( V[k+ini_2] < B[j] ){
			V[i+ini_1] = V[k+ini_2];
			if( otro ) X[i+ini_1] = X[k+ini_2];
			k++;
			i++;
		}
		else{
			V[i+ini_1] = B[j];
			if( otro ) X[i+ini_1] = Y[j];
			j++;
			i++;
		}
	}
	while( j <= fin_1 - ini_1 ){
		V[i+ini_1] = B[j];
		if( otro ) X[i+ini_1] = Y[j];
		i++;
		j++;
	}
	// Libera la memoria
	delete [] B;
	if( otro ) delete [] Y;
}
//..................................................................
template < class T >
void ordena<T>::burbuja( int n, T *V, bool otro, T *B ){
	// Ordena con respecto a V
	T aux;
	for( int i = 0; i < n; i++ ){
		for( int j = i + 1; j < n; j++ ){
			if( V[j] < V[i] ){
				aux  = V[i];
				V[i] = V[j];
				V[j] = aux;
				if( otro ){
					aux  = B[i];
					B[i] = B[j];
					B[j] = aux;
				}
			}
		}
	}
}
//..................................................................
//==================================================================
#endif
