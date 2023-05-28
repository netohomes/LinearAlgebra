#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <iomanip>
#include <cstdarg>  // TEMPORAL
//#include <typeinfo> // TEMPORAL
//==================================================================
//                             VECTOR
//==================================================================
//..................................................................
template <class T>
class miVector{
    bool mem; // Indica si val ya tiene memoria
	void pideMemoria( const int n1, const int n2 );
    public:
    // Columnas y renglones
    int  n;
    // Valores
    T *val;

    // Constructores
    miVector( int n );
    miVector( const miVector<T> &V );
    miVector( );

    // Destructor
    ~miVector( );

    // Inicializa el miVector
    void ini( T i_val );

	// Reinicia un miVector
	void reinicia( const int n );
	void limpia( void );

    // Suma y resta de miVectores
    miVector operator +( const miVector<T> &V );
    miVector operator -( const miVector<T> &V );

	// Producto de un vector por un escalar
	void mult_escalar( T );
	void mult_escalar( T, int, T* );

	// Producto contraido de dos vectores
	T produc_punto( int, T*, T* );

    // Asignación
    miVector &operator =( const miVector<T> &V );

	void copia_val( const int, T*, T* );

	// Mínimo y máximo
	T min( const int n, T *V );
	T max( const int n, T *V );

    // Imprime valores
    void print( ) const;
	void print( int num, ... ) const;
	void print( int n, T *val, int num, ... ) const;
};
//..................................................................
// Inicializa el vector a un valor dado
template <class T>
void miVector<T>::ini( T i_val ){
    for( int i = 0; i < n; i++ ) val[ i ] = i_val;
}
//..................................................................
template <class T>
void miVector<T>::reinicia( const int nn){
	// Modifica la memoria alojada
	pideMemoria( n, nn );
	n = nn;
}
template <class T>
void miVector<T>::limpia( void ){
	// Modifica la memoria alojada
	pideMemoria( n, 0 );
	n = 0;
}
//..................................................................
// Suma y resta de vectores
template <class T>
miVector<T> miVector<T>::operator +( const miVector<T> &V ){
    // Verifica que la suma se pueda realizar
    miVector<T> C( n );
    if( n != V.n ){
        //FALTA EXCEPCIÓN
        return C;
    }
    for( int i = 0; i < n; i++ )
            C.val[ i ] = val[ i ] + V.val[ i ];
    return C;
}
template <class T>
miVector<T> miVector<T>::operator -( const miVector<T> &V ){
    // Verifica que la resta se pueda realizar
    miVector<T> C( n );
    if( n != V.n ){
        //FALTA EXCEPCIÓN
        return C;
    }
    for( int i = 0; i < n; i++ )
        C.val[ i ] = val[ i ] - V.val[ i ];
    return C;
}
//..................................................................
template <class T>
void miVector<T>::mult_escalar( T factor ){
	for( int i = 0; i < n; i++ ) val[i] *= factor;
}
template <class T>
void miVector<T>::mult_escalar( T factor, int n, T *V ){
	for( int i = 0; i < n; i++ ) V[i] *= factor;
}
//..................................................................
template <class T>
T miVector<T>::produc_punto( const int n, T *V1, T *V2 ){
	// Calcula el producto contraido de dos miVectores
	T sum = 0.0;
	for( int i = 0; i < n; i++ ) sum += V1[i] * V2[i];
	return sum;
}
//..................................................................
template <class T>
void miVector<T>::copia_val( const int n, T *V1, T *V2 ){
	// Copia los valores de V2 a V1
	for( int i = 0; i < n; i++ ) V1[i] = V2[i];
}
//..................................................................
template <class T>
T miVector<T>::min( const int n, T *V ){
	// Devuelve el valor mínimo del miVector
	T min = V[0];
	for( int i = 0; i < n; i++ ){
		if( V[i] < min ) min = V[i];
	}
	return min;
}
template <class T>
T miVector<T>::max( const int n, T *V ){
	// Devuelve el valor maximo del miVector
	T max = V[0];
	for( int i = 0; i < n; i++ ){
		if( V[i] > max ) max = V[i];
	}
	return max;
}
//..................................................................
// Asignación
template <class T>
miVector<T> &miVector<T>::operator =( const miVector<T> &V ){
    // Aloca memoria para val
	pideMemoria( n, V.n );
	n = V.n;
    // Copia valores
    for( int i = 0; i < n; i++ )
        val[ i ] = V.val[ i ];
    return *this;
}
//..................................................................
// Imprime valores
template <class T>
void miVector<T>::print( ) const{
    std::cout << "MIVECTOR:" << std::endl;
    for( int i = 0; i < n; i++ )
        std::cout << val[ i ] << std::endl;
}
template <class T>
void miVector<T>::print( int num, ... ) const{
	va_list listaParams;
	va_start( listaParams, num );
	for( int i = 0; i < num; i++ )
		std::cout << va_arg( listaParams, char* ) << std::endl;
		//std::cout << typeid( va_arg( listaParams, char* ) ).name() << std::endl;
	va_end( listaParams );
    for( int i = 0; i < n; i++ ) std::cout << val[ i ] << std::endl;
}
template <class T>
void miVector<T>::print( int n, T *val, int num, ... ) const{
	va_list listaParams;
	va_start( listaParams, num );
	for( int i = 0; i < num; i++ )
		std::cout << va_arg( listaParams, char* ) << std::endl;
		//std::cout << typeid( va_arg( listaParams, char* ) ).name() << std::endl;
	va_end( listaParams );
    for( int i = 0; i < n; i++ ) std::cout << val[ i ] << std::endl;
}
//..................................................................
// Constructores
template <class T>
miVector<T>::miVector( int n ): n( n ){
	val = NULL;
	mem = false;
	pideMemoria( 0, n );
	ini( 0.0 );
}
template <class T>
miVector<T>::miVector( const miVector<T> &V ){
	val = NULL;
	mem = false;
	pideMemoria( n, V.n  );
	n = V.n;
    for( int i = 0; i < n; i++ ) val[ i ] = V.val[ i ];
}
template <class T>
miVector<T>::miVector( ): n( 0 ){ val = NULL; mem = false; }
//..................................................................
template <class T>
void miVector<T>::pideMemoria( const int n1, const int n2 ){
	// n1    Cantidad de memoria que se tiene
	// n2    Cantidad de memoria que se quiere tener
	if( n2 > 0 ){
		if( n1 != n2 ){
			if( val != NULL ) delete [] val;
			val = new T [n2];
			mem = true;
		}
	}
	else{
		if( val != NULL ) delete [] val;
		val = NULL;
		mem = false;
	}
}
//..................................................................
// Destructor
template <class T>
miVector<T>::~miVector( ){
    if( val != NULL ) delete [] val;
}
//..................................................................
//==================================================================
///                            MATRIZ
//==================================================================
//..................................................................
template <class T>
class matriz{
    bool mem; // Indica si val ya tiene memoria
    void pideMemoria( const int m1, const int n1, const int m2, const int n2 );

    public:
    // Columnas y renglones
    int  m, n;
    // Valores
    T **val;

    // Constructores
    matriz( int m, int n );
    matriz( const matriz &A );
    matriz( );

    // Destructor
    ~matriz( );

    // Inicializa la matriz a un valor
    void iniVal( T i_val );
	void iniVal( const int, const int, T**, const T );

	// Cambia las dimensiones de la matriz
	void reinicia( const int mm, const int nn );

    // Suma y resta de matrices
    matriz operator +( const matriz &A );
    matriz operator -( const matriz &A );

    // Multiplicación de matrices
    matriz operator *( const matriz &A );
	void mult_mat( const int, const int, const int, T**, T**,
                   T** );
	void mult3mat( const int l, const int m, const int n, const int o,
                   T **A, T **B, T **C, T **D );

    // Multiplicación de una matriz por un miVector
    miVector<T> operator |( const miVector<T> &V );
	void mult_mat_vec( const int, const int, T**,
	                      T*, T* );
	// Transponer
	void transpuesta( int, T** );

	// Matriz diagonal
	void diag( const int, T**, const T );

	// Copiar valores
	void copiaVal( const int, const int, T**, T**);

    // Asignación
    matriz &operator =( const matriz &A );

    // Imprime valores
    void print( ) const;
	void print( const int, const int, T** ) const;
	void print( int, ... ) const;

};
//..................................................................
// Inicializa la matriz a un valor dado
template <class T>
void matriz<T>::iniVal( T i_val ){
    for( int i = 0; i < m; i++ )
        for( int j = 0; j < n; j++ )
            val[ i ][ j ] = i_val;
}
template <class T>
void matriz<T>::iniVal( const int m, const int n, T **A,
                     const T val ){
	for( int i = 0; i < m; i++ )
		for( int j = 0; j < n; j++ )
			A[ i ][ j ] = val;
}
//..................................................................
template <class T>
void matriz<T>::reinicia( const int mm, const int nn ){
    pideMemoria( m , n, mm, nn );
    m = mm;
    n = nn;
}
//..................................................................
// Suma y resta de matrices
template <class T>
matriz<T> matriz<T>::operator +( const matriz<T> &A ){
    // Verifica que la suma se pueda realizar
    matriz<T> C( m, n );
    if( m != A.m || n != A.n ){
        //FALTA EXCEPCIÓN
        return C;
    }
    for( int i = 0; i < m; i++ )
        for( int j = 0; j < n; j++ )
            C.val[ i ][ j ] = val[ i ][ j ] + A.val[ i ][ j ];
    return C;
}
template <class T>
matriz<T> matriz<T>::operator -( const matriz<T> &A ){
    // Verifica que la resta se pueda realizar
    matriz<T> C( m, n );
    if( m != A.m || n != A.n ){
        //FALTA EXCEPCIÓN
        return C;
    }
    for( int i = 0; i < m; i++ )
        for( int j = 0; j < n; j++ )
            C.val[ i ][ j ] = val[ i ][ j ] - A.val[ i ][ j ];
    return C;
}
//..................................................................
// Multiplicación de dos matrices
template <class T>
matriz<T> matriz<T>::operator *( const matriz<T> &A ){
    matriz<T> C( m, A.n );
    // Verifica que es posible hacer la multiplicación
    if( n != A.m ){
        //FALTA EXCEPCIÓN
        return C;
    }
    T sum;
    for( int i = 0; i < m; i++ ){
        for( int j = 0; j < A.n; j++ ){
            sum = 0.0;
            for( int k = 0; k < n; k++ ){
                sum += val[ i ][ k ] * A.val[ k ][ j ];
            }
            C.val[ i ][ j ] = sum;
        }
    }
    return C;
}
template <class T>
void matriz<T>::mult_mat( const int m, const int n, const int p,
                       T **A, T **B, T **C ){
    T sum;
    for( int i = 0; i < m; i++ ){
        for( int j = 0; j < p; j++ ){
            sum = 0.0;
            for( int k = 0; k < n; k++ ){
                sum += A[ i ][ k ] * B[ k ][ j ];
            }
            C[ i ][ j ] = sum;
        }
    }
}
//...........................................................................
// Multiplicación de 3 matrices
template <class T>
void matriz<T>::mult3mat( const int l, const int m, const int n, const int o,
                       T **A, T **B, T **C, T **D ){
	// Realiza la multiplicación de 3 matrices [D] = [A][B][C] sin tener
	// que usar una matriz extra para guardar el producto [A][B] ó [B][C]
	// [D] es de tamaño l x o
    T sum;
	miVector<T> V( m );
	for( int j = 0; j < o; j++ ){
		for( int i = 0; i < l; i++ ){
			for( int p = 0; p < m; p++ ){
				sum = 0.0;
				for( int q = 0; q < n; q++ ){
					sum += B[p][q] * C[q][j];
				}
				V.val[p] = sum;
		    }
			sum = 0.0;
			for( int k = 0; k < m; k++ ) sum += A[i][k] * V.val[k];
			D[i][j] = sum;
		}
	}
}
//...........................................................................
// Multiplicación de una matriz por un miVector
template <class T>
miVector<T> matriz<T>::operator |( const miVector<T> &V ){
    miVector<T> C( V.n );
    // Verifica que es posible hacer la multiplicación
    if( n != V.n ){
        //FALTA EXCEPCIÓN
        return C;
    }
    T sum;
    for( int i = 0; i < m; i++ ){
        sum = 0.0;
        for( int j = 0; j < V.n; j++ )
            sum += val[ i ][ j ] * V.val[ j ];
        C.val[ i ] = sum;
    }
    return C;
}
template <class T>
void matriz<T>::mult_mat_vec( const int m, const int n, T **A,
                           T *v, T *b ){
	// El resultado lo guarda en b
	// m    renglones de la matriz
	// n    columnas de la matriz
    T sum;
    for( int i = 0; i < m; i++ ){
        sum = 0.0;
        for( int j = 0; j < n; j++ )
            sum += A[ i ][ j ] * v[ j ];
        b[ i ] = sum;
    }
}
//..................................................................
template <class T>
void matriz<T>::copiaVal( const int m, const int n, T **A, T **B){
	// Copia los valores de [B] a [A]
	for( int i = 0; i < m; i++ )
		for( int j = 0; j < n; j++ ) A[ i ][ j ] = B[ i ][ j ];
}
//..................................................................
// Asignación
template <class T>
matriz<T> &matriz<T>::operator =( const matriz<T> &A ){
    pideMemoria( m, n, A.m, A.n ); // Debe ir primero
    if( mem ){
        // Copia valores
        for( int i = 0; i < A.m; i++ )
            for( int j = 0; j < A.n; j++ )
                val[ i ][ j ] = A.val[ i ][ j ];
    }
    m = A.m;
    n = A.n;
    return *this;
}
//..................................................................
// Imprime valores
template <class T>
void matriz<T>::print( ) const{
    std::cout << "MATRIZ:" << std::endl;
    if( mem ){
        for( int i = 0; i < m; i++ ){
            for( int j = 0; j < n; j++ )
                std::cout << std::setprecision(5) << val[ i ][ j ] << "\t";
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
}
template <class T>
void matriz<T>::print( int num, ... ) const{
	va_list listaParams;
	va_start( listaParams, num );
	for( int i = 0; i < num; i++ )
		std::cout << va_arg( listaParams, char* ) << std::endl;
		//std::cout << typeid( va_arg( listaParams, char* ) ).name() << std::endl;
	va_end( listaParams );
    if( mem ){
        for( int i = 0; i < m; i++ ){
            for( int j = 0; j < n; j++ )
                std::cout << std::setprecision(5) << val[ i ][ j ] << "\t";
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
}
template <class T>
void matriz<T>::print( const int m, const int n, T **val ) const{
	std::cout << "MATRIZ:" << std::endl;
	for( int i = 0; i < m; i++ ){
		for( int j = 0; j < n; j++ )
			std::cout << val[ i ][ j ] << "\t";
		std::cout << std::endl;
	}
    std::cout << std::endl;
}
//..................................................................
// Transpuesta
template <class T>
void matriz<T>::transpuesta( int n, T **A ){
	// Transpone la matriz [A] gurdándola en el mismo espacio.
	// [A] debe ser cuadrada.
	T aux;
	for( int i = 0; i < n; i++ ){
		for( int j = 0; j <= i - 1; j++ ){
			aux = A[i][j];
			A[i][j] = A[j][i];
			A[j][i] = aux;
		}
	}
}
//..................................................................
template <class T>
void matriz<T>::diag( const int n, T **A, const T val ){
	// Hace la matriz !!cuadrada¡¡ [A] diagonal
	for( int i = 0; i < n; i++ ){
		for( int j = 0; j < n; j++ ){
			if( j == i ) A[i][j] = val;
			else         A[i][j] = 0.0;
		}
	}
}
//..................................................................
// Constructores
template <class T>
matriz<T>::matriz( int m, int n ): m( m ), n( n ){
    val = NULL;
    mem = false;
    pideMemoria( 0, 0, m, n );
    iniVal( 0.0 );
}
template <class T>
matriz<T>::matriz( const matriz<T> &A ){
    m = A.m;
    n = A.n;
    val = NULL;
    mem = false;
    pideMemoria( 0, 0, m, n );
	for( int i = 0; i < m; i++ )
		for( int j = 0; j < n; j++ )
			val[ i ][ j ] = A.val[ i ][ j ];
}
template <class T>
matriz<T>::matriz( ): m( 0 ), n( 0 ){ val = NULL; mem = false; }
//..................................................................
// Destructor
template <class T>
matriz<T>::~matriz( ){
    if( mem ){
        for( int i = 0; i < m; i++ ) delete []  val[ i ];
        delete [] val;
    }
}
//..................................................................
template <class T>
void matriz<T>::pideMemoria( const int m1, const int n1, const int m2,
                          const int n2 ){
    // Pide memoria para un arreglo de m x n
    // m1 y n1 representan la cantidad de memoria que se supone hay
    // en el **val y m2 y n2 son los valores que ahora se quieren
    // tener
    if( m2 > 0 && n2 > 0 ){
        if( m1 != m2 || n1 != n2 ){
            // Libera la que había
			if( val != NULL ){
	            for( int i = 0; i < m1; i++ ) delete []  val[ i ];
    	        delete [] val;
			}
            // Pide nueva
            val = new T* [ m2 ];
                for( int i = 0; i < m2; i++ )
                    val[ i ] = new T [ n2 ];
            mem = true;

   }
    }
    else{
        // Libera la que pudiera haber
        if( val != NULL ){
            for( int i = 0; i < m1; i++ ) delete []  val[ i ];
            delete [] val;
        }
		mem = false;
		val = NULL;
	}

}
//==================================================================
//                      MATRICES RALAS
//==================================================================
template <class T>
class matrizRala{
	int nReng;          // Número de renglones
	miVector <int> nCol;  // Número de columnas de cada renglón
    miVector <bool> memR; // Indica si el renglón "i" tiene memoria
    bool mem;           // Indica si "**val" tiene memoria
    int max;            // Numero máximo de valores distinto de
                        // cero de todos los renglones
 public:
    T **val;       // Valores de la matriz
    int **pos;          // Indica el número de columna que le
                        // corresponde a la entrada de la matriz
	matrizRala( void );
	matrizRala( const int );
	~matrizRala( void );
	void reinicia( const int );
	void reiniciaReng( const int, const int );
    void multVec( const int n, double *v, double *b );
    void imprime( void );
    int r( void );
};
//..................................................................
template <class T>
matrizRala<T>::matrizRala( void ): nReng(0), nCol( ), memR( ){
    mem = false;
    pos = NULL;
	val = NULL;
    max = 0;
}
template <class T>
matrizRala<T>::~matrizRala( void ){
	if( mem ){
        for( int i = 0; i < nReng; i++ ){
            if( memR.val[i] ){
                delete [] val[i];
                delete [] pos[i];
            }
        }
        delete [] val;
        delete [] pos;
	}
}
//..................................................................
template <class T>
matrizRala<T>::matrizRala( const int nr ): nReng( nr ), nCol( nr ),
                                        memR( nr ){
	val = new double* [nr];
	pos = new int* [nr];
	mem = true;
	for( int i = 0; i < nr; i++ ){
        memR.val[i] = false;
        val[i]  = NULL;
        pos[i]  = NULL;
	}
	max = 0;
}
//..................................................................
template <class T>
void matrizRala<T>::reinicia( const int nr ){
    // Reinicia todo
	if( mem ){
	    for( int i = 0; i < nReng; i++ ){
            if( memR.val[i] ){
                delete [] val[i];
                delete [] pos[i];
            }
	    }
		delete [] val;
		delete [] pos;
		nCol.limpia( ); memR.limpia( );
		nReng = 0;
		mem = false;
		val = NULL;
		pos = NULL;
        max = 0;
	}
    // Pide nueva memoria
	if( nr > 0 ){
		nReng = nr;
        nCol.reinicia( nr ); nCol.ini( 0 );
        memR.reinicia( nr ); memR.ini( false );
		val = new double* [nr];
		pos = new int* [nr];
		mem = true;
		for( int i = 0; i < nr; i++ ) pos[i] = NULL;
	}
}
//..................................................................
template <class T>
void matrizRala<T>::reiniciaReng( const int i, const int n ){
    val[i]  = new double [n];
    pos[i]  = new int [n];
    nCol.val[i] = n;
    memR.val[i] = true;
    if( n > max ) max = n;
}
//..................................................................
template <class T>
void matrizRala<T>::multVec( const int n, double *v, double *b ){
    // Multiplica la matriz rala por un Vector {v} y el resultado
    // lo guarda en {b}
    double sum;
    if( max > n ){
		std::cout << "El producto no se puede llevar a cabo";
		std::cout << std::endl;
        return;
    }
    for( int i = 0; i < nReng; i++ ){
        sum = 0.0;
        for( int j = 0; j < nCol.val[i]; j++ ){
            sum += val[i][j] * v[pos[i][j]];
        }
        b[i] = sum;
    }
}
//..................................................................
template <class T>
void matrizRala<T>::imprime( void ){
    for( int i = 0; i < nReng; i++ ){
        if( memR.val[i] )
			for( int j = 0; j < nCol.val[i]; j++ )
				std::cout << val[i][j] << "   ";
		std::cout << std::endl;
    }
}
//..................................................................
template <class T>
int matrizRala<T>::r( void ){  return nReng; }
//..................................................................
//==================================================================
#endif
