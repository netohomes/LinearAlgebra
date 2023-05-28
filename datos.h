#ifndef datos_h
#define datos_h
#include "matrix.h"
//============================================================================
class datos{
	std::string ruta;
	std::string archivo;
public:
    datos( );
    void existeArchivo( std::string );
    void leeMatriz( std::string, matriz<double>& );
	
	void leeMatrizSparse( std::string archivo, matriz<double> &A );
	
    void leeVector( std::string, miVector<double>& );
	
	void leeCoord2D( std::string, miVector<double>&, miVector<double>& );
	
	void leeVal_Rayleigh( std::string, int&, int&, matriz<double>&, miVector<double>& );
	
	void escribeMatriz( matriz<double>&, std::string, std::string = "", bool = false );
	
	void escribeVector( miVector<double>&, std::string, std::string = "", bool = false );
	
	void escribeCoord2D( miVector<double> &X, miVector<double> &Y, std::string nombre, std::string ruta,
                         bool siRuta );
};
//============================================================================
#endif
