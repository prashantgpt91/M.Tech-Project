/***********************************************************************************
// Matrix.h : matrix class for Graph class
//            
// 
// Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
//
// Contact : 
//			Dong-Qing Zhang :	 dongqing@gmail.com, 
//			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu
************************************************************************************/

#ifndef MATRIX_H
#define MATRIX_H

#include "BaseClass.h"

#define MAX_MATRIX_DIM	10

#define MATRIX_TYPE
#define prMATRIX_DOUBLE "D"
#define prMATRIX_FLOAT "F"
#define prMATRIX_INT "I"
#define prMATRIX_CHAR "C"

// Raw Matrix for low level operations: save, load, naming, structure
class PRRawMatrix : public MyBaseClass
{
public:
	char m_name[256];   // matrix name
	short m_cell_size;   // the size of a matrix cell
	int  m_datasize;
	char m_type[5];       // the data type

	void * m_data;
	short m_dim_arr[MAX_MATRIX_DIM];  // dimension array
	short m_dim_num;  // dimension array size

	void Init(short * dim_arr,int dim_num,const char * name,const char * type);
	void MakeItZero();

// functions 
	PRRawMatrix(short * dim_arr,int dim_num,const char * name,const char * type);
	PRRawMatrix();
	virtual ~PRRawMatrix();
};

class PRMatrix : public PRRawMatrix
{
	int m_row, m_col,m_depth;

	int * i;       // integer
	float *  fl;   // double 
	double * f;    // double 
	unsigned char * c;

	int DataPtr(short * coord,int dim);

public: 

	void GetSize(int * dim,int len);
	void SetName(const char * name);

	double get(int row,int col,int depth);
	double get(short * coord, int dim);
	double get(int row,int col);
	void set(int row,int col,double v);
	void set(int row,int col,int depth,double v);
	void set(short * coord,int dim,double v);

	PRMatrix * getv(short * coord, int dim);

// allocate a matrix, the type can be "D" - double, "F" - float, "I" - integer
	PRMatrix();
// 2D matrix
	PRMatrix(int row,int col, const char * name = "", const char * type = "D");
// 3D matrix
	PRMatrix(int row,int col, int depth, const char * name = "", const char * type = "D");
// any-D matrix 
	PRMatrix(short * dim_arr,int dim_num,const char * name,const char * type);
}; 

double VectorDist(double * v1,double * v2,int dim);
double VectorDot(double * v1,double * v2,int dim);
double VectorLength(double * v1,int dim);
double VectorCos(double *v1,double *v2,int dim);

#endif
