/***********************************************************************************
// Matrix.cpp : matrix class for Graph class
//            
// 
// Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
//
// Contact : 
//			Dong-Qing Zhang :	 dongqing@gmail.com, 
//			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu
************************************************************************************/

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#include "Matrix.h"

// 2-D matrix
PRRawMatrix::PRRawMatrix(short * dim_arr,int dim_num,const char * name,const char * type)
{
	Init(dim_arr,dim_num,name,type);
}

void PRRawMatrix::Init(short * dim_arr,int dim_num,const char * name,const char * type)
{
	int i;

	m_dim_num = dim_num;
	memcpy(m_dim_arr,dim_arr,dim_num*sizeof(short));

	if(type=="D")
		m_cell_size = sizeof(double);
	else if(type=="F")
		m_cell_size = sizeof(float);
	else if(type=="I")
		m_cell_size = sizeof(int);
	else if(type=="C")
		m_cell_size = sizeof(unsigned char);
	else
		ErrMsg("Error in PRRawMatrix::Init");

	strcpy(m_name,name);
	strcpy(m_type,type);

	m_datasize = 1;
	for(i=0;i<m_dim_num;i++)
		m_datasize*=m_dim_arr[i];
	m_datasize *= m_cell_size;

	if(m_datasize!=0)
	{
		m_data = (void *) MYMalloc(m_datasize);
		memset(m_data,0,m_datasize);
	}
}

void PRRawMatrix::MakeItZero()
{
	if(m_datasize!=0)
		memset(m_data,0,m_datasize);
}

PRRawMatrix::PRRawMatrix()
{
	m_data = NULL;
	m_cell_size = 0;
	m_datasize = 0;
	m_dim_num = 0;
	m_name[0] = 0;
}

PRRawMatrix::~PRRawMatrix()
{
	MYFree(m_data);
	m_data= NULL;
}

void PRMatrix::set(int row,int col,double v)
{
	short coord[2];

	coord[0] = row;
	coord[1] = col;
	set(coord,2,v);
}

double PRMatrix::get(int row,int col)
{
	short coord[2];

	coord[0] = row;
	coord[1] = col;
	return get(coord,2);
}


void PRMatrix::set(int row, int col, int depth, double v)
{
	short coord[3];

	coord[0] = row;
	coord[1] = col;
	coord[2] = depth;

	set(coord,3,v);
}

double PRMatrix::get(int row, int col, int depth)
{
	short coord[3];

	coord[0] = row;
	coord[1] = col;
	coord[2] = depth;

	return get(coord,3);
}

int PRMatrix::DataPtr(short *coord, int dim)
{
	int i,sum;

	sum = 0;
	for(i=dim-1;i>0;i--)
	{
		sum = (sum+(coord[i]-1))*m_dim_arr[i-1];
	}

	sum+= coord[0]-1;

	if(sum*m_cell_size <0 || sum*m_cell_size >= m_datasize)
	{
//		IBPause("error: out of range");
		return -1;
	}

	return sum;
}

// set value: coord -- coordinates
void PRMatrix::set(short *coord, int dim, double v)
{
	int ptr;
	
	ptr = DataPtr(coord,dim);
	if(ptr<0)
		return;

	if(m_type[0]=='D')
		f[ptr] = v;
	else if(m_type[0]=='F')
		fl[ptr] = (float)v;
	else if(m_type[0]=='I')
		i[ptr] = (int)v; 
	else if(m_type[0]=='C')
		c[ptr] = (unsigned char)v;
	else
		ErrMsg("datatype dosn't support");
}

// get vector or matrix
PRMatrix * PRMatrix::getv(short *coord, int dim)
{
	return NULL;
}

// get one value
double PRMatrix::get(short *coord, int dim)
{
	int ptr;
	
	ptr = DataPtr(coord,dim);
	if(ptr<0)
		return 0;

	if(m_type[0]=='D')
		return f[ptr];
	else if(m_type[0]=='F')
		return fl[ptr];
	else if(m_type[0]=='I')
		return i[ptr]; 
	else if(m_type[0]=='C')
		return c[ptr];
	else
		ErrMsg("datatype dosn't support");
	return 0;
}

PRMatrix::PRMatrix(int row,int col, const char * name , const char * type )
{
// init
	short dim_arr[2];
	dim_arr[0] = row;
	dim_arr[1] = col;
	Init(dim_arr,2,name,type);

	fl = (float *)m_data;
	f = (double *)m_data;
	i = (int *)m_data;
	c = (unsigned char *)m_data;
	
	m_row = row; m_col = col;
}

PRMatrix::PRMatrix(int row,int col, int depth,const char * name , const char * type )
{
// init
	short dim_arr[3];
	dim_arr[0] = row;
	dim_arr[1] = col;
	dim_arr[2] = depth;
	Init(dim_arr,3,name,type);
	
	fl = (float *)m_data;
	f = (double *)m_data;
	i = (int *)m_data;
	c = (unsigned char *)m_data;
	
	m_row = row; m_col = col; m_depth = depth;
}

PRMatrix::PRMatrix(short * dim_arr,int dim_num,const char * name,const char * type)
:PRRawMatrix(dim_arr,dim_num,name,type)
{
}

PRMatrix::PRMatrix()
{
}

void PRMatrix::SetName(const char *name)
{
	strcpy(m_name,name);
}

void PRMatrix::GetSize(int *dim, int len)
{
	int i;

	for(i=0;i<m_dim_num;i++)
		dim[i] = m_dim_arr[i];
}


// Vector distance
double VectorDist(double * v1,double * v2,int dim)
{
	int i;
	double sum;

	sum = 0;
	for(i=0;i<dim;i++)
		sum+= (v1[i]-v2[i])*(v1[i]-v2[i]);

	return sqrt(sum);
}

// dot product
double VectorDot(double * v1,double * v2,int dim)
{
	int i;
	double sum;

	sum = 0;
	for(i=0;i<dim;i++)
		sum+= v1[i]*v2[i];

	return sum;
}

double VectorLength(double * v1,int dim)
{
	int i;
	double sum;

	sum = 0;
	for(i=0;i<dim;i++)
		sum+= v1[i]*v1[i];
	return sqrt(sum);
}

// cosine distance, or correlation
double VectorCos(double *v1,double *v2,int dim)
{
	double len1,len2;

	len1 = VectorLength(v1,dim);
	len2 = VectorLength(v2,dim);

	if(len1<1.0e-10 || len2<1.0e-10)
		return 0;

	return VectorDot(v1,v2,dim)/(len1*len2);
}
