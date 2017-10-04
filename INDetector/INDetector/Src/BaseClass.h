/*****************************************************************************
// BaseClass.h : classes for I/O, string operations
// 
// Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
//
// Contact : 
//			Dong-Qing Zhang :	 dongqing@gmail.com, 
//			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu
*****************************************************************************/

#ifndef MYBASECLASS
#define MYBASECLASS

class MyBaseClass
{
public:
	void ErrMsg(const char * msg);
};

typedef char LineWord[200];
typedef struct tag_TextLine
{
	LineWord  * Word;
	int		  * WordInt;  // corresponding integer
	int			WordNo;  
} TextLine;

// for reading formatted text files
class LineDocument : public MyBaseClass
{
public:
	TextLine * m_TextLine;
	int		m_TextLineNo; 

public:
	void FreeMemory();

	int Find(const char * str,int * LineId,int * WordId);
	int WriteAppend(const char * filename);
	int Read(FILE * fp);
	int Read(const char * filename);

	int get(int lineid,int columnid,int * ele);
	char * get(int lineid,int columnid);

    LineDocument();
	~LineDocument();
};

typedef struct tag_SortEle
{
	int ind;   // index 1
	int ind2;  // index 2
	double d;
} SortEle;

class SortMachine
{
public:
	int    m_UnitNum;
	SortEle * m_Arr;

public:
	void Init(int UnitNum);
	void Sort();

	SortMachine();
	SortMachine(int UnitNum);
	SortMachine(int start,int UnitNum);
	~SortMachine();
};

#define STACK_MAX_SIZE 256*256*2
class MyStack
{
	int max_size;
	int * stack;
	int stk_ptr;

public:
	MyStack();
	~MyStack();
	MyStack(int size);

	int stack_push(int x,int y);
	int stack_pop(int * x,int * y);
};

// advanced string
class CAdvString
{
public:
	char * m_str;

public:
	char * c_str() { return m_str; }

	int GetLength();
	void TrimRight(const char * TgtStr = NULL);
	void TrimLeft(const char * TgtStr = NULL);
	CAdvString CutBefore(const char * ch);
	void Format(const char * format,const char * str,int n);

	const void operator =(const char *lpsz);
	const void operator =(const CAdvString & str);

	const char operator [](int n);
	const bool operator ==(const char *lpsz);
	const CAdvString operator +(const char *lpsz);
	const void operator +=(const char *lpsz);
	
	CAdvString();
	CAdvString(int n);
	CAdvString(const char * str);
	CAdvString(CAdvString const &a);

	virtual ~CAdvString();
};

const char * MyStrcat(const char * cat1,const char * cat2);
const char * MyStrcat(const char * cat1,const char * cat2,const char * cat3);
const char * MyStrcatExt(char * a,const char * cat1,const char * cat2);

#define MYMalloc(size)	 malloc(size)
#define MYFree(pointer)  free(pointer)

#endif
