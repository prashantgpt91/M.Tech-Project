/*****************************************************************************
// BaseClass.cpp : classes for I/O, string operations
// 
// Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
//
// Contact : 
//			Dong-Qing Zhang :	 dongqing@gmail.com, 
//			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu
*****************************************************************************/

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "assert.h"
#include "BaseClass.h"

void MyBaseClass::ErrMsg(const char * msg)
{
	printf("%s\n",msg);
}

MyStack::MyStack()
{
	stack = (int *)malloc(sizeof(int)*STACK_MAX_SIZE);
	if(stack==NULL)
	{
		printf("No enough memory\n");
		exit(0);
	}

	max_size = STACK_MAX_SIZE;
	stk_ptr = 0;
}

MyStack::~MyStack()
{
	free(stack);
}

MyStack::MyStack(int size)
{
	stack = (int *)malloc(sizeof(int)*size);
	if(stack==NULL)
	{
		printf("No enough memory\n");
		exit(0);
	}

	max_size = size;
	stk_ptr = 0;
}

int MyStack::stack_push(int x,int y)
{
	if(stk_ptr>=max_size)
	{
//		IBPause("Stack Overflow: heap.cpp");
		return -1;
	}
	stack[stk_ptr++] = x;
	stack[stk_ptr++] = y;
	return 1;
}

int MyStack::stack_pop(int * x,int * y)
{
	if(stk_ptr-1<0)
		return -1;

	*y = stack[--stk_ptr];
	*x = stack[--stk_ptr];
	return 1;
}

char g_Temp[10000];
const char * MyStrcat(const char * cat1,const char * cat2)
{
	strcpy(g_Temp,cat1);
	strcat(g_Temp,cat2);

	return g_Temp;
}

const char * MyStrcatExt(char * a,const char * cat1,const char * cat2)
{
	strcpy(a,cat1);
	strcat(a,cat2);

	return a;
}

const char * MyStrcat(const char * cat1,const char * cat2,const char * cat3)
{
	strcpy(g_Temp,cat1);
	strcat(g_Temp,cat2);
	strcat(g_Temp,cat3);

	return g_Temp;
}

CAdvString::CAdvString()
{
	m_str = NULL;
}

CAdvString::CAdvString(int n)
{
	m_str = (char *)malloc(n+1);
	if(m_str==NULL)
	{
		printf("No enough memory\n");
		exit(0);
	}
}

CAdvString::CAdvString(const char * str)
{
	m_str = NULL;
	*this = str;
}

CAdvString::CAdvString(CAdvString const &a)
{
	m_str = (char *)malloc(strlen(a.m_str)+1);
	if(m_str==NULL)
	{
		printf("No enough memory\n");
		exit(0);
	}

	strcpy(m_str,a.m_str);
}

CAdvString::~CAdvString()
{
	if(m_str)
	{
		free(m_str);
		m_str = NULL;
	}
}

int CAdvString::GetLength()
{
	return strlen(m_str);
}

void CAdvString::TrimLeft(const char *TgtStr)
{
	int p=0,L,i;
	char c, * stmp;

	if(m_str==NULL)
		return ;

	if(TgtStr == NULL)
		TgtStr = " \t";

	L = strlen(TgtStr);

	while((c=m_str[p])!=0)
	{
		for(i=0;i<L;i++)
			if(TgtStr[i]==c)
				break;
		
		if(i!=L)    // find
			p++;
		else
			break;  // couldn't find
	};

	if(p>0)
	{
		stmp = (char *)malloc(strlen(m_str)-p+1);
		if(stmp==NULL)
		{
			printf("No enough memory\n");
			exit(0);
		}

		strcpy(stmp,m_str+p);

		free(m_str);
		m_str = stmp;
	}
}

void CAdvString::TrimRight(const char *TgtStr)
{
	int p,L,i;
	char c;

	if(m_str==NULL)
		return ;

	if(TgtStr == NULL)
		TgtStr = " \t";

	L = strlen(TgtStr);
	p = strlen(m_str)-1;

	while(p!=0)
	{
		c=m_str[p];

		for(i=0;i<L;i++)
			if(TgtStr[i]==c)
				break;
		
		if(i!=L)    // find
			p--;
		else
			break;  // couldn't find
	};

	if(p!=((int)strlen(m_str)-1))
	{
		m_str[p+1]=0;
	}
}

// cut a sub string before ch
CAdvString CAdvString::CutBefore(const char * TgtStr)
{
	int p=0,L,i;
	char c, * stmp;
	CAdvString rstr;

	if(m_str==NULL)
	{
		rstr = "";
		return rstr;
	}

	if(TgtStr == NULL)
	{
		rstr = "";
		return rstr;
	}

	TrimLeft();

	L = strlen(TgtStr);

	while((c=m_str[p])!=0)
	{
		for(i=0;i<L;i++)
			if(TgtStr[i]==c)
				break;
		
		if(i==L)    // couldn't find
			p++;
		else
			break;    // find
	};

	if(p>0)
	{
		stmp = (char *)malloc(strlen(m_str)-p+1);

		strcpy(stmp,m_str+p);

		m_str[p]=0;
		rstr = m_str;

		free(m_str);
		m_str = stmp;

		rstr.TrimLeft();
		rstr.TrimRight();

		return rstr;
	}
	else
	{
		rstr = "";
		return rstr.c_str();
	}
}

const char CAdvString::operator [](int n)
{
	return m_str[n];
}

const bool CAdvString::operator ==(const char * str)
{
	return strcmp(m_str,str)==0;
}

const CAdvString CAdvString::operator +(const char * str)
{
	CAdvString tmp = m_str;

	realloc(tmp.m_str,strlen(m_str)+strlen(str)+1);
	strcpy(tmp.m_str+strlen(m_str),str);

	return tmp;
}

void CAdvString::Format(const char *format, const char *str, int n)
{
	char tmp[1000];

	sprintf(tmp,format,str,n);
	*this = tmp;
}

const void CAdvString::operator +=(const char * str)
{
	m_str = (char *)realloc(m_str,strlen(m_str)+strlen(str)+1);
	if(m_str==NULL)
	{
		printf("No enough memory\n");
		exit(0);
	}

	strcpy(m_str+strlen(m_str),str);
}

const void CAdvString::operator =(const CAdvString & str)
{
	*this = str.m_str;
}

const void CAdvString::operator =(const char * lpsz)
{
	if(m_str)
		free(m_str);

	m_str = (char *)malloc(strlen(lpsz)+1);
	if(m_str==NULL)
	{
		printf("No enough memory\n");
		exit(0);
	}

	strcpy(m_str,lpsz);
}

void SortMachine::Init(int UnitNum)
{
	m_UnitNum = UnitNum;
	m_Arr = (SortEle *)malloc(UnitNum*sizeof(SortEle));
}

SortMachine::SortMachine(int start,int UnitNum)
{
	int i;

	m_UnitNum = UnitNum;
	m_Arr = (SortEle *)malloc(UnitNum*sizeof(SortEle));

	for(i=0;i<UnitNum;i++)
		m_Arr[i].ind = start+i;
}

SortMachine::SortMachine()
{
	m_UnitNum = 0;
	m_Arr = NULL;
}

SortMachine::SortMachine(int UnitNum)
{
	int i;

	m_UnitNum = UnitNum;
	m_Arr = (SortEle *)malloc(UnitNum*sizeof(SortEle));

	for(i=0;i<UnitNum;i++)
		m_Arr[i].ind = i;
}

SortMachine::~SortMachine()
{
	if(m_Arr)
		free(m_Arr);
}

int SM_Compare( const void *arg1, const void *arg2 )
{
	SortEle  * ele1, * ele2;

	ele1 = (SortEle *)arg1;
	ele2 = (SortEle *)arg2;
	
	return ele1->d - ele2->d<0? 1:-1;
}

void SortMachine::Sort()
{
	qsort(m_Arr,m_UnitNum,sizeof(SortEle),SM_Compare);
}

int LineDocument::Read(FILE *fp)
{
	char * ret,line[100000];
	CAdvString str;
	CAdvString str2;
	int n;

	m_TextLineNo = 0;
	while(1)
	{
		ret = fgets(line,100000,fp);
		if(ret==NULL || line[0]=='#')
			break;

		str = line;
		str.TrimLeft();
		str.TrimRight();
		if(str.GetLength()==0)
			continue;

		n = 0;
		while(1)
		{
			if(str.GetLength()==0 || str[0]=='#')
				break;
			str.CutBefore(" \t");
			n++;
		}

		m_TextLine[m_TextLineNo].Word = (LineWord *)MYMalloc(n*sizeof(LineWord));
		m_TextLine[m_TextLineNo].WordNo = n;

		str = line;
		n = 0;
		while(1)
		{
			if(str.GetLength()==0 || str[0]=='#')
				break;
			
			str2 = str.CutBefore(" \t");
			str2.TrimRight(" \t\n");
			str2.TrimLeft(" \t\n");
			strcpy(m_TextLine[m_TextLineNo].Word[n],str2.c_str());

			n++;
		}

		m_TextLineNo++;
	};
	return 1;
}

void LineDocument::FreeMemory()
{
	int i;

	for(i=0;i<m_TextLineNo;i++)
	{
		MYFree(m_TextLine[i].Word);
		m_TextLine[i].Word = NULL;
	}
	m_TextLineNo = 0;

	if(m_TextLine)
	{
		MYFree(m_TextLine);
		m_TextLine = NULL;
	}
}

int LineDocument::Read(const char *filename)
{
	FILE * fp;
	char * ret,line[10000];
	CAdvString str;
	CAdvString str2,err;
	int n;

	fp = fopen(filename,"r");
	if(fp==NULL)
	{
		return 0;
	}

	m_TextLineNo = 0;
	while(1)
	{
		ret = fgets(line,10000,fp);
		if(ret==NULL || line[0]=='#')
			break;
		m_TextLineNo++;
	}
	
	m_TextLine = (TextLine *)MYMalloc(m_TextLineNo*sizeof(TextLine));
	m_TextLineNo = 0;
	rewind(fp);
	while(1)
	{
		ret = fgets(line,10000,fp);
		if(ret==NULL || line[0]=='#')
			break;

		str = line;
		str.TrimLeft();
		str.TrimRight();
		if(str.GetLength()==0)
			continue;

		n = 0;
		while(1)
		{
			if(str.GetLength()==0 || str[0]=='#')
				break;
			str.CutBefore(" \t");
			n++;
		}

		m_TextLine[m_TextLineNo].Word = (LineWord *)MYMalloc(n*sizeof(LineWord));
		m_TextLine[m_TextLineNo].WordNo = n;

		str = line;
		n = 0;
		while(1)
		{
			if(str.GetLength()==0 || str[0]=='#')
				break;
			
			str2 = str.CutBefore(" \t");
			str2.TrimRight(" \t\n\r");
			str2.TrimLeft(" \t\n\r");

			if(str2.GetLength()==0)
				continue;

			assert(strlen(str2.c_str())<200);
			strcpy(m_TextLine[m_TextLineNo].Word[n],str2.c_str());

			n++;
		}

		m_TextLineNo++;
	};

	fclose(fp);
	return 1;
}

LineDocument::LineDocument()
{
	m_TextLine = NULL;
	m_TextLineNo = 0;
}

LineDocument::~LineDocument()
{
	FreeMemory();
}

char * LineDocument::get(int lineid, int columnid)
{
	return m_TextLine[lineid].Word[columnid];
}

int LineDocument::get(int lineid, int columnid,int * ele)
{
	*ele = m_TextLine[lineid].WordInt[columnid];
	return m_TextLine[lineid].WordInt[columnid];
}

int LineDocument::WriteAppend(const char *filename)
{
	FILE * fp;
	int i,j;

	fp = fopen(filename,"a");
	for(i=0;i<m_TextLineNo;i++)
	{
		for(j=0;j<m_TextLine[i].WordNo;j++)
			fprintf(fp,"%s ",m_TextLine[i].Word[j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	return 1;
}

int LineDocument::Find(const char *str, int *LineId, int *WordId)
{
	int i,j;

	for(i=0;i<m_TextLineNo;i++)
	{
		for(j=0;j<m_TextLine[i].WordNo;j++)
		{
			if(strcmp(m_TextLine[i].Word[j],str)==0)
			{
				*LineId = i;  *WordId = j;
				return 1;
			}
		}
	}

	return 0;
}
