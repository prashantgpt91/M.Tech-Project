// T3DParameter.cpp: implementation of the T3DParameter class.
//
//////////////////////////////////////////////////////////////////////
/***********************************************************************************/
/* Revision history :
/*
/* Created by Dong-Qing Zhang, Corporate Research Burbank, Thomson Inc. Aug. 2007  */
/*   Added print(), Add(), Set(), Get(), GetStr(), Copy(), Read() 
/*														by D.-Q.Zhang , Aug. 2007
/* 
/* contact : Dong-Qing.Zhang@Thomson.net
/***********************************************************************************/

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "limits.h"
#include "float.h"
#include "assert.h"

#include "T3DParameter.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

T3DParameter::T3DParameter()
{
	m_FloatParaNum = 0;
	m_FloatPara = NULL;

	m_StringParaNum = 0;
	m_StringPara = NULL;
}

T3DParameter::~T3DParameter()
{
	if(m_FloatPara)
		MYFree(m_FloatPara);

	if(m_StringPara)
		MYFree(m_StringPara);
}

void T3DParameter::print()
{
	int i;

	printf("\nparameter list :\n");

	for(i=0;i<m_StringParaNum;i++)
		printf("%s = %s\n",m_StringPara[i].key, m_StringPara[i].v);

	for(i=0;i<m_FloatParaNum;i++)
		printf("%s = %.2f\n",m_FloatPara[i].key,m_FloatPara[i].v);

	printf("\n");
}

void T3DParameter::Add(const char *key, float default_v)
{
	m_FloatPara = (stc_FloatPara *)realloc(m_FloatPara,(++m_FloatParaNum)*sizeof(stc_FloatPara));
	m_FloatPara[m_FloatParaNum-1].v = default_v;
	strcpy(m_FloatPara[m_FloatParaNum-1].key,key);
}

void T3DParameter::Add(const char *key, const char *default_v)
{
	m_StringPara = (stc_StringPara *)realloc(m_StringPara,(++m_StringParaNum)*sizeof(stc_StringPara));
	strcpy(m_StringPara[m_StringParaNum-1].v, default_v);
	strcpy(m_StringPara[m_StringParaNum-1].key,key);
}

void T3DParameter::Set(const char *key, float v)
{
	for(int i=0;i<m_FloatParaNum;i++)
		if(strcmp(key,m_FloatPara[i].key)==0)
			m_FloatPara[i].v = v;
}

void T3DParameter::Set(const char *key, const char *v)
{
	for(int i=0;i<m_StringParaNum;i++)
		if(strcmp(key,m_StringPara[i].key)==0)
			strcpy(m_StringPara[i].v,v);
}

float T3DParameter::Get(const char *key)
{
	for(int i=0;i<m_FloatParaNum;i++)
		if(strcmp(key,m_FloatPara[i].key)==0)
			return m_FloatPara[i].v;

	return 0;
}

char * T3DParameter::GetStr(const char *key)
{
	for(int i=0;i<m_StringParaNum;i++)
		if(strcmp(key,m_StringPara[i].key)==0)
			return m_StringPara[i].v;

	return 0;
}

void T3DParameter::Copy(T3DParameter *para)
{
	if(m_FloatPara)
		MYFree(m_FloatPara);

	if(m_StringPara)
		MYFree(m_StringPara);

	m_FloatPara = (stc_FloatPara *)MYMalloc(para->m_FloatParaNum*sizeof(stc_FloatPara));
	memcpy(m_FloatPara,para->m_FloatPara,para->m_FloatParaNum*sizeof(stc_FloatPara));

	m_StringPara = (stc_StringPara *)MYMalloc(para->m_StringParaNum*sizeof(stc_StringPara));
	memcpy(m_StringPara,para->m_StringPara,para->m_StringParaNum*sizeof(stc_StringPara));
}

int T3DParameter::Read(const char *file_name)
{
	float v;
	char tmp_str[1000],tmp_str2[1000];

	FILE * fp = fopen(file_name,"r");
	if(fp==NULL)
		return -1;

	while(!feof(fp))
	{
		fscanf(fp,"%s",tmp_str);

		if(tmp_str[0]=='#' || tmp_str[0]=='\n' || tmp_str[0]=='/'
				|| tmp_str[0]==' ' || tmp_str[0]=='\t' || strlen(tmp_str)==0)  // comment line
		{
			fscanf( fp, "%[^\n]s", tmp_str);
			continue;
		}

//		fscanf(fp,"%f",&v);
		fscanf(fp,"%s",&tmp_str2);

		if(tmp_str2[0]>='0' && tmp_str2[0]<='9' || tmp_str2[0]=='-')
			Set(tmp_str, atof(tmp_str2));
		else
			Set(tmp_str, tmp_str2);

		fscanf( fp, "%[^\n]s", tmp_str);
	}

	fclose(fp);

	return 1;
}
