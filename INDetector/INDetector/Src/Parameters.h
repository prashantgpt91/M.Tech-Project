// T3DParameter.h: interface for the T3DParameter class.
//
//////////////////////////////////////////////////////////////////////
/***********************************************************************************/
/* Revision history :
/*
/* Created by Dong-Qing Zhang, Corporate Research Burbank, Thomson Inc. Aug. 2007  */
/*
/* contact : Dong-Qing.Zhang@Thomson.net
/***********************************************************************************/

#if !defined(AFX_T3DPARAMETER_H__AF489F7C_6E82_425D_8246_40706169ED2B__INCLUDED_)
#define AFX_T3DPARAMETER_H__AF489F7C_6E82_425D_8246_40706169ED2B__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

struct stc_FloatPara
{
	float v;
	char  key[100];   // keyword
};

struct stc_StringPara
{
	char  v[1000];    //
	char  key[100];   // keyword
};

class T3DParameter  
{
public:
	stc_FloatPara * m_FloatPara;
	int     m_FloatParaNum;

	stc_StringPara * m_StringPara;
	int     m_StringParaNum;

public:
	int Read(const char * file_name);
	void Copy(T3DParameter * para);
	float Get(const char * key);
	void Set(const char * key, float v);
	void Add(const char * key, float default_v);

public:
	virtual int ScanNum() { return 0; }
	virtual const char * ScanInfo(int scan,int &run_no) { return NULL; }   // get the run number of a scan
	virtual int RunScan(int scan,int run) { return run+1; }   // run the scan

public:
	void print();
	char * GetStr(const char * key);
	void Add(const char * key, const char * default_v);
	void Set(const char * key, const char * v);
	T3DParameter();
	virtual ~T3DParameter();
};

#endif // !defined(AFX_T3DPARAMETER_H__AF489F7C_6E82_425D_8246_40706169ED2B__INCLUDED_)
