/****************************************************************************
// DupDos.cpp : main functions for Image Near-Duplicate Detection
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
#include "math.h"
#include "malloc.h"
#include "assert.h"

#include "Configuration.h"

#include "Garg.h"
#include "GraphMatcher.h"

void   GargMatchlearnInit();
void   GargMatchLearn();
double GargMatching(const char * root,const char * file1,const char * file2);
void   BatchMatching(const char * root_dir,const char * list_file,const char * result_file);

int main(int argc, char* argv[])
{
	if(argc==1)
	{
		printf("\n");
		printf("DupDos -0|1 parameter_file -- called by MATLAB for learning\n");
		printf("DupDos -m parameter_file Arg1 Arg2 -- Match Arg1 and Arg2 and print the result\n");
		printf("DupDos -b parameter_file list_file result_file -- Batch matching. List_file is the file listing all Arg pairs for matching. Result_file is a file for saving results\n");

		exit(0);
	}

	ReadConfigFile(argv[2]);

	if(strcmp(argv[1],"-0")==0)   
		GargMatchlearnInit();   // learning init,   called from MATLAB

	if(strcmp(argv[1],"-1")==0)
		GargMatchLearn();		// learning E-step, called from MATLAB

	if(strcmp(argv[1],"-m")==0)
		GargMatching(g_GargRootPath.c_str(),argv[3],argv[4]);   // matching two ARGs, called from MATLAB or script

	if(strcmp(argv[1],"-b")==0)   // batch mode, generate all results
		BatchMatching(g_GargRootPath.c_str(),argv[3],argv[4]);

	return 0;
}

void GargMatchlearnInit()
{
	CAdvString tmp;
	LineDocument ld;
	int i;
	GARG  * garg1, * garg2;

	CAdvString GargDir = g_DataRootPath;
	GargDir += "Gargs/";

	CAdvString ResultFile = g_DataRootPath;
	ResultFile+="Results/";
	ResultFile +="result.txt";

	FILE * fp;
	CAdvString CorresFile1,CorresFile2;

	CAdvString CorresDir;
	CorresDir = g_DataRootPath;
	CorresDir += "Training/X/";

	CAdvString ParaDir;

	MatchMRFGraph * MatchMRF;

	fp = fopen(ResultFile.c_str(),"w");
	fclose(fp);

	MatchMRF = new MatchMRFGraph;

	tmp = g_DataRootPath;
	tmp +="Lists/training_pos.txt";
	ld.Read(tmp.c_str());
	
	ParaDir = g_DataRootPath;
	ParaDir += "Training/Parameters/";

	for(i=0;i<ld.m_TextLineNo;i++)
	{
		garg1 = GARG::LoadFromTextFile(GargDir.c_str(),ld.m_TextLine[i].Word[0],-1);
		garg2 = GARG::LoadFromTextFile(GargDir.c_str(),ld.m_TextLine[i].Word[1],-1);

		CorresFile1.Format("%sxiu_%04d.txt",CorresDir.c_str(),i);
		CorresFile2.Format("%sxijuv_%04d.txt",CorresDir.c_str(),i);

		MatchMRF->SaveBeliefInit(CorresFile1.c_str(),CorresFile2.c_str(),garg1,garg2);

		MatchMRF->FreeMemory();

		delete garg1;
		delete garg2;
	}
	
	delete 	MatchMRF;        // MRF for matching
}

void GargMatchLearn()
{
	CAdvString tmp;
	LineDocument ld;
	int i;
	GARG  * garg1, * garg2;
	double sim_garg;
	double sim_garg_sum=0;

	CAdvString GargDir = g_DataRootPath;
	GargDir+="Gargs/";

	CAdvString ResultFile = g_DataRootPath;
	ResultFile += "Results/";
	ResultFile += "result.txt";

	FILE * fp;
	CAdvString CorresFile1,CorresFile2;

	CAdvString CorresDir;
	CorresDir = g_DataRootPath;
	CorresDir += "Training/X/";

	CAdvString ParaDir;

	MatchMRFGraph * MatchMRF;

	fp = fopen(ResultFile.c_str(),"w");
	fclose(fp);

	MatchMRF = new MatchMRFGraph;

	tmp = g_DataRootPath;
	tmp += "Lists/training_pos.txt";
	ld.Read(tmp.c_str());

	ParaDir = g_DataRootPath;
	ParaDir += "Training/Parameters/";

	MatchMRF->LoadParameters(ParaDir.c_str());

	for(i=0;i<ld.m_TextLineNo;i++)
	{
		garg1 = GARG::LoadFromTextFile(GargDir.c_str(),ld.m_TextLine[i].Word[0],-1);
		garg2 = GARG::LoadFromTextFile(GargDir.c_str(),ld.m_TextLine[i].Word[1],-1);
		
		if(MatchMRF->FormAdjGraph(garg1,garg2,g_AdjNodeLimit)>0)
		{
			MatchMRF->CalculateCompatFunc();

			MatchMRF->Inference(4);

			CorresFile1.Format("%sxiu_%04d.txt",CorresDir.c_str(),i);
			CorresFile2.Format("%sxijuv_%04d.txt",CorresDir.c_str(),i);

			MatchMRF->SaveBelief(CorresFile1.c_str(),CorresFile2.c_str());

			sim_garg = MatchMRF->LogLikelihoodRatio();

			sim_garg_sum += sim_garg;

			fp = fopen(ResultFile.c_str(),"a");
			fprintf(fp,"%s->%s	%f\n",ld.m_TextLine[i].Word[0],
							ld.m_TextLine[i].Word[1],sim_garg);
			fclose(fp);
		}
		else 
			assert(0);

		MatchMRF->FreeMemory();

		delete garg1;
		delete garg2;
	}
	
	sim_garg_sum/=ld.m_TextLineNo;  // average

	delete 	MatchMRF;

	tmp = g_DataRootPath;
	tmp += "results/likelihood.txt";

	fp = fopen(tmp.c_str(),"a");
	fprintf(fp,"%f\n",sim_garg_sum);
	fclose(fp);
}

float EmdSimilarity(GARG *garg1, GARG *garg2);
double GargMatching(const char * root,const char * file1,const char * file2) 
{
	GARG  * garg1, * garg2;
	CAdvString GargDir = root,tmp;
	MatchMRFGraph    * MatchMRF;        // MRF for matching
	double sim_garg=0;

	garg1 = GARG::LoadFromTextFile(GargDir.c_str(),file1,-1);
	garg2 = GARG::LoadFromTextFile(GargDir.c_str(),file2,-1);
	
	if(g_InferenceAlgorithm!=3)  // The new algorithm
	{
		MatchMRF = new MatchMRFGraph;

		tmp = g_DataRootPath;
		tmp += "Training/Parameters/";
		MatchMRF->LoadParameters(tmp.c_str());

		if(MatchMRF->FormAdjGraph(garg1,garg2,g_AdjNodeLimit)>0)
		{
			if(MatchMRF->node_num >2)
			{
				MatchMRF->CalculateCompatFunc();

				MatchMRF->Inference(3);

				sim_garg = MatchMRF->LogLikelihoodRatio();

				printf("%f\n",sim_garg);
			}
			else
			{
				sim_garg = -100.0;

				printf("%f\n",-100.0);
			}
		}
		else 
			assert(0);

		MatchMRF->FreeMemory();
		delete 	MatchMRF;
	}
	else  // Earth Mover's Distance, only for BoP representation
	{
		sim_garg = EmdSimilarity(garg1,garg2);
		printf("%f\n",sim_garg);
	}
	
	delete garg1;
	delete garg2;

	return sim_garg;
}


void BatchMatching(const char * root_dir,const char * list_file,const char * result_file)
{
	CAdvString tmp;
	int i;
	float rst;
	
	LineDocument ld;
	ld.Read(list_file);

	FILE * fp = fopen(result_file,"w");
	fclose(fp);

	for(i=0;i<ld.m_TextLineNo;i++)
	{
		fp = fopen(result_file,"a");

		rst = (float)GargMatching(root_dir,ld.m_TextLine[i].Word[0],ld.m_TextLine[i].Word[1]);
		fprintf(fp,"%s %s %f\r\n",ld.m_TextLine[i].Word[0],ld.m_TextLine[i].Word[1],rst);

		fclose(fp);
	}
}
