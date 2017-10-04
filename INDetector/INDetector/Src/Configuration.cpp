/*****************************************************************************
// Configuration.cpp : read configuration file
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

#include "Configuration.h"
#include "BaseClass.h"

CAdvString g_DataRootPath="/";
CAdvString g_GargRootPath = "/";
int g_AdjNodeLimit = 50;    // the adjacent graph node limit
int g_SpatialPrefilterThreshold = 80;  // prefilter threshold with the spatial locations of parts
int g_InferenceAlgorithm = 0; // 0: Gibbs, 1: BP, 2: Mean field BP, 3: Earth Mover's Distance
int g_Representation = 0;     // 0: bag-of-parts, 1: Attributed Relation Graph
int g_Normalization = 0;       // normalization

void ReadConfigFile(const char * file)
{
	FILE * fp;
	CAdvString adv_str;
	CAdvString str,str2;
	char temp[1000],*ret;

	fp = fopen(file,"r");

	if(fp==NULL)
	{
		printf("\nCannot find config file, will use default parameters\n");
		return ;
	}

	while(1)
	{
		ret = fgets(temp,1000,fp);
		if(ret==NULL)
			break;

		adv_str = temp;
		str = adv_str.CutBefore(" \t");

		str.TrimLeft(); str.TrimRight();

		str2 = adv_str.CutBefore("# \n");

		if(str=="DATA_ROOT_PATH")
			g_DataRootPath = str2;

		if(str=="GARG_ROOT_PATH")
			g_GargRootPath = str2;

		if(str=="ADJ_NODE_LIMIT")
			g_AdjNodeLimit = atoi(str2.c_str());

		if(str=="SPATIAL_THRESHOLD")
			g_SpatialPrefilterThreshold = atoi(str2.c_str());

		if(str=="INFERENCE_ALG")
			g_InferenceAlgorithm = atoi(str2.c_str());

		if(str=="REPRESENTATION")
			g_Representation = atoi(str2.c_str());

		if(str=="NORMALIZATION")
			g_Normalization = atoi(str2.c_str());
	};
	
	fclose(fp);
}
