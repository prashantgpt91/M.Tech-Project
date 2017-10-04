/****************************************************************************
// Configuration.h : configuration parameters
// 
// Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
//
// Contact : 
//			Dong-Qing Zhang :	 dongqing@gmail.com, 
//			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu
*****************************************************************************/

#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include "BaseClass.h"

extern CAdvString g_DataRootPath;
extern CAdvString g_GargRootPath;

extern int g_AdjNodeLimit;
extern int g_SpatialPrefilterThreshold;
extern int g_InferenceAlgorithm; // 0: Gibbs, 1: BP, 2: Mean field BP
extern int g_Representation;     // 0: bag-of-parts, 1: Attributed Relation Graph
extern int g_Normalization;       // likelihood offset

void ReadConfigFile(const char * file);

#endif
