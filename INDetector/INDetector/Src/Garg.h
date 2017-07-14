/***********************************************************************************
// Garg.h : I/O functions for Attributed Relational Graph (ARG), and Bag-of-Parts
//            
// 
// Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
//
// Contact : 
//			Dong-Qing Zhang :	 dongqing@gmail.com, 
//			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu
************************************************************************************/

#ifndef GARG_H
#define GARG_H

#include "Graph.h"

// GARG node
class GGNode : public Node
{
public:
	double * feat_vec;
	int		 feat_dim;
};

// GARG edge
class GGEdge : public Edge
{
public:
	double * feat_vec;
	int		 feat_dim;
};

// Generalized Attributed Relational Graph or just Attributed Relational Graph 
class GARG : public Graph
{
public:
	GGNode * m_GGNodes;
	GGEdge * m_GGEdges;
	int m_gFeatDim, m_nFeatDim, m_eFeatDim,m_cFeatDim;

public:
	int GetSpatialLoc(int nodei,int * x,int * y);

	int GetEdgeFeature(double * feat,int ri,int dir);
	int GetNodeFeature(double * feat,int ri);

	static GARG * LoadFromTextFile(const char * garg_dir, const char * garg_name,int h=-1);

	void InitGraph(int node_num,int edge_num,int n_feat_d,int e_feat_d,int c_feat_d=0);

	GARG();
	virtual ~GARG();
};

#define GNODE_PTR(graph,i)  ((GGNode *)((unsigned char *)(graph->node_list)+i*(graph->m_node_size)))
#define GEDGE_PTR(graph,i)  ((GGEdge *)((unsigned char *)(graph->edge_list)+i*(graph->m_edge_size)))

#endif
