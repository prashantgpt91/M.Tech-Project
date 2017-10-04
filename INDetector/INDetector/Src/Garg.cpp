/***********************************************************************************
// Garg.cpp : I/O functions for Attributed Relational Graph (ARG), and Bag-of-Parts
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
#include "time.h"
#include "math.h"
#include "assert.h"

#include "Graph.h"
#include "Garg.h"

GARG::GARG()
{
}

GARG::~GARG()
{
	GGNode * node_list_ptr;
	GGEdge * edge_list_ptr;
	int i;

	for(i=0;i<node_num;i++)
	{
		node_list_ptr = (GGNode *)NODE_PTR(i);
		MYFree(node_list_ptr->feat_vec);
	}

	for(i=0;i<edge_num;i++)
	{
		edge_list_ptr = (GGEdge *)EDGE_PTR(i);
		MYFree(edge_list_ptr->feat_vec);
	}
}

void GARG::InitGraph(int node_num,int edge_num,int n_feat_d,int e_feat_d,int c_feat_d)
{
	GGNode * node_list_ptr;
	GGEdge * edge_list_ptr;
	int i;

	m_nFeatDim = n_feat_d;
	m_eFeatDim = e_feat_d;
	m_cFeatDim = c_feat_d;

	Graph::InitGraph(node_num,edge_num,sizeof(GGNode),sizeof(GGEdge));

	for(i=0;i<node_num;i++)
	{
		node_list_ptr = (GGNode *)NODE_PTR(i);
		node_list_ptr->feat_dim = n_feat_d;
		node_list_ptr->feat_vec = (double*)MYMalloc(sizeof(double)*n_feat_d);
	}

	for(i=0;i<edge_num;i++)
	{
		edge_list_ptr = (GGEdge *)EDGE_PTR(i);
		edge_list_ptr->feat_dim = e_feat_d;
		edge_list_ptr->feat_vec = (double*)MYMalloc(sizeof(double)*e_feat_d);
	}

	m_GGNodes = (GGNode * ) node_list;
	m_GGEdges = (GGEdge * ) edge_list;

}

static LineDocument ld1,ld2,ld_info;

GARG * GARG::LoadFromTextFile(const char *garg_dir, const char *garg_name,int h)
{
	int node_num,edge_num;
	int n_dim,e_dim;
	GARG * output;
	int i,k;
	GGNode * ggnode;
	GGEdge * ggedge;

	ld1.Read(MyStrcat(garg_dir,garg_name,".yi.txt"));  
	ld2.Read(MyStrcat(garg_dir,garg_name,".yij.txt"));
	ld_info.Read(MyStrcat(garg_dir,garg_name,".info.txt"));

	node_num = ld1.m_TextLineNo;
	edge_num = ld2.m_TextLineNo;

	if(node_num==0)
	{
		ld1.FreeMemory();
		ld2.FreeMemory();
		ld_info.FreeMemory();

		return NULL;
	}

	n_dim = ld1.m_TextLine[0].WordNo-1;
	e_dim = 2;

//	assert(edge_num==((node_num)*(node_num-1)/2));

	output = new GARG;

	output->InitGraph(node_num,edge_num,n_dim,e_dim);

// assign nodes
	for(i=0;i<output->node_num;i++)
	{
		ggnode = (GGNode *)GNODE_PTR(output,i);

// Extract one-node feature
		for(k=0;k<n_dim;k++)
			ggnode->feat_vec[k] = atof(ld1.m_TextLine[i].Word[1+k]);
	}

// assign edges
	for(i=0;i<output->edge_num;i++)
	{
		ggedge = (GGEdge *)GEDGE_PTR(output,i);

		ggedge->ni = (int)atoi(ld2.m_TextLine[i].Word[2]);
		ggedge->nj = (int)atoi(ld2.m_TextLine[i].Word[3]);
		ggedge->weight = 1.0;  // just for visualization

		for(k=0;k<e_dim;k++)
			ggedge->feat_vec[k] = atof(ld2.m_TextLine[i].Word[4+k]);
	}

	ld1.FreeMemory();
	ld2.FreeMemory();
	ld_info.FreeMemory();

	return output;	
}

int GARG::GetNodeFeature(double *feat,int ri)
{
	GGNode *nodei;
	int j;

//	nodei = (GGNode *)GNODE_PTR(this,ri);
	nodei = ((GGNode *)((unsigned char *)(node_list)+ri*(m_node_size)));

	for(j=0;j<m_nFeatDim;j++)
		feat[j] = nodei->feat_vec[j];

	return m_nFeatDim;
}

int GARG::GetEdgeFeature(double *feat,int ri,int dir)
{
	GGEdge * edge;

	edge = (GGEdge *)((unsigned char *)(edge_list)+ri*(m_edge_size));

	feat[0] = edge->feat_vec[0]*dir; // dx(r1)
	feat[1] = edge->feat_vec[1]*dir; // dy(r1)

	return m_eFeatDim;
}

int GARG::GetSpatialLoc(int nodei, int *x, int *y)
{
	*x = m_GGNodes[nodei].vx;
	*y = m_GGNodes[nodei].vy;

	return 0;
}
