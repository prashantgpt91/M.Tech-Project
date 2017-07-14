/***********************************************************************************
// graph.cpp : the base graph class for GARG, IGraph, and MatchMRFGraph
//            
// 
// Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
//
// Contact : 
//			Dong-Qing Zhang :	 dongqing@gmail.com, 
//			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu
************************************************************************************/

#ifndef GRAPH
#define GRAPH

#include "BaseClass.h"
#include "Matrix.h"

// general graph model
class Node
{
public:
	int i;		   // index of node
	short vx,vy;   //visualization coordinates x,y
	int   color;   // the color of the node, i.e. label
};

class Edge
{
public:
	int ni,nj;
	float weight;     // weighed graph
//	short direction;  // 0 = acyclic, 1 = ni->nj, -1 = nj->ni
};

class Clique
{
public:
	bool IsThis(int i,int j,int k);
	int ni[10];   // clique nodes
	int n;		  // clique order
};

class Graph : public MyBaseClass
{
	int _GetComponent(int i,int * row);
	int _GetComponentAff(double *A,int size,bool * node_list,int i,int * row);

public:
// graph node list
	Node * node_list;
	int node_num;
	int m_node_size;  // the data size of node 

// graph edge list
	Edge * edge_list;
	int edge_num;
	int m_edge_size;  // the data size of edge

	Clique * clique_list; // clique list, each one corresponding to k-th order clique
	int clique_num;		// clique number
	int m_clique_size;  // the data size of clique

public:
	void InitGraph(int node_num,int edge_num,int node_size = 0,int edge_size = 0);

	bool Reachable(int ni,int nj);

	int Ypos(int nodei);
	int Xpos(int nodei);

	int SearchClique(Clique * cliq_sc);

	int  FromColorVectorToNodes(PRMatrix * mat);
	int  GetConnectedComponent(int * result,int row,int column,int size_thr=0);
	int  GetCCFromAffinityMatrix(double * A,int size,int * result,int result_size);
	PRMatrix * FromEdgesToAffinityMatrix();   // from edges to affinity matrix

	void SetWeight(int nodei,int nodej,double weight);
	void FormCompleteGraph();

	void ProduceCliqueList(int clq_order,int data_size);
	void WeightFromAffinityMatrix(PRMatrix * mat);

// node, edge, clique functions
	virtual Edge * GetEdgeFromNodes(int ni,int nj,short *dir);   // get the edge from the start and end nodes
	virtual int    GetEdgeFromNodesInt(int ni,int nj,short * dir);
	void    GetNodeFromEdge(int EdgeId,int dir,int * ni,int * nj);
	int     GetNeighborClique(int nodei,int order,int result_n[]);

	virtual void ConstructGraph();
	virtual void ConstructNodes();
	virtual void ConstructEdges();

	virtual short GetNeighborNode(int nodei, int result[]);   
	virtual short GetNeighborNodeAndEdge(int nodei, int result_n[],int result_e[]);
	virtual short GetNeighborEdge(Node * node, int result[]);   
	virtual short GetNeighborEdge(int nodei, int result[]);   

	virtual short GetNeighborNode(Edge * edge, int result[]);   
	virtual short GetNeighborEdge(Edge * edge, int result[]);   

public:
	Graph(int node_num,int edge_num,int node_size = 0,int edge_size = 0);
	Graph() {  edge_list=0;  node_list=0; }; 
	virtual ~Graph();
};

// get the pointers
#define NODE_PTR(i)  ((Node *)((unsigned char *)node_list+i*m_node_size))
#define EDGE_PTR(i)  ((Edge *)((unsigned char *)edge_list+i*m_edge_size))
#define CLIQUE_PTR(i)  ((Clique *)((unsigned char *)clique_list+i*m_clique_size))

#endif
