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

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "time.h"

#include "Graph.h"

void Graph::InitGraph(int n_num,int e_num,int node_size,int edge_size)
{
	node_num = n_num;
	edge_num = e_num;
	clique_num = 0;

	if(edge_size == 0 && node_size == 0)
	{
		node_size = sizeof(Node);
		edge_size = sizeof(Edge);
	}

	edge_list = (Edge * )MYMalloc(edge_size*edge_num);
	memset(edge_list,0,edge_size*edge_num);
	node_list = (Node * )MYMalloc(node_size*node_num);
	memset(node_list,0,node_size*node_num);

	m_node_size = node_size;
	m_edge_size = edge_size;
}

Graph::Graph(int n_num,int e_num,int node_size,int edge_size)
{
	clique_num = 0;
	InitGraph(n_num,e_num,node_size,edge_size);
}

Graph::~Graph()
{
	if(edge_list)
		MYFree(edge_list);
	if(node_list)
		MYFree(node_list);
}

// default graph construction
void Graph::ConstructGraph()
{
	ConstructNodes();
	ConstructEdges();
}

void Graph::ConstructNodes()
{

}

void Graph::ConstructEdges()
{

}

void Graph::GetNodeFromEdge(int EdgeId, int dir,int * ni,int * nj)
{
	if(dir==1)
	{
		*ni = EDGE_PTR(EdgeId)->ni;
		*nj = EDGE_PTR(EdgeId)->nj;
	}
	else
	{
		*ni = EDGE_PTR(EdgeId)->nj;
		*nj = EDGE_PTR(EdgeId)->ni;
	}
}

// default neighborhood search 
short Graph::GetNeighborNode(int nodei, int result_n[])
{
	int i,n=0;
	int ni,nj;

	for(i=0;i<edge_num;i++)
	{
		ni = EDGE_PTR(i)->ni;
		nj = EDGE_PTR(i)->nj;

		if(ni == nodei)
			result_n[n++] = EDGE_PTR(i)->nj;
			
		if(nj == nodei)
			result_n[n++] = EDGE_PTR(i)->ni;
	}

	return n;
}

// default neighborhood search 
short Graph::GetNeighborNodeAndEdge(int nodei, int result_n[],int result_e[])
{
	int i,n=0;

	for(i=0;i<edge_num;i++)
	{
		if(EDGE_PTR(i)->ni == nodei)
		{
			result_n[n] = EDGE_PTR(i)->nj;
			result_e[n++] = i;
		}
			
		if(EDGE_PTR(i)->nj == nodei)
		{
			result_n[n] = EDGE_PTR(i)->ni;
			result_e[n++] = i;
		}
	}

	return n;
}

// get the neighborhood cliques
int Graph::GetNeighborClique(int nodei,int order,int result_n[])
{
	int n,i,j;
	Clique * clq;

	n = 0;
	for(i=0;i<clique_num;i++)
	{
		clq = CLIQUE_PTR(i);
		if(CLIQUE_PTR(i)->n!=order)
			continue;

		for(j=0;j<order;j++)
			if(CLIQUE_PTR(i)->ni[j]==nodei)
			{
				result_n[n++] = i;
				break;
			}
	}

	return n;
}

short Graph::GetNeighborEdge(int nodei, int result_e[])
{
	int i,n=0;

	for(i=0;i<edge_num;i++)
	{
		if(EDGE_PTR(i)->ni == nodei)
		{
			result_e[n++] = i;
		}
			
		if(EDGE_PTR(i)->nj == nodei)
		{
			result_e[n++] = i;
		}
	}

	return n;
}

short Graph::GetNeighborEdge(Node * node, int result[])
{
	return 0;
}

short Graph::GetNeighborNode(Edge * edge, int result[])
{
	return 0;
}

short Graph::GetNeighborEdge(Edge * edge, int result[])
{
	return 0;
}

int Graph::GetEdgeFromNodesInt(int ni,int nj,short * dir)
{
	int i;
	Edge * edge_list_ptr;

	for(i=0;i<edge_num;i++)
	{
		edge_list_ptr = EDGE_PTR(i);

		if((edge_list_ptr->ni == ni && edge_list_ptr->nj == nj))
		{
			if(dir!=NULL)
				*dir = 1;
			return i;
		}
		
		if((edge_list_ptr->nj == ni && edge_list_ptr->ni == nj))
		{
			if(dir!=NULL)
				*dir = -1;
			return i;
		}
	}

	*dir = 0;
	return -1;
}

Edge * Graph::GetEdgeFromNodes(int ni,int nj,short * dir)
{
	int i;
	Edge * edge_list_ptr;

	for(i=0;i<edge_num;i++)
	{
		edge_list_ptr = EDGE_PTR(i);

		if((edge_list_ptr->ni == ni && edge_list_ptr->nj == nj))
		{
			if(dir!=NULL)
				*dir = 1;
			return edge_list_ptr;
		}
		
		if((edge_list_ptr->nj == ni && edge_list_ptr->ni == nj))
		{
			if(dir!=NULL)
				*dir = -1;
			return edge_list_ptr;
		}
	}

	*dir = 0;
	return NULL; 
}

int Graph::Xpos(int nodei)
{
	Node * n1;
	n1 = NODE_PTR(nodei);
	return n1->vx;
}

int Graph::Ypos(int nodei)
{
	Node * n1;
	n1 = NODE_PTR(nodei);
	return n1->vy;
}

void Label2Color(int label,int * r,int *g,int *b)
{
	int color [6][3] = {
		255,0,0,
		0,255,0,
		0,0,255,
		255,255,0,
		0,255,255,
		255,0,255 };

	*r = color[label][0];
	*g = color[label][1];
	*b = color[label][2];
}

// import the weights from affinity matrix
void Graph::WeightFromAffinityMatrix(PRMatrix *mat)
{
	int i,ni,nj;
	Edge * edge_list_ptr;

	for(i=0;i<edge_num;i++)
	{
		edge_list_ptr = EDGE_PTR(i);

		ni = edge_list_ptr->ni;
		nj = edge_list_ptr->nj;
		edge_list_ptr->weight = (float)mat->get(ni+1,nj+1);
	}
	
}

// affinity matrix : the similarity between two nodes
PRMatrix * Graph::FromEdgesToAffinityMatrix()
{
	int i,ni,nj;
	Edge * edge_list_ptr;
	PRMatrix * temp;

	temp = new PRMatrix(node_num,node_num,"graph");

	for(i=0;i<node_num;i++)
	{
		temp->set(i+1,i+1,0);  // diagonal elements set to ZERO
	}

	for(i=0;i<edge_num;i++)
	{
		edge_list_ptr = EDGE_PTR(i);

		ni = edge_list_ptr->ni;
		nj = edge_list_ptr->nj;

		temp->set(ni+1,nj+1,edge_list_ptr->weight);
		temp->set(nj+1,ni+1,edge_list_ptr->weight);
	}

	return temp;
}

int Graph::FromColorVectorToNodes(PRMatrix * mat)
{
	int i;

	for(i=0;i<node_num;i++)
	{
		node_list[i].color = (int)floor(mat->get(i+1,1));
	}

	return 1;
}

// construct complete graph
// the edge number is n(n-1)/2
void Graph::FormCompleteGraph()
{
	int i,j,n=0;
	
	for(i=0;i<node_num;i++)
		for(j=0;j<node_num;j++)
		{
			if(i<=j)
				continue;

			edge_list[n].ni = i;
			edge_list[n].nj = j;
			edge_list[n].weight = 0;
			n++;
		}
}

// set weight of an edge
void Graph::SetWeight(int nodei, int nodej,double weight)
{
	Edge * edge;
	
	edge = GetEdgeFromNodes(nodei,nodej,NULL);
	edge->weight = (float)weight;
}

// reachability test, from ni -> nj, according to weight
bool Graph::Reachable(int ni, int nj)
{
	MyStack * stack;
	int tmp,result[1000],result_e[1000],neib_no=0,k;
	int nt,i;

	for(i=0;i<node_num;i++)
		NODE_PTR(i)->color = 0;

	stack = new MyStack;
	stack->stack_push(ni,0);
	NODE_PTR(i)->color = 1;

	while(stack->stack_pop(&nt,&tmp)!=-1)
	{
		neib_no = GetNeighborNodeAndEdge(nt, result,result_e);
		for(k=0;k<neib_no;k++)
		{
			if(result[k]==nj && EDGE_PTR(result_e[k])->weight > 0.5)
			{
				delete stack;
				return true;
			}
				
			if(NODE_PTR(result[k])->color==0 && EDGE_PTR(result_e[k])->weight > 0.5)  // hasn't been dried
			{
				stack->stack_push(result[k],0);
				NODE_PTR(result[k])->color=1;
			}
		}
	};
	delete stack;

	return false;
}

int Graph::_GetComponent(int i,int * row)
{
	MyStack * stack;
	int n=0,j,tmp,result[1000],result_e[1000],neib_no=0,k;

	stack = new MyStack;
	row[n++] = i;
	stack->stack_push(i,0);

	while(stack->stack_pop(&j,&tmp)!=-1)
	{
		neib_no = GetNeighborNodeAndEdge(j, result,result_e);
		for(k=0;k<neib_no;k++)
		{
			if(node_list[result[k]].color==0 && edge_list[result_e[k]].weight > 0.5)  // hasn't been dried
			{
				row[n++] = result[k];
				stack->stack_push(result[k],0);
				node_list[result[k]].color=1;
			}
		}
	};
	row[n] = -1;
	delete stack;
	return n;
}

// get the connected component of a graph, using edge weight thresholding
int Graph::GetConnectedComponent(int * result,int max_row,int column,int size_thr)
{
	int i,* row,n=0,size;

	row = result;

	for(i=0;i<node_num;i++)
		node_list[i].color = 0;

	for(i=0;i<node_num;i++)
	{
		if(node_list[i].color == 0)
		{
			node_list[i].color = 1;
			size = _GetComponent(i,row);

			if(size>size_thr)  // use size filter
			{
				row += column;
				n++;
			}
			if(n>=max_row)
				break;
		}
	}

	return n;
}

int Graph::_GetComponentAff(double *A,int size,bool * node_list,int i,int * row)
{
	MyStack * stack;
	int n=0,j,tmp,k;

	stack = new MyStack(size);
	row[n++] = i;
	stack->stack_push(i,0);

	while(stack->stack_pop(&j,&tmp)!=-1)
	{
		for(k=0;k<size;k++)
		{
			if(A[j*size+k]>0.0000001 && node_list[k]==false)
			{
				row[n++] = k;
				stack->stack_push(k,0);
				node_list[k] = true;
			}
		}
	};
	row[n] = -1;

	delete stack;
	return n+1;
}

// Get connected component from affinity matrix, instead of from graph
int Graph::GetCCFromAffinityMatrix(double *A, int size, int *result, int result_size)
{
	int i,* row,n=0,ccsize;
	bool * node_list;

	node_list = (bool *)MYMalloc(size*sizeof(bool));
	
	row = result;
	for(i=0;i<size;i++)
		node_list[i] = false;

	for(i=0;i<size;i++)
	{
		if(!node_list[i])
		{
			node_list[i] = true;
			ccsize = _GetComponentAff(A,size,node_list,i,row);
			row += ccsize;
			n++;
		}
	}

	MYFree(node_list);
	return n;
}

// search the clique
int Graph::SearchClique(Clique *cliq_sc)
{
	int i,j,k;

	for(i=0;i<clique_num;i++)
	{
		if(CLIQUE_PTR(i)->n != cliq_sc->n)
			continue;
		
		for(j=0;j<CLIQUE_PTR(i)->n;j++)
		{
			for(k=0;k<cliq_sc->n;k++)
			{
				if(CLIQUE_PTR(i)->ni[j] == cliq_sc->ni[k])
					break;
			}
			if(k==cliq_sc->n)
				break;
		}
		if(j==CLIQUE_PTR(i)->n)
			return i;
	}
	return -1;
}

// exhaustive search of cliques
// clq_order : the clique order
// data_size : the clique data size
void Graph::ProduceCliqueList(int clq_order,int data_size)
{
	int i,j,ni,nj;
	PRMatrix * mat;
	Edge * edge_list_ptr;
	Clique * clq;
	int cliq_limit;
	
	cliq_limit = 100;
	m_clique_size = data_size;

	clique_list = (Clique *)MYMalloc(data_size*cliq_limit);
	clique_num = 0;
	mat = FromEdgesToAffinityMatrix();

// build upon the edges of the graph
// finding 3-clique
	if(clq_order==3)
	{
		for(j=0;j<edge_num;j++)
		{
			edge_list_ptr = EDGE_PTR(j);
			ni = edge_list_ptr->ni;
			nj = edge_list_ptr->nj;
			
			for(i=0;i<node_num;i++)
			{
				if(mat->get(ni+1,i+1)>0 && mat->get(nj+1,i+1)>0)  // form a 3-clique
				{
					CLIQUE_PTR(clique_num)->n = 3; CLIQUE_PTR(clique_num)->ni[0] = ni;
					CLIQUE_PTR(clique_num)->ni[1] = nj; CLIQUE_PTR(clique_num)->ni[2] = i;
					
					if(SearchClique(CLIQUE_PTR(clique_num))==-1)
					{
						clq = CLIQUE_PTR(clique_num);
						clique_num++;
						if(clique_num>=cliq_limit)
						{
							cliq_limit +=100;
							clique_list = (Clique *)realloc(clique_list,data_size*cliq_limit);
						}
					}
				}
			}
		}
	}
	
	delete mat;
}

bool Clique::IsThis(int i, int j, int k)
{
	if(ni[0]==i && ni[1]==j && ni[2]==k) return true;
	if(ni[0]==i && ni[1]==k && ni[2]==j) return true;
	if(ni[0]==j && ni[1]==k && ni[2]==i) return true;
	if(ni[0]==j && ni[1]==i && ni[2]==k) return true;
	if(ni[0]==k && ni[1]==i && ni[2]==j) return true;
	if(ni[0]==k && ni[1]==j && ni[2]==i) return true;

	return false;
}
