/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: antsDijkstrasPathGraphFilter.hxx,v $
  Language:  C++
  Date:
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __antsDijkstrasPathGraphFilter_hxx
#define __antsDijkstrasPathGraphFilter_hxx

#include "antsDijkstrasPathGraphFilter.h"
#include "itkGraphToGraphFilter.h"
#include "vnl/vnl_numeric_traits.h"
#include <algorithm>
#include <utility>

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

using namespace boost;

namespace itk {
namespace ants {

template <class TGraph>
DijkstrasPathGraphFilter<TGraph>
::DijkstrasPathGraphFilter()
{
  this->m_SourceNodes.clear();
  this->m_TargetNodes.clear();
}

/** Generate the data */
template <class TGraph>
void
DijkstrasPathGraphFilter<TGraph>
::GenerateData()
{
  this->GenerateDistances();
}

template <class TGraph>
void
DijkstrasPathGraphFilter<TGraph>
::GenerateDistances()
{
  this->Initialize();


	// boost code
	//typedef adjacency_list < vecS, vecS, undirectedS, no_property, property < edge_weight_t, double > > BoostGraphType;
	//typedef graph_traits< BoostGraphType >::vertex_descriptor BoostVertexDescriptor;
	//typedef graph_traits< BoostGraphType >::edge_descriptor BoostEdgeDescriptor;
	//typedef std::pair<int, int> BoostEdgeType;

	BoostGraphType bGraph( m_Output->GetTotalNumberOfNodes() );
	std::vector<BoostEdgeType> edgeVec( m_Output->GetTotalNumberOfEdges() );

  property_map<BoostGraphType, edge_weight_t>::type weightmap = get(edge_weight, bGraph);

  for (unsigned long i=0; i<m_Output->GetTotalNumberOfEdges(); i++)
    {
		BoostEdgeDescriptor e; 
		bool inserted;
		
		EdgePointerType graphEdge = m_Output->GetEdgePointer( i );
		int curNodeNum = graphEdge->SourceIdentifier;
		int toNodeNum = graphEdge->TargetIdentifier;
		
		tie(e, inserted) = add_edge(curNodeNum, toNodeNum, bGraph);
		
		if (graphEdge->Weight < 0)
		  {
		  std::cout << "Negative edge weight detected.. exiting" << std::endl;
		  return;
      }
		
		weightmap[e] = graphEdge->Weight;
    }

  int startNodeNum = this->m_SourceNodes[0];    
  this->m_PreviousVertexList.reserve( num_vertices(bGraph) );
  
	std::vector<int> distance(num_vertices(bGraph));
	BoostVertexDescriptor s = vertex(startNodeNum, bGraph);
	typename property_map<BoostGraphType, vertex_index_t>::type indexmap = get(vertex_index, bGraph);
	
  std::cout << "Start boost..." << std::flush;
	dijkstra_shortest_paths(bGraph, s, &(this->m_PreviousVertexList[0]), &distance[0], weightmap, indexmap, 
                                std::less<float>(), closed_plus<float>(), 
                        	(std::numeric_limits<int>::max)(), 0,
                        	default_dijkstra_visitor());
  std::cout << "Done." << std::endl;
  
  
  int endNodeNum = this->m_TargetNodes[0];

  
  for (unsigned long i=0; i<m_Output->GetTotalNumberOfNodes(); i++)
    {
    NodePointerType node = m_Output->GetNodePointer( i );
    if (distance[i] < (std::numeric_limits<int>::max)() )
      {
      node->Weight = distance[i];
      node->Visited = true;
      }  
    else
      {
      node->Weight = 0;
      node->Visited = false;
      }
    }
  
  /*
  BoostVertexDescriptor c = vertex(endNodeNum, bGraph);
	while (indexmap[c] != indexmap[s])
	{
		int tmp_c = indexmap[c];

    std::cout << tmp_c << " (" << distance[tmp_c] << ") ";  

    c  = vertex(indexmap[previousNodeMap[c]], bGraph);
		
	}
  */
  return;
  
  
  

  /* NON-BOOST IMPLEMENTATION
  this->m_CurrentNode = m_Output->GetNodePointer( m_SourceNodes[0] );
  this->m_CurrentNode->IsActive = true;
  this->m_CurrentNode->Weight = 0;
  EdgeListType neighbors = this->GetEdgesToUnvisitedNeighbors( this->m_CurrentNode );
 
  while ( neighbors.size() > 0 )
    {
    
    typename EdgeListType::const_iterator edgeIt = neighbors.begin();
    
    while ( edgeIt != neighbors.end() )
      {
      float currentDistance = this->m_CurrentNode->Weight;
      float edgeDistance = (*edgeIt)->Weight;
      
      if ( (this->m_Output->GetNodePointer( (*edgeIt)->TargetIdentifier )->Weight < 0) ||
           (this->m_Output->GetNodePointer( (*edgeIt)->TargetIdentifier )->Weight > (currentDistance + edgeDistance)) )
        {
        NodeType * neighborNode = this->m_Output->GetNodePointer( (*edgeIt)->TargetIdentifier );
        neighborNode->Weight = currentDistance + edgeDistance;
        if (neighborNode->Paths.size() < 1)
          {
          neighborNode->Paths.push_back( (*edgeIt)->Identifier);
          }
        else 
          {
          neighborNode->Paths[0] = (*edgeIt)->Identifier;
          }
        }     

      ++edgeIt;   
      }
    
    neighbors.clear();
    this->m_CurrentNode->Visited = true;
    NodeIdentifierType lastNode = this->m_CurrentNode->Identifier;
    
    if (this->AdvanceCurrentNode())
      {
      this->m_CurrentNode->HasPath = true;
      neighbors = this->GetEdgesToUnvisitedNeighbors( this->m_CurrentNode );
      }
    else
      {
      std::cout << "Done advancing" << std::endl;
      }
    }
    */
}

/*
template <class TGraph>
bool
DijkstrasPathGraphFilter<TGraph>
::AdvanceCurrentNode()
{
  float minDistance = -1;
  bool foundNode = false;
  unsigned long nextIdx;
  unsigned long unVisited = 0;
  
  for (unsigned long idx=0; idx < m_Output->GetTotalNumberOfNodes(); idx++)
    {
               
    NodePointerType testNode = m_Output->GetNodePointer(idx);
    
    if (!testNode->Visited) ++unVisited;
    
    if ( (!testNode->Visited) && (testNode->Weight > 0) )
      {
      
      EdgeListType neighbors = this->GetEdgesToUnvisitedNeighbors( testNode );
      if (neighbors.size() == 0)
        {
        testNode->Visited = true;
        }
      else 
        {
        if (!foundNode)
          {
          nextIdx = idx;
          minDistance = testNode->Weight;
          foundNode = true;
          //std::cout << "FOUND NODE: " << nextIdx << std::endl;
          }
        else
          {
          if (testNode->Weight < minDistance)
            {
            nextIdx = idx;
            minDistance = testNode->Weight;
            }
          }
        }
      }
    }  
  
  if ((unVisited % 100) == 0) std::cout << "Unvisited nodes = " << unVisited << std::endl;
  
  if (!foundNode) return false;
  
  this->m_CurrentNode = m_Output->GetNodePointer(nextIdx);
  //std::cout << "Returning node = " << this->m_CurrentNode->Identifier << std::endl;
  return true;
}
*/


template <class TGraph>
typename DijkstrasPathGraphFilter<TGraph>::NodeListType
DijkstrasPathGraphFilter<TGraph>
::BackTrace(NodeIdentifierType target)
{
  
  NodeListType nodeList;
  
  NodePointerType current = m_Output->GetNodePointer( target );
  nodeList.push_front(current);
  while (!current->IsSource && current->HasPath)
    {
 
      EdgeType * pathEdge = this->m_Output->GetEdgePointer( current->Paths[0] );

      //std::cout << current->Identifier << " ( " << pathEdge->Identifier << " ) " << std::flush;

      current = this->m_Output->GetNodePointer( pathEdge->SourceIdentifier );
      nodeList.push_front(current);
            
    }
  //std::cout << current->Identifier << std::endl;
  
  if (!current->IsSource)
    {
    std::cout << "Warning (You dun goofed) -> No possible pathway found for given target" << std::endl;
  }
  
  return nodeList;
  
}


template <class TGraph>
typename DijkstrasPathGraphFilter<TGraph>::EdgeListType
DijkstrasPathGraphFilter<TGraph>
::GetEdgesToUnvisitedNeighbors(NodePointerType node)
{
  EdgeListType edges;
  typename EdgeIdentifierContainerType::const_iterator edgeIt = node->OutgoingEdges.begin();
  
  while ( edgeIt != node->OutgoingEdges.end() )
    {
    if ( !m_Output->GetNode( m_Output->GetEdge(*edgeIt).TargetIdentifier ).Visited )
      {
      edges.push_back( m_Output->GetEdgePointer( m_Output->GetEdge(*edgeIt).Identifier ) );
      }
      ++edgeIt;
    }
    
  return edges;
}

template <class TGraph>
void
DijkstrasPathGraphFilter<TGraph>
::Initialize()
{
  
  if (this->GetDebug())
    {
      std::cout << "Initialize()" << std::endl;
      std::cout << "# sources = " << this->m_SourceNodes.size() << std::endl;
      std::cout << "# targets = " << this->m_TargetNodes.size() << std::endl; 
    }
  
  this->AllocateOutputs();
  this->m_Output = this->GetOutput();

  /** Set other node parameters */
  NodeIteratorType It( this->m_Output );
  NodePointerType node;
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    node = It.GetPointer();
    node->IsActive = false;
    node->TimeStamp = 0;
    node->IsSource = false;
    node->IsTarget = false;
    node->Visited = false;
    node->Weight = -1.0;
    node->HasPath = false;
     
    // Order edges by minimum weights
    //this->SortEdgeIdentifiers( & node->IncomingEdges );
    //this->SortEdgeIdentifiers( & node->OutgoingEdges );
      
    }
    
  for (unsigned int i=0; i<this->m_SourceNodes.size(); i++)
    {
    this->m_Output->GetNode( this->m_SourceNodes[i] ).IsSource = true;    
    }
  for (unsigned int i=0; i<this->m_TargetNodes.size(); i++)
    {
    this->m_Output->GetNode( this->m_TargetNodes[i] ).IsTarget = true;  
    }
}

/*
template <class TGraph>
void
DijkstrasPathGraphFilter<TGraph>
::SortEdgeIdentifiers( EdgeIdentifierContainerType * edgeIds )
{
  
  std::vector< EdgeType > edges(edgeIds->size());
  
  for (unsigned int i=0; i<edges.size(); i++)
    {
    edges[i] = m_Output->GetEdge( (*edgeIds)[i] );
    }
  
  std::sort( edges.begin(), edges.end(), GraphTraitsType::EdgeSortPredicate );
  for (unsigned int i=0; i<edges.size(); i++)
    {
    (*edgeIds)[i] = edges[i].Identifier;
    }
  
}
*/

template <class TGraph>
void
DijkstrasPathGraphFilter<TGraph>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

}
} // end namespace itk

#endif
