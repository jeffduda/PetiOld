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

template <class TGraph, class TImage>
DijkstrasPathGraphFilter<TGraph, TImage>
::DijkstrasPathGraphFilter()
{
  this->m_SourceNodes.clear();
  this->m_TargetNodes.clear();
  this->m_Labels = NULL;
}

/** Generate the data */
template <class TGraph, class TImage>
void
DijkstrasPathGraphFilter<TGraph, TImage>
::GenerateData()
{

	// Get region info
	typename LabelStatsType::Pointer labelStats = LabelStatsType::New();
	labelStats->SetInput( this->m_Labels );
	labelStats->Update();
	
	typename LabeledImageType::ValueType minLabel = 1;
	typename LabeledImageType::ValueType maxLabel = 2;
	
	if (labelStats->GetMinimum() > minLabel)
		{
		minLabel = labelStats->GetMinimum();
		}
	if (labelStats->GetMaximum() > maxLabel)
		{
		maxLabel = labelStats->GetMaximum();
		}

	LabelType nLabels = maxLabel - minLabel + 1;

	this->m_PathImage = PathImageType::New();
	typename PathImageType::RegionType pathImageRegion;
	typename PathImageType::RegionType::SizeType pathImageSize;
	pathImageSize[0] = nLabels;
	pathImageSize[1] = nLabels;
	pathImageRegion.SetSize( pathImageSize );
	this->m_PathImage->SetRegions( pathImageRegion );
	this->m_PathImage->Allocate();

	this->m_Tree.clear();
	
	for (LabelType lab1=minLabel; lab1<=maxLabel; lab1++)
		for (LabelType lab2=minLabel; lab2<=maxLabel; lab2++)
			{
				if (lab1 != lab2)
					{
					this->SetSourceLabel( lab1 );
					this->SetTargetLabel( lab2 );
					this->InitializeNodeLists();

					for (unsigned int srcN=0; srcN<this->m_SourceNodes.size(); srcN++)
						{
						this->GenerateDistancesFromSource( this->m_SourceNodes[srcN] );
						this->GrowTree();
						}
			
					typename PathImageType::IndexType pathImageIndex;
					pathImageIndex[0] = lab1-minLabel;
					pathImageIndex[1] = lab2-minLabel;
					//this->m_PathImage->SetPixel(pathImageIndex, this->m_Tree);
					//this->m_Tree.clear();
					}				
			}
	
}

template <class TGraph, class TImage>
void
DijkstrasPathGraphFilter<TGraph, TImage>
::GenerateDistancesFromSource(unsigned int sourceNodeId)
{

	// Reset node values
	this->Initialize();
	this->m_Output->GetNodePointer(sourceNodeId)->IsSource = true;

	BoostGraphType bGraph( m_Output->GetTotalNumberOfNodes() );
	std::vector<BoostEdgeType> edgeVec( m_Output->GetTotalNumberOfEdges() );

	property_map<BoostGraphType, edge_weight_t>::type weightmap = get(edge_weight, bGraph);

	for (unsigned long i=0; i<m_Output->GetTotalNumberOfEdges(); i++)
	  {
	  
	  // Check 
	  EdgePointerType graphEdge = m_Output->GetEdgePointer( i );
	  // FIXME - throw an exception here
		if (graphEdge->Weight < 0)
			{
			std::cerr << "Negative edge weight detected.. exiting" << std::endl;
			return;
	    }
	    	  
		BoostEdgeDescriptor e; 
		bool inserted;
		int curNodeNum = graphEdge->SourceIdentifier;
		int toNodeNum = graphEdge->TargetIdentifier;
		tie(e, inserted) = add_edge(curNodeNum, toNodeNum, bGraph);
	  weightmap[e] = graphEdge->Weight;
	  }

	//int startNodeNum = this->m_SourceNodes[sourceNodeId];    
	this->m_PreviousVertexList.clear();
	this->m_PreviousVertexList.reserve( num_vertices(bGraph) );
	
	std::vector<int> distance(num_vertices(bGraph));
	BoostVertexDescriptor s = vertex(sourceNodeId, bGraph);
	typename property_map<BoostGraphType, vertex_index_t>::type indexmap = get(vertex_index, bGraph);

	dijkstra_shortest_paths_no_color_map(bGraph, s, &(this->m_PreviousVertexList[0]), &distance[0], weightmap, indexmap, 
	                              std::less<float>(), closed_plus<float>(), 
	                      	(std::numeric_limits<int>::max)(), 0,
	                      	default_dijkstra_visitor());
	 
	for (unsigned long i=0; i<m_Output->GetTotalNumberOfNodes(); i++)
	  {
	  NodePointerType node = m_Output->GetNodePointer( i );
	  if (distance[i] < (std::numeric_limits<int>::max)() )
	    {
	    node->Weight = distance[i];
	    node->Visited = true;
	    node->PreviousNodes.clear();     
	    node->PreviousNodes.push_back( indexmap[this->m_PreviousVertexList[vertex(i,bGraph)]] );
	    }  
	  else
	    {
	    node->Weight = 0;
	    node->Visited = false;
	    }
	  }
}


template <class TGraph, class TImage>
typename DijkstrasPathGraphFilter<TGraph, TImage>::NodeListType
DijkstrasPathGraphFilter<TGraph, TImage>
::BackTrace(NodeIdentifierType target)
{

  NodeListType nodeList;  
  NodePointerType current = m_Output->GetNodePointer( target );
  nodeList.push_front(current);
  while (!current->IsSource && (current->PreviousNodes.size() > 0) )
    {
    current = this->m_Output->GetNodePointer( current->PreviousNodes[0] );
    nodeList.push_front(current);        
    }
  
  // No path found
  if (!current->IsSource)
    {
    nodeList.clear();
  	}
  
  return nodeList;
  
}

template <class TGraph, class TImage>
void
DijkstrasPathGraphFilter<TGraph, TImage>
::GrowTree( )
{  
  for (unsigned int i=0; i<this->m_TargetNodes.size(); i++)
    {
    NodeListType path = this->BackTrace( this->m_TargetNodes[i] );
    if (path.size() > 0)
			{
			this->m_Tree.push_back(path);  
			}
    }
}


template <class TGraph, class TImage>
void
DijkstrasPathGraphFilter<TGraph, TImage>
::InitializeNodeLists( )
{
	this->m_SourceNodes.clear();
	this->m_TargetNodes.clear();

  /** Set other node parameters */
  NodeIteratorType It( this->GetInput() );
  NodePointerType node;
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    node = It.GetPointer();
    
    if (this->m_Labels->GetPixel(node->ImageIndex) == this->m_SourceLabel)
    	{
    	//node->IsSource = true;
    	this->m_SourceNodes.push_back( node->Identifier );
    	}
    if (this->m_Labels->GetPixel(node->ImageIndex) == this->m_TargetLabel)
    	{
    	//node->IsTarget = true;
    	this->m_TargetNodes.push_back( node->Identifier );
    	}    
    }

}

template <class TGraph, class TImage>
void
DijkstrasPathGraphFilter<TGraph, TImage>
::Initialize( )
{
 
  this->AllocateOutputs();
  this->m_Output = this->GetOutput();

  /** Set other node parameters */
  NodeIteratorType It( this->m_Output );
  NodePointerType node;
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    node = It.GetPointer();
    node->TimeStamp = 0;
    node->IsSource = false;
    node->IsTarget = false;
    node->Visited = false;
    node->Weight = -1.0;
    node->PreviousNodes.clear();
    }

}

template <class TGraph, class TImage>
void
DijkstrasPathGraphFilter<TGraph, TImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

}
} // end namespace itk

#endif
