/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit (ITK)
  Module:    $RCSfile: antsDijkstrasAlgorithm.cxx,v $
  Language:  C++
  Date:      $Date: 2009/10/20 20:16:43 $
  Version:   $Revision: 1.1 $

=========================================================================*/
#ifndef _antsDijkstrasAlgorithm_hxx_
#define _antsDijkstrasAlgorithm_hxx_

#include "antsDijkstrasAlgorithm.h"

namespace itk {
namespace ants {

template<class TGraphSearchNode > 
DijkstrasAlgorithm<TGraphSearchNode>::DijkstrasAlgorithm()
{
  m_Graph = GraphType::New();
  m_QS = DijkstrasAlgorithmQueue<TGraphSearchNode>::New();
  m_MaxCost=vnl_huge_val(m_MaxCost); // not defined for unsigned char
  this->m_TotalCost=0;
};


template<class TGraphSearchNode > 
void DijkstrasAlgorithm<TGraphSearchNode>::SetGraphSize(GraphSizeType Sz)
{
  for (int i=0; i < GraphDimension; i++)
  {
    m_GraphSize[i]=Sz[i]; 
  }
}


template<class TGraphSearchNode > 
void DijkstrasAlgorithm<TGraphSearchNode>::InitializeGraph()
{
  //bool InBorder=false;
  m_GraphRegion.SetSize( m_GraphSize );
  //m_Graph->SetLargestPossibleRegion( m_GraphRegion ); 
  //m_Graph->SetRequestedRegion( m_GraphRegion );  
  //m_Graph->SetBufferedRegion( m_GraphRegion );
  m_Graph->SetRegions( m_GraphRegion );

  m_Graph->Allocate(); 
  GraphIteratorType GraphIterator( m_Graph, m_GraphRegion );
  GraphIterator.GoToBegin();
  NodeLocationType loc;	
  while(  !GraphIterator.IsAtEnd()  )
  {

    typename GraphSearchNode<PixelType,CoordRep,GraphDimension>::Pointer G=NULL;
    GraphIterator.Set(G);  
    ++GraphIterator; 
    /* 
    m_GraphIndex = GraphIterator.GetIndex();  
    //std::cout << " allocating "  << m_GraphIndex << std::endl;  
	///GraphSearchNode<PixelType,CoordRep,GraphDimension>::Pointer G=
    G=TGraphSearchNode::New();
    G->SetUnVisited();
    G->SetTotalCost(m_MaxCost);
    for (int i=0; i<GraphDimension; i++) loc[i]=m_GraphIndex[i];
    G->SetLocation(loc);
    G->SetPredecessor(NULL);
    m_Graph->SetPixel(m_GraphIndex,G);*/
    m_Graph->SetPixel( GraphIterator.GetIndex(),NULL);  // USE IF POINTER IMAGE defines visited

  }
  m_SearchFinished=false;
}

template<class TGraphSearchNode > 
void DijkstrasAlgorithm<TGraphSearchNode>::InitializeQueue()
{
  int n = m_QS->m_SourceNodes.size();
  GraphIteratorType GraphIterator( m_Graph, m_GraphRegion );
  m_GraphSize=m_Graph->GetLargestPossibleRegion().GetSize();
  GraphIterator.GoToBegin();
  m_GraphIndex = GraphIterator.GetIndex();
  NodeLocationType loc;
  // make sure the graph contains the right pointers 
  for (int i=0; i<n; i++)
  {
    typename GraphSearchNode<PixelType,CoordRep,GraphDimension>::Pointer G = m_QS->m_SourceNodes[i];
	  G->SetPredecessor(G);
  	m_QS->m_Q.push(G);
  	loc=G->GetLocation();
  	for (int d=0;d<GraphDimension;d++) m_GraphIndex[d]=(long int)(loc[d]+0.5);
	m_Graph->SetPixel(m_GraphIndex,G);
  }    
  for (unsigned int i=0; i<m_QS->m_SinkNodes.size(); i++)
  {
    typename GraphSearchNode<PixelType,CoordRep,GraphDimension>::Pointer G = m_QS->m_SinkNodes[i];
	  G->SetPredecessor(NULL);
	  loc=G->GetLocation();
	  for (int d=0;d<GraphDimension;d++) m_GraphIndex[d]=(long)loc[d];
	    m_Graph->SetPixel(m_GraphIndex,G);
  }   
  m_SearchFinished=false;
  if (m_EdgeTemplate.empty()) InitializeEdgeTemplate();
}

template<class TGraphSearchNode > 
void DijkstrasAlgorithm<TGraphSearchNode>::InitializeEdgeTemplate
  (vector < unsigned int > UserEdgeTemplate,unsigned int R)
{ 
  for (int i=0; i<GraphDimension; i++)
  {
	m_Radius[i] = R;
  }

  m_EdgeTemplate=UserEdgeTemplate;
}

template<class TGraphSearchNode > 
void DijkstrasAlgorithm<TGraphSearchNode>::InitializeEdgeTemplate()
{ 
  int MaxIndex=1;
  for (int i=0; i<GraphDimension; i++)
  {
	m_Radius[i]=1;
  }
  for (int i=0; i<GraphDimension; i++)
  {
	MaxIndex=MaxIndex*(2*m_Radius[i]+1);
  }
  MaxIndex=MaxIndex-1;
//  int Center = MaxIndex/2;
  for (int i=0; i<=MaxIndex; i++) 
  { 
  //  if (i != Center) 
	  m_EdgeTemplate.insert(m_EdgeTemplate.end(),i);
  }

}

/** 
*  Compute the local cost using Manhattan distance.
*/
template<class TGraphSearchNode > 
typename DijkstrasAlgorithm<TGraphSearchNode>::
PixelType DijkstrasAlgorithm<TGraphSearchNode>::LocalCost() 
{
  return 1.0; // manhattan distance
};

template<class TGraphSearchNode >
bool DijkstrasAlgorithm<TGraphSearchNode>::TerminationCondition()
{	
  if (!m_QS->m_SinkNodes.empty())
  {
    if (m_NeighborNode == m_QS->m_SinkNodes[0] && !m_SearchFinished  ) 
    {
      m_SearchFinished=true;
	    m_NeighborNode->SetPredecessor(m_CurrentNode);
  	}
  }
  else m_SearchFinished=true;
  return m_SearchFinished;
}


template<class TGraphSearchNode >
void DijkstrasAlgorithm<TGraphSearchNode>::SearchEdgeSet() 
{
  //std::cout << " SES 0 " << std::endl;
  int i=0;
  //int j=0;
  GraphNeighborhoodIteratorType GHood(m_Radius, m_Graph,m_Graph->GetRequestedRegion());
  GraphNeighborhoodIndexType	GNI;
  

  
  //std::cout << " SES 1 " << std::endl;
  for (i=0; i < GraphDimension; i++)
    {   
      //std::cout << " SES 2 " << std::endl;
      GNI[i]=(long)(m_CurrentNode->GetLocation()[i]+0.5);
    } 
  
  //std::cout << " SES 3 " << std::endl;
  GHood.SetLocation(GNI);

  //bool inimage=true;
  for (unsigned int dd=0; dd<GraphDimension; dd++)
    {
      if ( GNI[dd] < 2 || GNI[dd] > (long)(m_GraphSize[dd]-2)) return;
    }


  for (i = 0; i < (int) m_EdgeTemplate.size(); i++)
    {    
    //std::cout << " SES 4 " << std::endl;
    //std::cout << " ET " << m_EdgeTemplate[i]  <<  " RAD " << m_Radius << " ind " << GHood.GetIndex(m_EdgeTemplate[i]) << std::endl;
    if (!GHood.GetPixel(m_EdgeTemplate[i])) //std::cout << " OK " << std::endl;
      ///else
      {
	//	std::cout << " NOT OK  " <<std::endl;
	GraphNeighborhoodIndexType ind = GHood.GetIndex(m_EdgeTemplate[i]);    
	typename TGraphSearchNode::Pointer G=TGraphSearchNode::New();
	G->SetUnVisited();
	G->SetTotalCost(m_MaxCost);
	NodeLocationType loc;
	for (int i=0; i<GraphDimension; i++) loc[i]=ind[i];
	G->SetLocation(loc);
	G->SetPredecessor(m_CurrentNode);
	m_Graph->SetPixel(ind,G);  
      }
    m_NeighborNode=GHood.GetPixel(m_EdgeTemplate[i]);
    //std::cout << "DA  i " << i << " position " << m_NeighborNode->GetLocation() << endl;
    this->TerminationCondition();
    if (!m_SearchFinished && m_CurrentNode != m_NeighborNode &&
	!m_NeighborNode->GetDelivered())
      {
        m_NewCost = m_CurrentCost + LocalCost();
	CheckNodeStatus();
      }
  }
}


template<class TGraphSearchNode > 
void DijkstrasAlgorithm<TGraphSearchNode>::CheckNodeStatus() 
// checks a graph neighbor's status 
{
  if (!m_NeighborNode->GetVisited()  )
    {
      // set the cost and put into the queue
      m_NeighborNode->SetTotalCost(m_NewCost);
      m_NeighborNode->SetPredecessor(m_CurrentNode);
      m_NeighborNode->SetVisited();
      m_QS->m_Q.push(m_NeighborNode);	
    }
  else if (m_NewCost < m_NeighborNode->GetTotalCost() )
    {
      m_NeighborNode->SetTotalCost(m_NewCost);
      m_NeighborNode->SetPredecessor(m_CurrentNode);
      m_QS->m_Q.push(m_NeighborNode);	
    }
  
}

template<class TGraphSearchNode > 
void DijkstrasAlgorithm<TGraphSearchNode>::FindPath()
{               
  if (m_QS->m_SourceNodes.empty())
  {
    std::cout << "ERROR !! DID NOT SET SOURCE!!\n";
    return;
  }
  
  //std::cout << "DA start find path " << " Q size " << m_QS->m_Q.size() << " \n";
  
  while ( !m_SearchFinished && !m_QS->m_Q.empty()  )
  {
    m_CurrentNode=m_QS->m_Q.top();
    m_CurrentCost=m_CurrentNode->GetTotalCost();
    m_QS->m_Q.pop();
    
    if (!m_CurrentNode->GetDelivered())
    {
	    m_QS->IncrementTimer();
	    //std::cout << " searching " << m_CurrentNode->GetLocation()   << " \n";
	    this->SearchEdgeSet();	
	    this->m_TotalCost += m_CurrentNode->GetTotalCost();
	    //std::cout << " " << this->m_TotalCost << std::flush;
	    //if ( (m_CurrentNode->GetTimer() % 1.e5 ) == 0) 
	    // std::cout << " searched  " << m_CurrentNode->GetTimer()   << " \n";
	  }
    m_CurrentNode->SetDelivered();
      
  }  // end of while
  
  m_NumberSearched = (unsigned long) m_QS->GetTimer();
  std::cout << "Done with find path " << " Q size " << m_QS->m_Q.size() << " num searched " << m_NumberSearched << " \n";
  
  return;
  
}

}
} // end namespace itk

#endif

