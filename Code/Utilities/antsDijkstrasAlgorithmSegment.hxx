/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit (ITK)
  Module:    $RCSfile: antsDijkstrasAlgorithmSegment.hxx,v $
  Language:  C++
  Date:      $Date: 2009/10/20 20:16:43 $
  Version:   $Revision: 1.1 $

=========================================================================*/
#ifndef _antsDijkstrasAlgorithmSegment_hxx_
#define _antsDijkstrasAlgorithmSegment_hxx_

#include "antsDijkstrasAlgorithmSegment.h"

namespace itk {
namespace ants {

/*
 * compute the local cost
 */
template<class TGraphSearchNode, typename TCostImageType > 
typename DijkstrasAlgorithmSegment<TGraphSearchNode,TCostImageType>::PixelType 
DijkstrasAlgorithmSegment<TGraphSearchNode,TCostImageType>::LocalCost() 
{
  
  if (!m_LocalCostFP) 
  {
    NodeLocationType loc1=this->m_CurrentNode->GetLocation();
    NodeLocationType loc2=this->m_NeighborNode->GetLocation();
    NodeLocationType stepsize=loc2-loc1;
    float mag=0;
    for (int i=0; i<GraphDimension; i++) 
    {
	    this->m_CostIndex[i]=loc2[i];
	    mag+=stepsize[i]*stepsize[i];
    }
    return sqrt(mag);
  }
  
  else return this->m_LocalCostFP(this);
  
};


template<class TGraphSearchNode, typename TCostImageType > 
bool
DijkstrasAlgorithmSegment<TGraphSearchNode,TCostImageType>::TerminationCondition() 
{
  
  typename TensorImageType::IndexType ind;
  NodeLocationType loc=this->m_CurrentNode->GetLocation();
  for (int i=0; i<GraphDimension; i++)  ind[i]=loc[i];

  this->m_SearchFinished=false;
  
  //  if the current node is a sink node, then you are done 
//  if (this->m_LabelImage)
  //  if ( (unsigned int) this->m_LabelImage->GetPixel(ind) == this->m_SinkLabel)
  //   this->m_SearchFinished=true;

/*
  for (int i=0; i<this->m_QS->m_SinkNodes.size(); i++)
    {
    typename GraphSearchNode<PixelType,CoordRep,GraphDimension>::Pointer G = this->m_QS->m_SinkNodes[i];
    if (this->m_CurrentNode == G) this->m_SearchFinished=true;
    }
  */
  
  if ( this->m_CurrentNode->GetTotalCost() >= this->m_MaxCost) this->m_SearchFinished=true;
  
  return this->m_SearchFinished;
  
}

}
} // end namespace itk

#endif
