/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit (ITK)
  Module:    $RCSfile: itkDijkstrasAlgorithmSegment.cxx,v $
  Language:  C++
  Date:      $Date: 2009/10/20 20:16:43 $
  Version:   $Revision: 1.1 $

=========================================================================*/
#ifndef _itkDijkstrasAlgorithmSegment_cxx_
#define _itkDijkstrasAlgorithmSegment_cxx_

#include "itkDijkstrasAlgorithmSegment.h"
#include "TensorFunctions.h"

namespace itk {

/*
 * compute the local cost
 */
template<class TGraphSearchNode, typename TCostImageType > 
typename DijkstrasAlgorithmSegment<TGraphSearchNode,TCostImageType>::PixelType 
DijkstrasAlgorithmSegment<TGraphSearchNode,TCostImageType>::LocalCost() 
{
  
  /* New testing metric */
  if (m_TensorImage)
  {
    typename TensorImageType::SpacingType spacing=m_TensorImage->GetSpacing();
    typename TensorImageType::IndexType ind1,ind2,ind3;
    NodeLocationType loc1=this->m_CurrentNode->GetPredecessor()->GetLocation();
    NodeLocationType loc2=this->m_CurrentNode->GetLocation();
    NodeLocationType loc3=this->m_NeighborNode->GetLocation();
    
    typedef itk::Vector<float,3> VectorType;
    itk::Vector<float,3> directionOfPath;
    itk::Vector<float,3> NewDeigv;
    itk::Vector<float,3> CurDeigv;
    
    VectorType currentEigenvector;
    VectorType prevDirectionOfPath;
    
    float cost = 0.0;
    float stepmag=0;
    float ndderiv[GraphDimension];
    
    for (int i=0; i<GraphDimension; i++) 
  	{
  	  ind1[i]=loc1[i];
  	  ind2[i]=loc2[i];
  	  ind3[i]=loc3[i];
  	  this->m_CostIndex[i]=loc2[i];
  	  directionOfPath[i]=((float)loc3[i]-(float)loc2[i])*spacing[i];  // direction (current --> next)
  	  prevDirectionOfPath[i]=((float)loc2[i]-(float)loc1[i])*spacing[i]; // direction (prev --> current)
  	  stepmag += (float)directionOfPath[i]*(float)directionOfPath[i];
  	}
  	
    stepmag = sqrt(stepmag);
    if (stepmag == 0)
      stepmag=1;
    directionOfPath /= stepmag;  //  make unit vector (direction from current to next)
    float lengthCost = stepmag;

  


    float fa = this->m_GraphCostImage->GetPixel( ind3 );
    if (fa < 0) fa = 0.0;
    if (fa > 1) fa = 1.0;
    
    float faScale = 3.0;
    float faCost = cost + faScale*(1.0-fa);
    
    float dirScale = 3.0;

    currentEigenvector = this->m_VectorImage->GetPixel(ind2);
    
    float pathToCurrentDirDP = currentEigenvector*directionOfPath;
    float pathCurvatureDP = directionOfPath * prevDirectionOfPath;

    pathToCurrentDirDP = vcl_fabs(pathToCurrentDirDP);
    float dirCost = cost + dirScale*(1.0-pathToCurrentDirDP);
    
    float curvatureScale = 3.0;
    float curvatureCost = curvatureScale*(-pathCurvatureDP);
    //if (pathCurvatureDP < 0.5) return this->m_MaxCost;

    float fddScale = 3.0;
    float fdd = DiffusionCoefficient<TensorType,VectorType>(this->m_TensorImage->GetPixel(ind2),directionOfPath,true);
    float fddCost = fddScale*(1.0-fdd);
    
    cost = lengthCost + faCost + dirCost;
    
    if (this->m_PriorImage)
    {
      float priorValue = this->m_PriorImage->GetPixel(ind3);
      if (priorValue == 0)
        cost = this->m_MaxCost;    
        
    }

    // FIXME -- possibly to strict here
    if ( this->m_GraphCostImage->GetPixel( ind3 ) < this->m_FAThreshold ) 
  	{
  	  cost = this->m_MaxCost;  //  cost very HIGH if below FA threshold
  	}
    

    // FIXME - what is this?
    this->m_TotalCostImage->SetPixel(ind2,this->m_CurrentNode->GetTotalCost());

    return cost;
   
    
    // VectorImage is principle eigenvector field
    VectorType e1,e2,e3;
    e1=this->m_VectorImage->GetPixel(ind1);
    e2=this->m_VectorImage->GetPixel(ind2);
    e3=this->m_VectorImage->GetPixel(ind3);
    
    
    float kappa2=0;
    float newdeigmag=0;
    float curdeigmag=0;
    itk::Vector<float,3> deigv2;
    
    
    for (int i=0; i<GraphDimension; i++) 
  	{
  	  NewDeigv[i]=e3[i];  //*spacing[i];
  	  CurDeigv[i]=e2[i];  //*spacing[i];
  	  newdeigmag+=NewDeigv[i]*NewDeigv[i];
  	  curdeigmag+=CurDeigv[i]*CurDeigv[i];
  	  ndderiv[i]=2.0*(float)e2[i]-(float)e1[i]-(float)e3[i];
  	  kappa2+=ndderiv[i]*ndderiv[i];
  	}
    newdeigmag=sqrt(newdeigmag);
    curdeigmag=sqrt(curdeigmag);
    if (newdeigmag == 0) newdeigmag=1;
    if (curdeigmag == 0) curdeigmag=1;
    NewDeigv/=newdeigmag; // make unit length  next WS principal eigen vector
    CurDeigv/=curdeigmag;// make unit length  current WS principal eigen vector
    kappa2=sqrt(kappa2);
    
    
    // get dot products of eigenvectors and path dir
    // account for antiparallel vectors
    float newunitvdp=0;
    float newunitvdp2=0;
    float curunitvdp=0;
    float curunitvdp2=0;
    for (int i=0; i<GraphDimension; i++) 
  	{
  	  float uv1=directionOfPath[i];
  	  float uv2=NewDeigv[i];
  	  newunitvdp+=uv1*uv2;
  	  newunitvdp2+=uv1*(-1.0*uv2);
  	  uv1=CurDeigv[i];
  	  curunitvdp+=uv1*uv2;
  	  curunitvdp2+=uv1*(-1.0*uv2);
  	}

    if (newunitvdp2 > newunitvdp) newunitvdp=newunitvdp2;
    if (curunitvdp2 > curunitvdp) curunitvdp=curunitvdp2;

    float w1=this->m_Weight1;
    float w2=this->m_Weight2;
    float w3=this->m_Weight3;
    float ex1=this->m_Exp1;
    float ex2=this->m_Exp2;
    float ex3=this->m_Exp3;

    float faval = (1.0-this->m_GraphCostImage->GetPixel(ind3));
    //if (faval < 0.01) faval=0.0;
    if (faval < 0.1) faval = 0.1;

    //float fdd = GetFiberDirectedDiffusion<TensorType,VectorType>( m_TensorImage->GetPixel(ind2), directionOfPath, true );

    float multfac= pow((double)fabs(faval),(double)ex1);
    float term2 = pow((double)(1.0-curunitvdp),(double)ex2);
    float term3 = pow((double)(1.0-newunitvdp),(double)ex3);
    cost  = w1 * multfac * ( w2*term2   + w3* term3 );
    
    //cost = w1*multfac + w2*term2 + w3*term3;

    // FIXME - this makes much of the above work pointless
    //if (w1 > 0) cost = TensorFunctions<TensorType>GetMetricTensorCost( directionOfPath, this->m_TensorImage->GetPixel(ind3), w1);
    
    // FIXME - not sure why this is happening
    if (this->m_CurrentNode) this->BackTrack(this->m_CurrentNode);
    float psz=(float)this->GetPathSize()+1;
    cost=cost/fabs(psz)*1.e2;

    if (this->m_PriorImage)
  	{
    	float priorval=this->m_PriorImage->GetPixel(ind3);
    	bool opt1=false;
    	if (opt1)
    	{	
    	  if (priorval < this->m_FAThreshold) return this->m_MaxCost;
    	  priorval=pow(priorval,2);
    	  cost/=priorval;
    	}
    	else
    	{
    	  float diff=(priorval-0.75);
    	  if ( diff < 0 ) 
    	  {
    	    diff*=diff;
    	    priorval= exp(-1.0*diff/0.1);
    	  }
    	if (priorval <  0.01 ) return this->m_MaxCost;
    	cost/=priorval;
    	}
  	}
      
    return cost;
    
  }  //   if (m_TensorImage)



  /*
  if (m_TensorImage)
  {
    typename TensorImageType::SpacingType spacing=m_TensorImage->GetSpacing();
    typename TensorImageType::IndexType ind1,ind2,ind3;
    NodeLocationType loc1=this->m_CurrentNode->GetPredecessor()->GetLocation();
    NodeLocationType loc2=this->m_CurrentNode->GetLocation();
    NodeLocationType loc3=this->m_NeighborNode->GetLocation();
    
    typedef itk::Vector<float,3> VectorType;
    itk::Vector<float,3> directionOfPath;
    itk::Vector<float,3> NewDeigv;
    itk::Vector<float,3> CurDeigv;
    
    float stepmag=0;
    float ndderiv[GraphDimension];
    
    for (int i=0; i<GraphDimension; i++) 
  	{
  	  ind1[i]=loc1[i];
  	  ind2[i]=loc2[i];
  	  ind3[i]=loc3[i];
  	  this->m_CostIndex[i]=loc2[i];
  	  directionOfPath[i]=((float)loc3[i]-(float)loc2[i])*spacing[i];  // direction (current --> next)
  	  stepmag+=(float)directionOfPath[i]*(float)directionOfPath[i];
  	}
  	
    stepmag=sqrt(stepmag);
    if (stepmag == 0)stepmag=1;
    directionOfPath/=stepmag;  //  make unit vector (direction from current to next)
    
    // FIXME -- possibly to strict here
    //if ( this->m_GraphCostImage->GetPixel(  ind3 ) <   this->m_FAThreshold ) 
  	//{
  	//  return this->m_MaxCost;  //  cost very HIGH if below FA threshold
  	//}


    // FIXME - what is this?
    this->m_TotalCostImage->SetPixel(ind2,this->m_CurrentNode->GetTotalCost());


    
    // VectorImage is principle eigenvector field
    VectorType e1,e2,e3;
    e1=this->m_VectorImage->GetPixel(ind1);
    e2=this->m_VectorImage->GetPixel(ind2);
    e3=this->m_VectorImage->GetPixel(ind3);
    
    
    float kappa2=0;
    float newdeigmag=0;
    float curdeigmag=0;
    itk::Vector<float,3> deigv2;
    
    
    for (int i=0; i<GraphDimension; i++) 
  	{
  	  NewDeigv[i]=e3[i];  //*spacing[i];
  	  CurDeigv[i]=e2[i];  //*spacing[i];
  	  newdeigmag+=NewDeigv[i]*NewDeigv[i];
  	  curdeigmag+=CurDeigv[i]*CurDeigv[i];
  	  ndderiv[i]=2.0*(float)e2[i]-(float)e1[i]-(float)e3[i];
  	  kappa2+=ndderiv[i]*ndderiv[i];
  	}
    newdeigmag=sqrt(newdeigmag);
    curdeigmag=sqrt(curdeigmag);
    if (newdeigmag == 0) newdeigmag=1;
    if (curdeigmag == 0) curdeigmag=1;
    NewDeigv/=newdeigmag; // make unit length  next WS principal eigen vector
    CurDeigv/=curdeigmag;// make unit length  current WS principal eigen vector
    kappa2=sqrt(kappa2);
    
    
    // get dot products of eigenvectors and path dir
    // account for antiparallel vectors
    float newunitvdp=0;
    float newunitvdp2=0;
    float curunitvdp=0;
    float curunitvdp2=0;
    for (int i=0; i<GraphDimension; i++) 
  	{
  	  float uv1=directionOfPath[i];
  	  float uv2=NewDeigv[i];
  	  newunitvdp+=uv1*uv2;
  	  newunitvdp2+=uv1*(-1.0*uv2);
  	  uv1=CurDeigv[i];
  	  curunitvdp+=uv1*uv2;
  	  curunitvdp2+=uv1*(-1.0*uv2);
  	}

    if (newunitvdp2 > newunitvdp) newunitvdp=newunitvdp2;
    if (curunitvdp2 > curunitvdp) curunitvdp=curunitvdp2;

    float w1=this->m_Weight1;
    float w2=this->m_Weight2;
    float w3=this->m_Weight3;
    float ex1=this->m_Exp1;
    float ex2=this->m_Exp2;
    float ex3=this->m_Exp3;

    float faval = (1.0-this->m_GraphCostImage->GetPixel(ind3));
    //if (faval < 0.01) faval=0.0;
    if (faval < 0.1) faval = 0.1;

    //float fdd = GetFiberDirectedDiffusion<TensorType,VectorType>( m_TensorImage->GetPixel(ind2), directionOfPath, true );

    float multfac= pow((double)fabs(faval),(double)ex1);
    float term2 = pow((double)(1.0-curunitvdp),(double)ex2);
    float term3 = pow((double)(1.0-newunitvdp),(double)ex3);
    float cost  = w1 * multfac * ( w2*term2   + w3* term3 );
    
    //cost = w1*multfac + w2*term2 + w3*term3;

    // FIXME - this makes much of the above work pointless
    //if (w1 > 0) cost = TensorFunctions<TensorType>GetMetricTensorCost( directionOfPath, this->m_TensorImage->GetPixel(ind3), w1);
    
    // FIXME - not sure why this is happening
    if (this->m_CurrentNode) this->BackTrack(this->m_CurrentNode);
    float psz=(float)this->GetPathSize()+1;
    cost=cost/fabs(psz)*1.e2;

    if (this->m_PriorImage)
  	{
    	float priorval=this->m_PriorImage->GetPixel(ind3);
    	bool opt1=false;
    	if (opt1)
    	{	
    	  if (priorval < this->m_FAThreshold) return this->m_MaxCost;
    	  priorval=pow(priorval,2);
    	  cost/=priorval;
    	}
    	else
    	{
    	  float diff=(priorval-0.75);
    	  if ( diff < 0 ) 
    	  {
    	    diff*=diff;
    	    priorval= exp(-1.0*diff/0.1);
    	  }
    	if (priorval <  0.01 ) return this->m_MaxCost;
    	cost/=priorval;
    	}
  	}
      
    return cost;
    
  }  //   if (m_TensorImage)
  */

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

} // end namespace itk

#endif
