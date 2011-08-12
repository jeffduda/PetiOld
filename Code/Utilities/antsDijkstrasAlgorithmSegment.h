/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit (ITK)
  Module:    $RCSfile: antsDijkstrasAlgorithmSegment.h,v $
  Language:  C++
  Date:      $Date: 2009/10/20 20:16:43 $
  Version:   $Revision: 1.1 $

=========================================================================*/
#ifndef _antsDijkstrasAlgorithmSegment_h_
#define _antsDijkstrasAlgorithmSegment_h_

#include "itkImageFileReader.h" 
#include "itkImageFileWriter.h"
#include "antsDijkstrasAlgorithm.h"
//#include "itkTensorDeviatoricEuclideanDistanceRegistrationFunction.h"

namespace itk {
namespace ants {

template<class TGraphSearchNode, typename TCostType>
class DijkstrasAlgorithmSegment : 
public DijkstrasAlgorithm<TGraphSearchNode>
{
public: 
  typedef DijkstrasAlgorithmSegment             Self;                   
  typedef DijkstrasAlgorithm<TGraphSearchNode> Superclass;         
  typedef SmartPointer<Self>       Pointer;   
  typedef SmartPointer<const Self> ConstPointer;  
  itkTypeMacro(DijkstrasAlgorithmSegment,DijkstrasAlgorithm);
  itkNewMacro(Self);  

  typedef typename Superclass::NodeLocationType NodeLocationType;
  typedef typename Superclass::NodeListType NodeListType;

  enum{GraphDimension=TGraphSearchNode::GraphDimension}; /* dimension of the graph */

  typedef itk::Image<TCostType,GraphDimension>      GraphCostImageType; 
  typedef itk::Vector<float, 6> TensorType;
  typedef itk::Vector<float, 3> VectorType;
  typedef itk::Image<TensorType, GraphDimension> TensorImageType;
  typedef   typename TensorImageType::Pointer TensorImagePointer;
  typedef itk::Image<VectorType, GraphDimension> VectorImageType;
  typedef   typename VectorImageType::Pointer VectorImagePointer;


  typedef typename TGraphSearchNode::PixelType        PixelType; /*  pixel type for the cost */
  typedef typename TGraphSearchNode::CoordRep         CoordRep;  /* coordinate type */
  typedef itk::NeighborhoodIterator<GraphCostImageType>  CostImageNeighborhoodIteratorType; 
  typedef typename GraphCostImageType::IndexType                      CostImageIndexType;

  // FILE INPUT 
  typedef   itk::ImageFileReader<GraphCostImageType> ReaderType;  
  typename ReaderType::Pointer CostImageReader;

  /* re-implement local cost */
  PixelType LocalCost() ; /* computes the local cost */
  
  inline std::vector< itk::Vector<float,3u> > GetPath() 
  { 
    std::vector< itk::Vector<float,3u> > path;
    path.clear();
    int psz = this->GetPathSize();
    for(int i=0; i< psz; i++)
    {
	    path.push_back(this->GetPathAtIndex(i)->GetLocation());
    }
    return path;
  }
  
  inline void SetSinkLabel( unsigned int  i ) { this->m_SinkLabel=i; }
  inline typename GraphCostImageType::Pointer GetGraphCostImage(){return m_GraphCostImage;}
  
  inline typename GraphCostImageType::Pointer GetTotalCostImage(){return m_TotalCostImage;}
  
   void  SetPriorImage(  typename GraphCostImageType::Pointer pr ){  this->m_PriorImage=pr;}
   void  SetLabelImage(  typename GraphCostImageType::Pointer pr ){  this->m_LabelImage=pr;}
  
  inline void SetGraphCostImage( GraphCostImageType* GCI) 
  {  
	  m_GraphCostImage=GCI;	

	  this->m_TotalCostImage=GraphCostImageType::New();

	  this->m_TotalCostImage->SetLargestPossibleRegion(  GCI->GetLargestPossibleRegion()  ); 
	  this->m_TotalCostImage->SetRequestedRegion(  GCI->GetLargestPossibleRegion() );  
	  this->m_TotalCostImage->SetBufferedRegion(  GCI->GetLargestPossibleRegion() );
	  this->m_TotalCostImage->SetSpacing(GCI->GetSpacing());
	  this->m_TotalCostImage->SetOrigin(GCI->GetOrigin());
	  this->m_TotalCostImage->SetDirection(GCI->GetDirection());
	  this->m_TotalCostImage->Allocate(); 
	  this->m_TotalCostImage->FillBuffer(0); 

  }
  inline void SetTensorImage( TensorImageType* GCI) 
  {  
	  m_TensorImage=GCI;	
  }
  inline void SetVectorImage( VectorImageType* GCI) 
  {  
	  m_VectorImage=GCI;	
  }
 
  inline void SetLocalCostFP(PixelType (*FP)(Pointer))
  {
    m_LocalCostFP=FP;
  }  
  inline void SetWeights( float w1, float w2, float w3 )
  {
    this->m_Weight1=w1; //std::cout << "w1 = " << w1 << std::endl;
    this->m_Weight2=w2; //std::cout << "w2 = " << w2 << std::endl;
    this->m_Weight3=w3; //std::cout << "w3 = " << w3 << std::endl;
  }  
  
  inline void SetExponents( float e1, float e2, float e3 )
  {
    this->m_Exp1=e1; //std::cout << "e1 = " << e1 << std::endl;
    this->m_Exp2=e2; //std::cout << "e2 = " << e2 << std::endl;
    this->m_Exp3=e3; //std::cout << "e3 = " << e3 << std::endl;
  }  
  
  virtual bool TerminationCondition();  /** decides when the algorithm stops */

  //  void writeimage();

  void SetFAThreshold( float f ) { this->m_FAThreshold=f; }

protected:
  DijkstrasAlgorithmSegment()
    { 
      m_PriorImage=NULL;
      m_LabelImage=NULL;
      m_LocalCostFP=NULL;  m_TensorImage=NULL; m_VectorImage=NULL; m_GraphCostImage=NULL; 
      m_Weight1=m_Weight2=m_Weight3=1; 
      m_Exp1=m_Exp2=m_Exp3=1;
      this->m_DTMetric=NULL;
      this->m_FAThreshold=0.25;
      this->m_MaxCost=1.e9;
    }
  ~DijkstrasAlgorithmSegment(){}

private:
  DijkstrasAlgorithmSegment(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented  

  unsigned int m_SinkLabel;
  CostImageIndexType  m_CostIndex; 
  GraphCostImageType* m_GraphCostImage;  /* if we are searching an image, set this */
  PixelType (*m_LocalCostFP)(Pointer);  // local cost function pointer 
  TensorImagePointer  m_TensorImage;
  VectorImagePointer  m_VectorImage;
  typename GraphCostImageType::Pointer  m_TotalCostImage;
  typename GraphCostImageType::Pointer  m_PriorImage;
  typename GraphCostImageType::Pointer  m_LabelImage;
  float m_Weight1;
  float m_Weight2;
  float m_Weight3;
  float m_Exp1;
  float m_Exp2;
  float m_Exp3;
  float m_FAThreshold;

};

}
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "antsDijkstrasAlgorithmSegment.hxx"
#endif

#endif
