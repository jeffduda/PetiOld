/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit (ITK)
  Module:    $RCSfile: itkGeodesicPaths.h,v $
  Language:  C++
  Date:      $Date: 2009/10/20 20:16:43 $
  Version:   $Revision: 1.1 $

=========================================================================*/
#ifndef _itkGeodesicPaths_h_
#define _itkGeodesicPaths_h_

#include "itkImageFileReader.h" 
#include "itkImageFileWriter.h"
#include <itkRawImageIO.h>
#include "itkDijkstrasAlgorithm.h"
namespace itk {

template<class TGraphSearchNode, typename TCostImageType=float>
class GeodesicPaths : 
public DijkstrasAlgorithm<TGraphSearchNode>
{
public: 
  typedef GeodesicPaths             Self;                   
  typedef DijkstrasAlgorithm<TGraphSearchNode> Superclass;         
  typedef SmartPointer<Self>       Pointer;   
  typedef SmartPointer<const Self> ConstPointer;  
  itkTypeMacro(GeodesicPaths,DijkstrasAlgorithm);
  itkNewMacro(Self);  
  
  typedef typename Superclass::GraphNeighborhoodIteratorType 
     GraphNeighborhoodIteratorType;
  
  typedef typename Superclass::GraphIteratorType 
     GraphIteratorType;

  typedef typename Superclass::GraphNeighborhoodIndexType 
     GraphNeighborhoodIndexType;

  

  typedef typename Superclass::NodeLocationType NodeLocationType;
  enum{GraphDimension=TGraphSearchNode::GraphDimension}; /* dimension of the graph */

  typedef itk::Image<TCostImageType,GraphDimension>      GraphCostImageType; 
  typedef typename GraphCostImageType::PixelType      GraphCostImagePixelType; 
  typedef typename TGraphSearchNode::PixelType        PixelType; /*  pixel type for the cost */
  typedef typename TGraphSearchNode::CoordRep         CoordRep;  /* coordinate type */
  typedef itk::NeighborhoodIterator<GraphCostImageType>  CostImageNeighborhoodIteratorType; 
  typedef typename GraphCostImageType::IndexType                      CostImageIndexType;


  typedef itk::Vector<float,3>         VectorType;
  typedef itk::Image<VectorType,3>     FieldType;
  typedef typename FieldType::Pointer FieldTypePointer;


  
  // FILE INPUT 
  typedef   itk::ImageFileReader<GraphCostImageType> ReaderType;  
  typename ReaderType::Pointer CostImageReader;

  /* re-implement local cost */
  PixelType LocalCost() ; /* computes the local cost */
  
  inline typename GraphCostImageType::Pointer GetGraphCostImage(){return this->m_GraphCostImage;}
  
  inline void SetGraphCostImage( typename GraphCostImageType::Pointer GCI) 
  {  
	  this->m_GraphCostImage=GCI;	
    this->SetGraphSize(GCI->GetLargestPossibleRegion().GetSize()); 
    PixelType tt=0;
    this->SetMaxCost(vnl_huge_val(tt));	 
  }
 
  inline void SetLocalCostFP(PixelType (*FP)(Pointer))
  {
    this->m_LocalCostFP=FP; 
  }  
  
  void SearchEdgeSet();

  void writeimage(const char* fn);

  void SetSurfaceLabel(GraphCostImagePixelType p) { this->m_SurfaceLabel=p;}
  void SetSourceLabel(GraphCostImagePixelType p) { this->m_SourceLabel=p;}

  void SetRadius(unsigned int r) 
    { for (int i=0; i<GraphDimension; i++) this->m_Radius[i]=r; }

  void MeanGeodesic(unsigned long = 0, float frac=0.1);

  void ConnectedSurface(float frac);
  void LabelConnectedSurface(float frac=0.95);
  
  void GenerateEquivalentSphere(unsigned int num, float rad);

  void PathVariation();
    
  void CopyTotalCostToCostImage();

  void MapEndCostToSource();

  // this function outputs at every point the vector that 
  // points to the deepest part of the surface
  FieldTypePointer MapToVectorField();

protected: 
  GeodesicPaths()
    { 
      this->m_LocalCostFP=NULL;   
      this->m_SurfaceLabel=2; 
      this->m_MeanGeodesic=0;
      for (int i=0; i < GraphDimension; i++) this->m_Radius[i]=1;
    }
  ~GeodesicPaths(){}

private:
  GeodesicPaths(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented  
  CostImageIndexType  m_CostIndex; 
  typename GraphCostImageType::Pointer m_GraphCostImage;  /* if we are searching an image, set this */
  PixelType (*m_LocalCostFP)(Pointer);  // local cost function pointer 
  GraphCostImagePixelType    m_SurfaceLabel;
  GraphCostImagePixelType    m_SourceLabel;

  double                     m_MeanGeodesic;
  double                     m_SphereGeodesic;


};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGeodesicPaths.cxx"
#endif

#endif
