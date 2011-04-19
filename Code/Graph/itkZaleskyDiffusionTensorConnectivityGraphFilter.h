/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkZaleskyDiffusionTensorConnectivityGraphFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/11/11 11:51:34 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkZaleskyDiffusionTensorConnectivityGraphFilter_h
#define __itkZaleskyDiffusionTensorConnectivityGraphFilter_h

#include "itkGraphSource.h"
#include "itkObjectFactory.h"
#include "itkDefaultImageToGraphFunctor.h"
#include "itkZaleskyDiffusionTensorConnectivityGraphFunctor.h"

namespace itk
{

/** \class ZaleskyDiffusionTensorConnectivityGraphFilter
 * \brief 
 *
 * ZaleskyDiffusionTensorConnectivityGraphFilter is the base class for all process objects that output
 * Graph data and require image data as input. Specifically, this class
 * defines the SetInput() method for defining the input to a filter.
 *
 * \ingroup ImageFilters
 */
template <class TInputImage, class TTensorImage, class TOutputGraph>
class ITK_EXPORT ZaleskyDiffusionTensorConnectivityGraphFilter : public ImageToGraphFilter<TInputImage,TOutputGraph>
{
public:
  /** Standard class typedefs. */
  typedef ZaleskyDiffusionTensorConnectivityGraphFilter            Self;
  typedef ImageToGraphFilter<TInputImage,TOutputGraph>     Superclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;
  
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro( ZaleskyDiffusionTensorConnectivityGraphFilter, ImageToGraphFilter );

  /** Some Image related typedefs. */
  typedef   TInputImage                        ImageType;
  typedef   typename ImageType::Pointer        ImagePointer;
  typedef   typename ImageType::ConstPointer   ImageConstPointer;
  typedef   typename ImageType::RegionType     RegionType; 
  typedef   typename ImageType::PixelType      PixelType; 
  typedef   typename ImageType::IndexType      IndexType; 
  
  typedef   TTensorImage                       TensorImageType;
  typedef typename TensorImageType::PixelType  TensorType;
  typedef typename TensorType::EigenValuesArrayType EigenValuesArrayType;
  typedef typename TensorType::EigenVectorsMatrixType EigenVectorsMatrixType;
  
  /** Some Graph related typedefs. */
  typedef TOutputGraph                              GraphType;
  typedef typename GraphType::GraphTraitsType       GraphTraitsType; 
  typedef typename GraphType::Pointer               GraphPointer;
  typedef typename GraphType::NodeType              NodeType; 
  typedef typename GraphType::NodePointerType       NodePointerType; 
  typedef typename GraphType::NodeIdentifierType    NodeIdentifierType; 
  typedef typename GraphTraitsType::NodeWeightType  NodeWeightType;
  typedef typename GraphType::EdgeType              EdgeType; 
  typedef typename GraphType::EdgePointerType       EdgePointerType; 
  typedef typename GraphType::EdgeIdentifierType    EdgeIdentifierType; 
  typedef typename GraphTraitsType::EdgeWeightType  EdgeWeightType;
  
  /** Abstract ZaleskyDiffusionTensorConnectivityGraphFunctorType */
  typedef ImageToGraphFunctor<ImageType, GraphType> ImageToGraphFunctorType;
  typedef typename ImageToGraphFunctorType::Pointer ImageToGraphFunctorPointer;
  
  typedef ZaleskyDiffusionTensorConnectivityGraphFunctor<ImageType, GraphType> ZaleskyFunctorType;
  typedef typename ZaleskyFunctorType::Pointer ZaleskyFunctorPointer;
  
  /** Set/Get ZaleskyDiffusionTensorConnectivityGraphFunctor */
  itkGetObjectMacro( ImageToGraphFunctor, ImageToGraphFunctorType );
  itkSetObjectMacro( ImageToGraphFunctor, ImageToGraphFunctorType );
  
  /* Set/Get MinimumAngleThreshold for defining connectivity between adjacent tensors */
  itkGetMacro( MinimumAngleThreshold, float );
  itkSetMacro( MinimumAngleThreshold, float );

  /** Set the tensor input image */
  void SetDiffusionTensorImage( TensorImageType * );
  const TensorImageType * GetDiffusionTensorImage( void );

  /** Prepare the output */
  void GenerateOutputInformation( void );

protected:
  ZaleskyDiffusionTensorConnectivityGraphFilter();
  ~ZaleskyDiffusionTensorConnectivityGraphFilter();
  void PrintSelf( std::ostream& os, Indent indent ) const;
  
  void GenerateData();
 
private:
  ZaleskyDiffusionTensorConnectivityGraphFilter( const ZaleskyDiffusionTensorConnectivityGraphFilter& ); //purposely not implemented
  void operator=( const ZaleskyDiffusionTensorConnectivityGraphFilter& ); //purposely not implemented

  ImageToGraphFunctorPointer m_ImageToGraphFunctor;
  
  float m_MinimumAngleThreshold;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkZaleskyDiffusionTensorConnectivityGraphFilter.txx"
#endif

#endif
