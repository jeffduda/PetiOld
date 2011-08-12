/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkZaleskyDiffusionTensorConnectivityGraphFunctor.h,v $
  Language:  C++
  Date:      $Date: 2008/11/11 11:51:34 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkZaleskyDiffusionTensorConnectivityGraphFunctor_h
#define __itkZaleskyDiffusionTensorConnectivityGraphFunctor_h

#include "itkProcessObject.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkArray.h"
#include "itkDefaultImageToGraphFunctor.h"

#include <list>

namespace itk
{

/** \class ZaleskyDiffusionTensorConnectivityGraphFunctor
 * \brief Class which defines node/edge weighting in constructing a
 *        graph from an image.
 **/
template<typename TInputImage, typename TOutputGraph>
class ZaleskyDiffusionTensorConnectivityGraphFunctor
: public ImageToGraphFunctor<TInputImage, TOutputGraph>
{
public:
  /** Standard class typedefs. */
  typedef ZaleskyDiffusionTensorConnectivityGraphFunctor Self;
  typedef ImageToGraphFunctor<TInputImage, TOutputGraph> Superclass;
  typedef SmartPointer<Self>                             Pointer;
  typedef SmartPointer<const Self>                       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ZaleskyDiffusionTensorConnectivityGraphFunctor, ImageToGraphFunctor);

  typedef TOutputGraph OutputGraphType;

  typedef typename Superclass::IndexType            IndexType;
  typedef typename Superclass::NodeType             NodeType;
  typedef typename Superclass::EdgeType             EdgeType;
  typedef typename Superclass::NodeIteratorType     NodeIteratorType;
  typedef typename Superclass::EdgeIteratorType     EdgeIteratorType;
  typedef typename Superclass::NodeWeightType       NodeWeightType;
  typedef typename Superclass::EdgeWeightType       EdgeWeightType;
  typedef typename Superclass::NodeImageType        NodeImageType;
  typedef typename Superclass::NodePointerType      NodePointerType;
  typedef typename Superclass::EdgePointerType      EdgePointerType;
  typedef typename Superclass::EdgeIdentifierContainerType
                                                    EdgeIdentifierContainerType;

  virtual bool IsPixelANode(IndexType idx)
    { return ( !this->m_ExcludeBackground ||
      ( this->GetInput()->GetPixel( idx ) != this->m_BackgroundValue ) ); }
  virtual EdgeWeightType GetEdgeWeight(IndexType idx1, IndexType idx2 )
      { return ( static_cast<EdgeWeightType>( 1 ) ); }
  virtual NodeWeightType GetNodeWeight( IndexType idx )
      { return ( static_cast<NodeWeightType>( 1 ) ); }
  virtual void NormalizeGraph( NodeImageType *im, OutputGraphType *g ) {}
  
protected:
  ZaleskyDiffusionTensorConnectivityGraphFunctor() {}
  ~ZaleskyDiffusionTensorConnectivityGraphFunctor() {}
  void PrintSelf(std::ostream& os, Indent indent) const
     { Superclass::PrintSelf( os, indent ); }

private:
  ZaleskyDiffusionTensorConnectivityGraphFunctor(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk

//#ifndef ITK_MANUAL_INSTANTIATION
//#include "itkZaleskyDiffusionTensorConnectivityGraphFunctor.hxx"
//#endif

#endif
