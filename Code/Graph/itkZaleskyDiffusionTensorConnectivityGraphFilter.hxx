/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkZaleskyDiffusionTensorConnectivityGraphFilter.hxx,v $
  Language:  C++
  Date:      $Date: 2008/11/11 11:51:34 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkZaleskyDiffusionTensorConnectivityGraphFilter_hxx
#define __itkZaleskyDiffusionTensorConnectivityGraphFilter_hxx

#include "itkZaleskyDiffusionTensorConnectivityGraphFilter.h"
#include "itkDefaultImageToGraphFunctor.h"
#include "itkZaleskyDiffusionTensorConnectivityGraphFunctor.h"
#include "itkShapedNeighborhoodIterator.h"

namespace itk
{

/**
 *
 */
template <class TInputImage, class TTensorImage, class TOutputGraph>
ZaleskyDiffusionTensorConnectivityGraphFilter<TInputImage, TTensorImage, TOutputGraph>
::ZaleskyDiffusionTensorConnectivityGraphFilter()
{
  this->ProcessObject::SetNumberOfRequiredInputs( 2 );

  GraphPointer output
    = dynamic_cast<GraphType*>( this->MakeOutput( 0 ).GetPointer() );

  this->ProcessObject::SetNumberOfRequiredOutputs( 1 );
  this->ProcessObject::SetNthOutput( 0, output.GetPointer() );

  ZaleskyFunctorPointer ZaleskyGraphFunctor = ZaleskyFunctorType::New();
  this->m_ImageToGraphFunctor = ZaleskyGraphFunctor;
  
  /* Default to threshold of 20 degrees */
  this->m_MinimumAngleThreshold = 2*0.174532925;
  //this->m_MinimumAngleThreshold = 10.0;
}

/**
 *
 */
template <class TInputImage, class TTensorImage, class TOutputGraph>
ZaleskyDiffusionTensorConnectivityGraphFilter<TInputImage, TTensorImage, TOutputGraph>
::~ZaleskyDiffusionTensorConnectivityGraphFilter()
{
}

template <class TInputImage, class TTensorImage, class TOutputGraph>
void
ZaleskyDiffusionTensorConnectivityGraphFilter<TInputImage, TTensorImage, TOutputGraph>
::SetDiffusionTensorImage( TensorImageType* image )
{
  this->SetNthInput( 1, const_cast<TensorImageType *>( image ) );
}

template <class TInputImage, class TTensorImage, class TOutputGraph>
const typename ZaleskyDiffusionTensorConnectivityGraphFilter<TInputImage, TTensorImage, TOutputGraph>::TensorImageType *
ZaleskyDiffusionTensorConnectivityGraphFilter<TInputImage, TTensorImage, TOutputGraph>
::GetDiffusionTensorImage( void )
{
  return dynamic_cast<const TensorImageType*>
    ( this->ProcessObject::GetInput( 1 ) );
}

/**
 *
 */
template <class TInputImage, class TTensorImage, class TOutputGraph>
void
ZaleskyDiffusionTensorConnectivityGraphFilter<TInputImage, TTensorImage, TOutputGraph>
::GenerateData()
{
  ImageConstPointer input;
  input = static_cast<const ImageType *>( this->GetInput( 0 ) );

  this->m_ImageToGraphFunctor->SetInput( input );
  this->m_ImageToGraphFunctor->ActivateAllNeighbors();

  typedef Image<NodeIdentifierType, ImageType::ImageDimension> NodeImageType;
  typename NodeImageType::Pointer nodes = NodeImageType::New();
  nodes->SetRegions( input->GetBufferedRegion() );
  nodes->Allocate();

  ShapedNeighborhoodIterator<NodeImageType> It
    ( this->m_ImageToGraphFunctor->GetRadius(), nodes,
    nodes->GetBufferedRegion() );
  It.ClearActiveList();
  typename ImageToGraphFunctorType::IndexListType::const_iterator it;
  for( it = m_ImageToGraphFunctor->GetActiveIndexList().begin();
    it != m_ImageToGraphFunctor->GetActiveIndexList().end(); ++it )
    {
    It.ActivateOffset( It.GetOffset( *it ) );
    }

  GraphPointer output = this->GetOutput();

  /** Create the graph. */
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    NodePointerType node;
    IndexType idx = It.GetIndex();
    if( this->m_ImageToGraphFunctor->IsPixelANode( idx ) )
      {
      node = output->CreateNewNode(
        this->m_ImageToGraphFunctor->GetNodeWeight( idx ) );
      //node->ImageIndex = idx;
      for (unsigned int i=0; i<ImageType::ImageDimension; i++)
        {
        node->ImageIndex[i] = idx[i];  
        }
      It.SetCenterPixel( node->Identifier );
      }
    }

  /** Create the edges between neighboring pixels. */
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    IndexType idx = It.GetIndex();
    if ( this->m_ImageToGraphFunctor->IsPixelANode( idx ) )
      {
      TensorType tensor = this->GetDiffusionTensorImage()->GetPixel(idx);
      EigenValuesArrayType eigenValues;
      EigenVectorsMatrixType eigenVectors;
      tensor.ComputeEigenAnalysis(eigenValues, eigenVectors);
      
      typename ShapedNeighborhoodIterator<NodeImageType>::Iterator it;
      for( it = It.Begin(); !it.IsAtEnd(); it++ )
        {
        unsigned int i = it.GetNeighborhoodIndex();
        bool IsInBounds;
        It.GetPixel( i, IsInBounds );
        IndexType idx_i = It.GetIndex( i );
        if( !IsInBounds || idx == idx_i )
          {
          continue;
          }
            
        // If both nodes are in labeled white matter
        if( this->m_ImageToGraphFunctor->IsPixelANode( idx_i ) )
          {
          TensorType tensor_i = this->GetDiffusionTensorImage()->GetPixel(idx_i);
          EigenValuesArrayType eigenValues_i;
          EigenVectorsMatrixType eigenVectors_i;
          tensor_i.ComputeEigenAnalysis(eigenValues_i, eigenVectors_i);
               
          float minimumAngleCriteria = 0;
          for (unsigned int j=0; j<EigenVectorsMatrixType::ColumnDimensions; j++)
            {
            minimumAngleCriteria += eigenVectors[2][j] * eigenVectors_i[2][j];     
            }
          minimumAngleCriteria = fabs(minimumAngleCriteria);

          if ( vcl_acos(minimumAngleCriteria) < this->m_MinimumAngleThreshold )
            {           
            output->CreateNewEdge(It.GetCenterPixel(), It.GetPixel( i ), 1.0);
            }

          }
        }
      }      
    }
}

/**
 *
 */
template <class TInputImage, class TTensorImage, class TOutputGraph>
void
ZaleskyDiffusionTensorConnectivityGraphFilter<TInputImage, TTensorImage, TOutputGraph>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
}

/**
 * copy information from first input to all outputs
 * This is a void implementation to prevent the
 * ProcessObject version to be called
 */
template <class TInputImage, class TTensorImage, class TOutputGraph>
void
ZaleskyDiffusionTensorConnectivityGraphFilter<TInputImage, TTensorImage, TOutputGraph>
::GenerateOutputInformation()
{
}


} // end namespace itk

#endif
