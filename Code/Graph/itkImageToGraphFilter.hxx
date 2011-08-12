/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImageToGraphFilter.hxx,v $
  Language:  C++
  Date:      $Date: 2008/11/11 11:51:34 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkImageToGraphFilter_hxx
#define __itkImageToGraphFilter_hxx

#include "itkImageToGraphFilter.h"
#include "itkDefaultImageToGraphFunctor.h"
#include "itkShapedNeighborhoodIterator.h"

namespace itk
{

/**
 *
 */
template <class TInputImage, class TOutputGraph>
ImageToGraphFilter<TInputImage, TOutputGraph>
::ImageToGraphFilter()
{
  this->ProcessObject::SetNumberOfRequiredInputs( 1 );

  GraphPointer output
    = dynamic_cast<GraphType*>( this->MakeOutput( 0 ).GetPointer() );

  this->ProcessObject::SetNumberOfRequiredOutputs( 1 );
  this->ProcessObject::SetNthOutput( 0, output.GetPointer() );

  typedef DefaultImageToGraphFunctor<ImageType,
    GraphType> DefaultImageToGraphFunctorType;
  typename DefaultImageToGraphFunctorType::Pointer DefaultImageToGraphFunctor
     = DefaultImageToGraphFunctorType::New();
  this->m_ImageToGraphFunctor = DefaultImageToGraphFunctor;
}

/**
 *
 */
template <class TInputImage, class TOutputGraph>
ImageToGraphFilter<TInputImage, TOutputGraph>
::~ImageToGraphFilter()
{
}


/**
 *   Make Ouput
 */
template <class TInputImage, class TOutputGraph>
DataObject::Pointer
ImageToGraphFilter<TInputImage, TOutputGraph>
::MakeOutput( unsigned int )
{
  GraphPointer  outputGraph = GraphType::New();
  return dynamic_cast< DataObject *>( outputGraph.GetPointer() );
}

template <class TInputImage, class TOutputGraph>
void
ImageToGraphFilter<TInputImage, TOutputGraph>
::SetInput( ImageType* image )
{
  this->SetNthInput( 0, const_cast<ImageType *>( image ) );
}

template <class TInputImage, class TOutputGraph>
void
ImageToGraphFilter<TInputImage, TOutputGraph>
::SetInput( unsigned int idx, const ImageType *input )
{
  // process object is not const-correct, the const_cast
  // is required here.
  this->ProcessObject::SetNthInput( idx,
    const_cast< ImageType * >( input ) );
}

template <class TInputImage, class TOutputGraph>
const typename ImageToGraphFilter<TInputImage, TOutputGraph>::ImageType *
ImageToGraphFilter<TInputImage, TOutputGraph>
::GetInput( unsigned int idx )
{
  return dynamic_cast<const ImageType*>
    ( this->ProcessObject::GetInput( idx ) );
}

template <class TInputImage, class TOutputGraph>
typename ImageToGraphFilter<TInputImage, TOutputGraph>::GraphType *
ImageToGraphFilter<TInputImage, TOutputGraph>
::GetOutput( void )
{
  return dynamic_cast<GraphType*>
    ( this->ProcessObject::GetOutput( 0 ) );
}


/**
 *
 */
template <class TInputImage, class TOutputGraph>
void
ImageToGraphFilter<TInputImage, TOutputGraph>
::GenerateData()
{
  ImageConstPointer input;
  input = static_cast<const ImageType *>( this->GetInput( 0 ) );

  this->m_ImageToGraphFunctor->SetInput( input );

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
    typename ShapedNeighborhoodIterator<NodeImageType>::Iterator localIt;
    
    if ( this->m_ImageToGraphFunctor->IsPixelANode( idx ) )
      {    

       for( localIt = It.Begin(); !localIt.IsAtEnd(); localIt++ )
         {

         unsigned int i = localIt.GetNeighborhoodIndex();
         bool IsInBounds;
         It.GetPixel( i, IsInBounds );       
         IndexType idx_i = It.GetIndex( i );

         if( this->m_ImageToGraphFunctor->IsPixelANode( idx ) &&
           this->m_ImageToGraphFunctor->IsPixelANode( idx_i ) && (idx_i != idx) && IsInBounds )
           {
           output->CreateNewEdge( It.GetCenterPixel(), It.GetPixel( i ),
           this->m_ImageToGraphFunctor->GetEdgeWeight( idx, idx_i ) );
           }
           
         }
       }
    }
  this->m_ImageToGraphFunctor->NormalizeGraph( nodes, output );
}

/**
 *
 */
template <class TInputImage, class TOutputGraph>
void
ImageToGraphFilter<TInputImage, TOutputGraph>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
}

/**
 * copy information from first input to all outputs
 * This is a void implementation to prevent the
 * ProcessObject version to be called
 */
template <class TInputImage, class TOutputGraph>
void
ImageToGraphFilter<TInputImage, TOutputGraph>
::GenerateOutputInformation()
{
}


} // end namespace itk

#endif
