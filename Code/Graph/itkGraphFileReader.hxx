/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGraphFileReader.hxx,v $
  Language:  C++
  Date:      $Date: 2008/11/11 11:51:34 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGraphFileReader_hxx
#define __itkGraphFileReader_hxx

#include "itkGraphFileReader.h"

namespace itk
{

/**
 *
 */
template <class TInputImage, class TOutputGraph>
GraphFileReader<TInputImage, TOutputGraph>
::GraphFileReader()
{
  this->ProcessObject::SetNumberOfRequiredInputs( 0 );

  GraphPointer output
    = dynamic_cast<GraphType*>( this->MakeOutput( 0 ).GetPointer() );

  this->ProcessObject::SetNumberOfRequiredOutputs( 1 );
  this->ProcessObject::SetNthOutput( 0, output.GetPointer() );

}

/**
 *
 */
template <class TInputImage, class TOutputGraph>
GraphFileReader<TInputImage, TOutputGraph>
::~GraphFileReader()
{
}


/**
 *   Make Ouput
 */
template <class TInputImage, class TOutputGraph>
DataObject::Pointer
GraphFileReader<TInputImage, TOutputGraph>
::MakeOutput( unsigned int )
{
  GraphPointer  outputGraph = GraphType::New();
  return dynamic_cast< DataObject *>( outputGraph.GetPointer() );
}

template <class TInputImage, class TOutputGraph>
typename GraphFileReader<TInputImage, TOutputGraph>::GraphType *
GraphFileReader<TInputImage, TOutputGraph>
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
GraphFileReader<TInputImage, TOutputGraph>
::GenerateData()
{

  if( this->m_FileName == "" )
    {
    itkExceptionMacro( "No FileName" );
    return;
    }

  //
  // Read output file
  //
  std::ifstream inputFile( m_FileName.c_str() );

  if( !inputFile.is_open() )
    {
    itkExceptionMacro( "Unable to open file\n"
        "inputFilename= " << m_FileName );
    return;
    }
  else
    {
    inputFile.close();
    }

  /**
   * Get filename extension
   */
  std::string::size_type pos = this->m_FileName.rfind( "." );
  std::string extension( this->m_FileName, pos+1, this->m_FileName.length()-1 );

  if( extension == "vtk" )
    {
    this->ReadVtkPolyData();
    }
  else
    {
    itkExceptionMacro( "Unknown extension: " << extension );
    }
}

template <class TInputImage, class TOutputGraph>
void
GraphFileReader<TInputImage, TOutputGraph>
::ReadVtkPolyData()
{
  
  ImageConstPointer input;
  input = static_cast<const ImageType *>( this->GetInput( 0 ) );
  GraphPointer output = this->GetOutput();
  
  typename VtkReaderType::Pointer reader = VtkReaderType::New();
  reader->SetFileName( this->m_FileName );
  reader->Update();
  
  for (unsigned int i=0; i<reader->GetOutput()->GetNumberOfPoints(); i++)
    {
    NodePointerType node = output->CreateNewNode();
    for (unsigned int j=0; j<ImageType::ImageDimension; j++)
      {
      node->ImageIndex[j] = reader->GetOutput()->GetPoint( i )[j];
      }
      
    //if (nodeWeights != NULL)
    //  {
    //  node->Weight = nodeWeights->GetComponent(i,0);
    //  }
    
    }
  
  typename VtkReaderType::LineSetType::Pointer lines = reader->GetLines();  
  for (unsigned int j=0; j<lines->Size(); j++)
    {
    typename VtkReaderType::LineType line = lines->GetElement(j);
    if (line.Size() == 2)
      {
      EdgePointerType edge = output->CreateNewEdge();
      edge->SourceIdentifier = line[0];
      edge->TargetIdentifier = line[1];
      }
      
    //if (edgeWeights != NULL)
    //  {
    //  edge->Weight = edgeWeights->GetComponent(j,0);
    //  }  
      
    }

}

/**
 *
 */
template <class TInputImage, class TOutputGraph>
void
GraphFileReader<TInputImage, TOutputGraph>
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
GraphFileReader<TInputImage, TOutputGraph>
::GenerateOutputInformation()
{
}


} // end namespace itk

#endif
