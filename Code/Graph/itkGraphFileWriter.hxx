/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGraphFileWriter.hxx,v $
  Language:  C++
  Date:      $Date: 2009/03/04 23:10:58 $
  Version:   $Revision: 1.18 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGraphFileWriter_hxx
#define __itkGraphFileWriter_hxx

#include "itkGraphFileWriter.h"
#include <fstream>

namespace itk {

//
// Constructor
//
template<class TInputGraph, class TInputImage>
GraphFileWriter<TInputGraph,TInputImage>
::GraphFileWriter()
{
  this->m_Input = NULL;
  this->m_FileName = "";
}

//
// Destructor
//
template<class TInputGraph, class TInputImage>
GraphFileWriter<TInputGraph,TInputImage>
::~GraphFileWriter()
{
}

//
// Set the input mesh
//
template<class TInputGraph, class TInputImage>
void
GraphFileWriter<TInputGraph,TInputImage>
::SetInput(InputGraphType * input)
{
  this->m_Input = input;
}

//
// Write the input mesh to the output file
//
template<class TInputGraph, class TInputImage>
void GraphFileWriter<TInputGraph,TInputImage>
::Update()
{
  this->GenerateData();
}

//
// Write the input mesh to the output file
//
template<class TInputGraph, class TInputImage>
void GraphFileWriter<TInputGraph,TInputImage>
::Write()
{
  this->GenerateData();
}

template<class TInputGraph, class TInputImage>
void
GraphFileWriter<TInputGraph,TInputImage>
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
  std::ofstream outputFile( m_FileName.c_str() );

  if( !outputFile.is_open() )
    {
    itkExceptionMacro( "Unable to open file\n"
        "outputFilename= " << m_FileName );
    return;
    }
  else
    {
    outputFile.close();
    }

  /**
   * Get filename extension
   */

  std::string::size_type pos = this->m_FileName.rfind( "." );
  std::string extension( this->m_FileName, pos+1, this->m_FileName.length()-1 );

  if( extension == "vtk" )
    {
    this->WriteVtkPolyData();
    }
  else if(( extension == "dot" ) || (extension == "gv") )
    {
    this->WriteDot();
    }  
  else
    {
    itkExceptionMacro( "Unknown extension: " << extension );
    }
}

template<class TInputGraph, class TInputImage>
void
GraphFileWriter<TInputGraph,TInputImage>
::WriteVtkPolyData()
{

  // convert to itk mesh
  typename MeshType::Pointer mesh = MeshType::New();
  mesh->GetPoints()->Reserve( this->m_Input->GetTotalNumberOfNodes() );

  for (unsigned int i=0; i<this->m_Input->GetTotalNumberOfNodes(); i++)
    { 
    PointType meshPoint;
        
    // FIXME - add option to pass reference image for conversion to spatial points
    for (unsigned int j=0; j<this->m_Input->GetNode(i).ImageIndex.GetIndexDimension(); j++)
      {
      meshPoint[j] = this->m_Input->GetNode(i).ImageIndex[j];
      }

    mesh->SetPoint( this->m_Input->GetNode(i).Identifier, meshPoint );
    
    }
    
    
  typename VtkWriterType::LineSetType::Pointer lines = VtkWriterType::LineSetType::New();
  lines->Reserve( this->m_Input->GetTotalNumberOfEdges() );
  for (unsigned int i=0; i<this->m_Input->GetTotalNumberOfEdges(); i++)
    { 
    typename VtkWriterType::LineType line(2);
    line[0] = this->m_Input->GetEdge(i).SourceIdentifier;
    line[1] = this->m_Input->GetEdge(i).TargetIdentifier;
    lines->InsertElement( i, line );
    }
    
    
  typename VtkWriterType::Pointer writer = VtkWriterType::New();
  writer->SetFileName( this->m_FileName );
  writer->SetInput( mesh );
  writer->SetLines( lines );
  writer->Update();
  
}

template<class TInputGraph, class TInputImage>
void
GraphFileWriter<TInputGraph,TInputImage>
::WriteDot()
{
  std::ofstream outputFile( m_FileName.c_str() );
  outputFile << "graph " << m_FileName << " {" << std::endl;
  for (unsigned int i=0; i<m_Input->GetTotalNumberOfEdges(); i++)
    {
    unsigned int a = m_Input->GetEdge(i).SourceIdentifier;
    unsigned int b = m_Input->GetEdge(i).TargetIdentifier;
    outputFile << "N" << a << " -- N" << b << ";" << std::endl;
    }
  outputFile << "}" << std::endl;
  outputFile.close();
}



template<class TInputGraph, class TInputImage>
void
GraphFileWriter<TInputGraph,TInputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "FileName: " << this->m_FileName << std::endl;
}


} //end of namespace itk

#endif
