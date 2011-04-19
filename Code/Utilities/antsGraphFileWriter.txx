/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: antsGraphFileWriter.txx,v $
  Language:  C++
  Date:      $Date: 2009/03/04 23:10:58 $
  Version:   $Revision: 1.18 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __antsGraphFileWriter_txx
#define __antsGraphFileWriter_txx

#include "antsGraphFileWriter.h"

#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkTriangleFilter.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkFloatArray.h"

#include <fstream>

namespace itk {
namespace ants {

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
  vtkPoints * pointSet = vtkPoints::New();
  pointSet->Initialize();
  pointSet->SetNumberOfPoints( this->m_Input->GetTotalNumberOfNodes() );  
  
  //vtkPointData * nodeWeights = vtkPointData::New();
  vtkFloatArray * nodeWeights = vtkFloatArray::New();
  nodeWeights->SetNumberOfTuples(1);
  nodeWeights->SetNumberOfValues( this->m_Input->GetTotalNumberOfNodes() ); // this allocates memory
  nodeWeights->SetName( "node-weights" );
  
  for (unsigned int i=0; i<this->m_Input->GetTotalNumberOfNodes(); i++)
    {
    float pt[this->m_Input->GetNode(i).ImageIndex.GetIndexDimension()];
    
    for (unsigned int j=0; j<this->m_Input->GetNode(i).ImageIndex.GetIndexDimension(); j++)
      {
      pt[j] = this->m_Input->GetNode(i).ImageIndex[j];
      }
    pointSet->SetPoint(i,pt);
    nodeWeights->SetValue(i,this->m_Input->GetNode(i).Weight);
    }
    
  vtkCellArray * edges = vtkCellArray::New();
  edges->Initialize();
  
  vtkFloatArray * edgeWeights = vtkFloatArray::New();
  edgeWeights->SetNumberOfTuples(1);
  edgeWeights->SetNumberOfValues( this->m_Input->GetTotalNumberOfEdges() ); // this allocates memory
  edgeWeights->SetName( "edge-weights" );
  
  vtkIdType * pts = new vtkIdType [ 2 ];
  for (unsigned int i=0; i<m_Input->GetTotalNumberOfEdges(); i++)
    {
    pts[0] = this->m_Input->GetEdge(i).SourceIdentifier;
    pts[1] = this->m_Input->GetEdge(i).TargetIdentifier;
    edges->InsertNextCell(2,pts);
    edgeWeights->SetValue(i,this->m_Input->GetEdge(i).Weight);
    }
  delete [] pts;
  
  std::cout << nodeWeights->GetNumberOfTuples() << std::endl;
  
  vtkPolyData * outputData = vtkPolyData::New();
  outputData->Initialize();
  outputData->SetPoints( pointSet );
  outputData->SetLines( edges );
  outputData->GetPointData()->AddArray( nodeWeights );
  outputData->GetCellData()->AddArray( edgeWeights );
  
  vtkPolyDataWriter * writer = vtkPolyDataWriter::New();
  writer->SetInput( outputData );
  writer->SetFileName( this->m_FileName.c_str() );
  //writer->SetFileTypeToBinary();
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

}
} //end of namespace itk

#endif
