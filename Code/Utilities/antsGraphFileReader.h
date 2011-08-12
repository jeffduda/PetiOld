/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: antsVtkPolyDataFileReader.h,v $
  Language:  C++
  Date:      $Date: 2009/03/04 23:10:58 $
  Version:   $Revision: 1.17 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __antsGraphFileReader_h
#define __antsGraphFileReader_h

#include "itkGraph.h"
#include "itkGraphSource.h"

#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkTriangleFilter.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataReader.h"

namespace itk {
namespace ants {

/** \class VtkPolyDataFileReader
 * \brief
 * Writes an itkMesh to a file in various txt file formats.
 *
 */
template < class TOutputGraph, class TOutputImage=Image<float,TOutputGraph::GraphTraitsType::ImageDimension> >
class GraphFileReader: public itk::GraphSource<TOutputGraph>
{
public:
  /** Standard "Self" typedef. */
  typedef GraphFileReader               Self;
  typedef Object                                  Superclass;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;

  /** Method for creation through the object factory */
  itkNewMacro( Self );

  /** Write the Input mesh to the Output file.
   * Use either Update() or Write(). */
  void Update( void );
  void Read( void );

  /** Extract dimension from the output mesh. */
  //itkStaticConstMacro( Dimension, unsigned int,
  //                     TInputMesh::PointType::Dimension );


  /** Run-time type information (and related methods). */
  itkTypeMacro( GraphFileReader, GraphSource );

  /** Hold on to the type information specified by the template parameters. */
  typedef TOutputGraph                              OutputGraphType;
  typedef typename TOutputGraph::Pointer            OutputGraphPointer;
  typedef typename OutputGraphType::GraphTraitsType GraphTraitsType;
  typedef typename GraphTraitsType::IndexType       IndexType;
  
  typedef TOutputImage                              OutputImageType;
  typedef typename OutputImageType::IndexType       ImageIndexType;
  
  
  OutputGraphType * GetOutput( void );
  OutputGraphType * GetOutput( unsigned int idx );

  /** Set/Get the name of the file where data are written. */
  itkSetStringMacro( FileName );
  itkGetStringMacro( FileName );

protected:
  GraphFileReader();
  virtual ~GraphFileReader();

  virtual void GenerateData();

  std::string                          m_FileName;

  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  GraphFileReader(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  void ReadVtkPolyData();

};

}
} // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "antsGraphFileReader.hxx"
#endif

#endif
