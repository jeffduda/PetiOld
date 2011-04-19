/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGraphFileReader.h,v $
  Language:  C++
  Date:      $Date: 2008/11/11 11:51:34 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGraphFileReader_h
#define __itkGraphFileReader_h

#include "itkGraphSource.h"
#include "itkObjectFactory.h"
#include "itkDefaultImageToGraphFunctor.h"

namespace itk
{

/** \class GraphFileReader
 * \brief 
 *
 * GraphFileReader is the base class for all process objects that output
 * Graph data and require image data as input. Specifically, this class
 * defines the SetInput() method for defining the input to a filter.
 *
 * \ingroup ImageFilters
 */
template < class TOutputGraph, class TInputImage=Image<float,TOutputGraph::GraphTraitsType::ImageDimension> >
class ITK_EXPORT GraphFileReader : public GraphSource<TOutputGraph>
{
public:
  /** Standard class typedefs. */
  typedef GraphFileReader            Self;
  typedef GraphSource<TOutputGraph>     Superclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;
  
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro( GraphFileReader, GraphSource );

  /** Create a valid output. */
  DataObject::Pointer  MakeOutput( unsigned int idx );

  /** Some Image related typedefs. */
  typedef   TInputImage                        ImageType;
  typedef   typename ImageType::Pointer        ImagePointer;
  typedef   typename ImageType::ConstPointer   ImageConstPointer;
  typedef   typename ImageType::RegionType     RegionType; 
  typedef   typename ImageType::PixelType      PixelType; 
  typedef   typename ImageType::IndexType      IndexType; 

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

  /** Set/Get the name of the file where data are written. */
  itkSetStringMacro( FileName );
  itkGetStringMacro( FileName );
  
  /** Get the output Graph of this process object.  */
  GraphType * GetOutput( void );
  
  /** Prepare the output */
  void GenerateOutputInformation( void );



protected:
  GraphFileReader();
  ~GraphFileReader();
  void PrintSelf( std::ostream& os, Indent indent ) const;
  
  void GenerateData();

  void ReadVtkPolyData();
  
  std::string m_FileName;
 
private:
  GraphFileReader( const GraphFileReader& ); //purposely not implemented
  void operator=( const GraphFileReader& ); //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGraphFileReader.txx"
#endif

#endif
