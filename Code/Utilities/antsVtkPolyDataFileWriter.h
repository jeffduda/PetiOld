/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: antsVtkPolyDataFileWriter.h,v $
  Language:  C++
  Date:      $Date: 2009/03/04 23:10:58 $
  Version:   $Revision: 1.17 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __antsVtkPolyDataFileWriter_h
#define __antsVtkPolyDataFileWriter_h

#include "itkMesh.h"

#include "itkArray.h"
#include "itkImage.h"
#include "itkVectorContainer.h"

namespace itk {
namespace ants {

/** \class VtkPolyDataFileWriter
 * \brief
 * Writes an itkMesh to a file in various txt file formats.
 *
 */
template <class TInputMesh>
class VtkPolyDataFileWriter : public Object
{
public:
  /** Standard "Self" typedef. */
  typedef VtkPolyDataFileWriter               Self;
  typedef Object                                  Superclass;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;

  /** Method for creation through the object factory */
  itkNewMacro( Self );

  /** Write the Input mesh to the Output file.
   * Use either Update() or Write(). */
  void Update( void );
  void Write( void );

  /** Extract dimension from the output mesh. */
  itkStaticConstMacro( Dimension, unsigned int,
                       TInputMesh::PointType::Dimension );


  /** Run-time type information (and related methods). */
  itkTypeMacro( VtkPolyDataFileWriter, Object );

  /** Hold on to the type information specified by the template parameters. */
  typedef TInputMesh                             InputMeshType;
  typedef typename TInputMesh::Pointer           InputMeshPointer;
  typedef typename InputMeshType::MeshTraits     MeshTraits;
  typedef typename InputMeshType::Superclass     PointSetType;
  typedef typename InputMeshType::PointType      PointType;
  typedef typename MeshTraits::PixelType         PixelType;
  typedef Array<PixelType>                       MultiComponentScalarType;
  typedef Array<unsigned long>                   LineType;
  typedef VectorContainer<long, 
    MultiComponentScalarType>                    MultiComponentScalarSetType; 
  typedef VectorContainer<long,
    SmartPointer<MultiComponentScalarSetType> >   MultiComponentScalarMultiSetType;  
  typedef VectorContainer<long, std::string >     MultiComponentScalarSetNamesType;
    
  typedef VectorContainer<long, LineType>        LineSetType; 
  typedef Image<PixelType,
    itkGetStaticConstMacro( Dimension )>         VtkPolyDataImageType;
  typedef typename
    VtkPolyDataImageType::SizeType           ImageSizeType;
  typedef typename
    VtkPolyDataImageType::PointType          ImageOriginType;
  typedef typename
    VtkPolyDataImageType::SpacingType        ImageSpacingType;
  typedef typename
    VtkPolyDataImageType::DirectionType      ImageDirectionType;


  /** Set the Input */
  void SetInput( InputMeshType * input );

  /** Set/Get the name of the file where data are written. */
  itkSetStringMacro( FileName );
  itkGetStringMacro( FileName );

  /** Specify other attributes */
  itkSetMacro( Lines, typename LineSetType::Pointer );
  
  itkSetMacro( Polygons, typename LineSetType::Pointer );

  itkSetMacro( MultiComponentScalarSets, 
    typename MultiComponentScalarMultiSetType::Pointer );
    
  itkSetMacro( MultiComponentScalarSetNames,
    typename MultiComponentScalarSetNamesType::Pointer );

  /** Specify image attributes if output is an image. */
  itkSetMacro( ImageSize, ImageSizeType );
  itkGetConstMacro( ImageSize, ImageSizeType );

  itkSetMacro( ImageOrigin, ImageOriginType );
  itkGetConstMacro( ImageOrigin, ImageOriginType );

  itkSetMacro( ImageSpacing, ImageSpacingType );
  itkGetConstMacro( ImageSpacing, ImageSpacingType );

  itkSetMacro( ImageDirection, ImageDirectionType );
  itkGetConstMacro( ImageDirection, ImageDirectionType );

protected:
  VtkPolyDataFileWriter();
  virtual ~VtkPolyDataFileWriter();

  virtual void GenerateData();

  std::string                         m_FileName;
  InputMeshPointer                    m_Input;

  typename MultiComponentScalarMultiSetType::Pointer   m_MultiComponentScalarSets;
  typename LineSetType::Pointer                        m_Lines;
  typename LineSetType::Pointer                        m_Polygons;
  
  /**
   * If output is an image type, the attributes must be specified.
   */
  ImageSizeType                       m_ImageSize;
  ImageSpacingType                    m_ImageSpacing;
  ImageOriginType                     m_ImageOrigin;
  ImageDirectionType                  m_ImageDirection;
  
  typename MultiComponentScalarSetNamesType::Pointer m_MultiComponentScalarSetNames;

  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  VtkPolyDataFileWriter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  void WritePointsToAvantsFile();
  void WritePointsToImageFile();


  void WriteVTKFile();
  void WritePointsToVTKFile();
  void WritePointDataToVTKFile();
  void WriteCellDataToVTKFile();
  void WriteLinesToVTKFile();
  
  void WriteVTKPolygons();

};

}
} // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "antsVtkPolyDataFileWriter.txx"
#endif

#endif
