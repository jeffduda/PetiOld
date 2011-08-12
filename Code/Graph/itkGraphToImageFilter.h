/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGraphToImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/11/11 11:51:34 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGraphToImageFilter_h
#define __itkGraphToImageFilter_h

#include "itkImageSource.h"

namespace itk
{

/** \class GraphToImageFilter
 * \brief Base class for filters that take a Graph as input and
 *  produce an image as output.
 *
 * \par
 * Base class for filters that take a Graph as input and produce an image as
 * output. By default, if the user does not specify the size of the output
 * image, the maximum size of the graphs's bounding box is used.  The default
 * spacing is currently assumed internally to be 1.0.
 */
template <class TInputGraph, class TOutputImage>
class ITK_EXPORT GraphToImageFilter : public ImageSource<TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef GraphToImageFilter            Self;
  typedef ImageSource<TOutputImage>     Superclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( GraphToImageFilter, ImageSource );

  /** Some convenient typedefs. */
  typedef typename Superclass::OutputImageRegionType  OutputImageRegionType;
  typedef          TInputGraph                        InputGraphType;
  typedef typename InputGraphType::Pointer            InputGraphPointer;
  typedef typename InputGraphType::ConstPointer       InputGraphConstPointer;
  typedef          TOutputImage                       OutputImageType;
  typedef typename OutputImageType::Pointer           OutputImagePointer;
  typedef typename OutputImageType::SizeType          SizeType;
  typedef typename OutputImageType::ValueType         ValueType;

  /** ImageDimension constants */
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Set/Get the graph input of this process object.  */
  virtual void SetInput(const InputGraphType *graph);
  virtual void SetInput(unsigned int, const TInputGraph *graph);
  const InputGraphType * GetInput(void);
  const InputGraphType * GetInput(unsigned int idx);

  /** Spacing (size of a pixel) of the output image. The
   * spacing is the geometric distance between image samples.
   * It is stored internally as double, but may be set from
   * float. \sa GetSpacing() */
  virtual void SetSpacing( const double spacing[OutputImageDimension] );
  virtual void SetSpacing( const float spacing[OutputImageDimension] );
  virtual const double* GetSpacing() const;

  /** Set/Get the value for pixels on and off the path.
  * By default, this filter will return a "0" image with path pixels set to 1 */
  itkSetMacro( BackgroundValue, ValueType );
  itkGetMacro( BackgroundValue, ValueType );

  /** The origin of the output image. The origin is the geometric
   * coordinates of the index (0,0,...,0).  It is stored internally
   * as double but may be set from float.
   * \sa GetOrigin() */
  virtual void SetOrigin( const double origin[OutputImageDimension] );
  virtual void SetOrigin( const float origin[OutputImageDimension] );

  /** Set/Get Size */
  itkSetMacro( Size, SizeType );
  itkGetMacro( Size, SizeType );

protected:
  GraphToImageFilter();
  ~GraphToImageFilter();

  virtual void GenerateOutputInformation(){}; // do nothing
  virtual void GenerateData();


  SizeType     m_Size;
  double       m_Spacing[OutputImageDimension];
  double       m_Origin[OutputImageDimension];
  ValueType    m_BackgroundValue;

  virtual void PrintSelf( std::ostream& os, Indent indent ) const;

private:
  GraphToImageFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented


};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGraphToImageFilter.hxx"
#endif

#endif
