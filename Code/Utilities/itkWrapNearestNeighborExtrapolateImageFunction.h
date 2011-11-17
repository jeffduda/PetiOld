/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkNearestNeighborExtrapolateImageFunction_h
#define __itkNearestNeighborExtrapolateImageFunction_h

#include "itkExtrapolateImageFunction.h"

namespace itk
{
/** \class NearestNeighborExtrapolateImageFunction
 * \brief Nearest neighbor extrapolation of a scalar image.
 *
 * NearestNeighborExtrapolateImageFunction extrapolate image intensity at
 * a specified point, continuous index or index by copying the intensity
 * of the nearest neighbor within the image buffer.
 *
 * This class is templated
 * over the input image type and the coordinate representation type
 * (e.g. float or double).
 *
 * \ingroup ImageFunctions
 * \ingroup ITKImageFunction
 */
template< class TInputImage, class TCoordRep = float >
class ITK_EXPORT WrapNearestNeighborExtrapolateImageFunction:
  public ExtrapolateImageFunction< TInputImage, TCoordRep >
{
public:
  /** Standard class typedefs. */
  typedef WrapNearestNeighborExtrapolateImageFunction        Self;
  typedef ExtrapolateImageFunction< TInputImage, TCoordRep > Superclass;
  typedef SmartPointer< Self >                               Pointer;
  typedef SmartPointer< const Self >                         ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(WrapNearestNeighborExtrapolateImageFunction,
               InterpolateImageFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** OutputType typedef support. */
  typedef typename Superclass::OutputType OutputType;

  /** InputImageType typedef support. */
  typedef typename Superclass::InputImageType InputImageType;

  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Index typedef support. */
  typedef typename Superclass::IndexType      IndexType;
  typedef typename IndexType::IndexValueType  IndexValueType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** Evaluate the function at a ContinuousIndex position
   *
   * Returns the extrapolated image intensity at a
   * specified position by returning the intensity of the
   * nearest neighbor within the image buffer.
   *
   */
  virtual OutputType EvaluateAtContinuousIndex(
    const ContinuousIndexType & index) const
  {

    std::cout << "EvaluateAtContinuousIndex" << std::endl;
    IndexType nindex;

    for ( unsigned int j = 0; j < ImageDimension; j++ )
      {
      nindex[j] = Math::RoundHalfIntegerUp< IndexValueType >(index[j]);
      int indexRange = this->GetEndIndex()[j] - this->GetStartIndex()[j] + 1;

      if ( nindex[j] < this->GetStartIndex()[j] )
        {
        std::cout << " too low " << std::endl;
        nindex[j] = this->GetEndIndex()[j] - ( nindex[j] % indexRange );
        }
      else if ( nindex[j] > this->GetEndIndex()[j] )
        {
        std::cout << " too high " << std::endl;
        nindex[j] = this->GetStartIndex()[j] + ( nindex[j] % indexRange );
        }
      }
    return static_cast< OutputType >( this->GetInputImage()->GetPixel(nindex) );
  }

  /** Evaluate the function at a ContinuousIndex position
   *
   * Returns the extrapolated image intensity at a
   * specified position by returning the intensity of the
   * nearest neighbor within the image buffer.
   *
   */
  virtual OutputType EvaluateAtIndex(
    const IndexType & index) const
  {
    std::cout << "EvaluateAtIndex" << std::endl;

    IndexType nindex;

    for ( unsigned int j = 0; j < ImageDimension; j++ )
      {
      int indexRange = this->GetEndIndex()[j] - this->GetStartIndex()[j] + 1;
      if ( index[j] < this->GetStartIndex()[j] )
        {
        nindex[j] = this->GetEndIndex()[j] - ( nindex[j] % indexRange );
        }
      else if ( index[j] > this->GetEndIndex()[j] )
        {
        nindex[j] = this->GetStartIndex()[j] + ( nindex[j] % indexRange );
        }
      else
        {
        nindex[j] = index[j];
        }
      }
    return static_cast< OutputType >( this->GetInputImage()->GetPixel(nindex) );
  }

protected:
  WrapNearestNeighborExtrapolateImageFunction(){}
  ~WrapNearestNeighborExtrapolateImageFunction(){}
  void PrintSelf(std::ostream & os, Indent indent) const
  { Superclass::PrintSelf(os, indent); }
private:
  WrapNearestNeighborExtrapolateImageFunction(const Self &); //purposely not
                                                         // implemented
  void operator=(const Self &);                          //purposely not

  // implemented
};
} // end namespace itk

#endif
