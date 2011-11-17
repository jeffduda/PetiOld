#include <limits.h>
#include <iomanip>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkTranslationTransform.h"
#include "itkTransformFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkWrapNearestNeighborExtrapolateImageFunction.h"

/* 
 * Main progam
 */
int main( int argc, char *argv[] )
{
  
  char * inputfile = argv[1];
  unsigned int wrapDim = atoi( argv[2] );
  char * outputfile = argv[3];

  typedef itk::Image< int, 3 > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typedef itk::TransformFileWriter TransformWriterType;

  typedef itk::ConnectedComponentImageFilter< ImageType, ImageType > ConnectedFilterType;
  typedef itk::RelabelComponentImageFilter< ImageType, ImageType > RelabelFilterType;
  typedef itk::LabelStatisticsImageFilter< ImageType, ImageType > LabelStatsFilter;

  typedef itk::TranslationTransform< double, 3> TransformType;
  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
  typedef itk::WrapNearestNeighborExtrapolateImageFunction<ImageType, double> ExtrapolatorType;
  

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName (inputfile);
  reader->Update();

  ImageType::RegionType region = reader->GetOutput()->GetLargestPossibleRegion();
  ImageType::RegionType::SizeType size = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
  unsigned long origDim = size[wrapDim];
  size[wrapDim] = 2*size[wrapDim];

  std::cout << "Wrapping dimension = " << wrapDim << std::endl;
  std::cout << "Switching to size = " << size << std::endl;

  ImageType::Pointer wideImage = ImageType::New();
  region.SetSize( size );
  wideImage->SetRegions( region );
  wideImage->SetDirection( reader->GetOutput()->GetDirection() );
  wideImage->SetSpacing( reader->GetOutput()->GetSpacing() );
  wideImage->Allocate();

  //ImageType::PointType origin = reader->GetOutput()->GetOrigin();
  //origin[wrapDim] = origin[wrapDim] - (origDim/2)*wideImage->GetSpacing()[wrapDim];
  //wideImage->SetOrigin( origin );

  itk::ImageRegionIteratorWithIndex<ImageType> iter( wideImage, wideImage->GetLargestPossibleRegion() );
  while ( !iter.IsAtEnd() )
    {
    ImageType::IndexType idx = iter.GetIndex();
    ImageType::IndexType origIdx = idx;

    int wrapIdx = idx[wrapDim] - origDim/2;
    if (wrapIdx < 0)
      {
      wrapIdx = origDim + wrapIdx;
      }
    else if ( wrapIdx >= origDim )
      {
      wrapIdx = wrapIdx - origDim;
      }
    
    origIdx[wrapDim] = wrapIdx;
    iter.Set( reader->GetOutput()->GetPixel( origIdx ) );
    ++iter;
    }

  ConnectedFilterType::Pointer filter = ConnectedFilterType::New();
  RelabelFilterType::Pointer relabel = RelabelFilterType::New();

  filter->SetInput( wideImage );
  filter->SetFullyConnected( 0 );
  relabel->SetInput( filter->GetOutput() );
  relabel->SetMinimumObjectSize( 20 );
  relabel->Update();

  itk::ImageRegionIteratorWithIndex<ImageType> iter2( relabel->GetOutput(), relabel->GetOutput()->GetLargestPossibleRegion() );
  bool foundSubRegions = false;

  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();

  while ( !iter2.IsAtEnd() )
    {
    if ( iter2.Value() > 1 )
      {
      iter2.Set( 0 );
      foundSubRegions = true;
      }
    ++iter2;
    }

  if (!foundSubRegions)
    {
    std::cout << "Object was continuous, unwrapping failed. Try eroding the mask first." << std::endl;
    }
  else
    {

    std::cout << "Get Label Stats" << std::endl;
    LabelStatsFilter::Pointer stats = LabelStatsFilter::New();
    stats->SetInput( relabel->GetOutput() );
    stats->SetLabelInput( relabel->GetOutput() );
    stats->Update();
    
    LabelStatsFilter::BoundingBoxType box = stats->GetBoundingBox( 1 );
    unsigned int low = box[2*wrapDim];
    unsigned int high = box[2*wrapDim + 1];
      
    unsigned int objSize = high-low+1;
    unsigned int totalmargin = origDim - objSize;
    unsigned int offset = totalmargin - totalmargin/2;

    if (objSize > origDim)
      {
      std::cout << "Object is bigger than original image space." << std::endl;
      }
    else
      {      
      int shift = low - origDim/2;
      std::cout << "offset = " << offset << ", shift = " << shift << std::endl; 

      TransformType::OutputVectorType offsetVector;
      offsetVector.Fill( 0 );
      offsetVector[wrapDim] = (shift - offset)*reader->GetOutput()->GetSpacing()[wrapDim];
      offsetVector = reader->GetOutput()->GetDirection() * offsetVector;

      transform->SetOffset( offsetVector );
      std::cout << transform << std::endl;

      ExtrapolatorType::Pointer wrap = ExtrapolatorType::New();

      ResampleFilterType::Pointer resample = ResampleFilterType::New();
      //resample->SetOutputSpacing( reader->GetOutput()->GetSpacing() );
      //resample->SetOutputDirection( reader->GetOutput()->GetDirection() );
      //resample->SetOutputOrigin( reader->GetOutput()->GetOrigin() );
      //resample->SetOutputSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
      resample->SetOutputParametersFromImage( reader->GetOutput() );
      resample->SetInput( reader->GetOutput() );
      resample->SetExtrapolator( wrap );
      resample->SetTransform( transform );

      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( outputfile );
      writer->SetInput( resample->GetOutput() );
      //writer->Update();
      
      }
    }
    
  TransformWriterType::Pointer transWriter = TransformWriterType::New();
  transWriter->SetFileName( outputfile );
  transWriter->SetInput( transform );
  transWriter->Update();

} 
