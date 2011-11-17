#include <limits.h>
#include <iomanip>

#include "itkVectorImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTranslationTransform.h"
#include "itkTransformFileReader.h"
#include "itkResampleImageFilter.h"
#include "itkWrapNearestNeighborExtrapolateImageFunction.h"
#include "itkTransformFactory.h"

/* 
 * Main progam
 */
int main( int argc, char *argv[] )
{
  
  char * inputfile = argv[1];
  char * transformfile = argv[2];
  char * outputfile = argv[3];

  typedef itk::VectorImage< int, 3 > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typedef itk::TransformFileReader TransformReaderType;

  typedef itk::TranslationTransform< double, 3> TransformType;
  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
  typedef itk::WrapNearestNeighborExtrapolateImageFunction<ImageType, double> ExtrapolatorType;
  
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName (inputfile);
  reader->Update();

  std::cout << "read image" << std::endl;

  std::cout << "reading " << transformfile << std::endl;
  itk::TransformFactory<TransformType>::RegisterTransform();
  TransformReaderType::Pointer tReader = TransformReaderType::New();
  tReader->SetFileName( transformfile );
  tReader->Update();
  TransformType::Pointer transform = 
    dynamic_cast< TransformType* >( ( tReader->GetTransformList() )->front().GetPointer() );
     
  std::cout << transform << std::endl;

  ExtrapolatorType::Pointer wrap = ExtrapolatorType::New();
  ResampleFilterType::Pointer resample = ResampleFilterType::New();
  resample->SetInput( reader->GetOutput() );
  resample->SetOutputParametersFromImage( reader->GetOutput() );
  resample->SetExtrapolator( wrap );
  resample->SetTransform( transform );
  try
    {
    std::cout << "resample" << std::endl;
    resample->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Caught exception:" << std::endl << e << std::endl;
    }
  
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputfile );
  writer->SetInput( resample->GetOutput() );
  writer->Update();
} 
