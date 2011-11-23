#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <utility>


#include "itkVectorImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTranslationTransform.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkWrapNearestNeighborExtrapolateImageFunction.h"
#include "itkTransformFactory.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"

#include "antsCommandLineParser.h"
#include "antsCommandLineOption.h"


void ConvertToLowerCase( std::string& str )
{
  std::transform( str.begin(), str.end(), str.begin(), tolower );
// You may need to cast the above line to (int(*)(int))
// tolower - this works as is on VC 7.1 but may not work on
// other compilers
} 

template <unsigned int Dimension, class PixelType>
int ImageTransformOps( itk::ants::CommandLineParser::OptionType *option,
  itk::ants::CommandLineParser::OptionType *outputOption = NULL )
{

  typedef itk::Image<PixelType, Dimension>      ImageType;
  typedef itk::ImageFileReader< ImageType >     ImageReaderType;
  typedef itk::ImageFileWriter< ImageType >     ImageWriterType;

  typedef itk::Image<unsigned long, Dimension>    LabelImageType;
  typedef itk::ImageFileReader< LabelImageType >  LabelReaderType;
  typedef itk::ImageFileWriter< LabelImageType >  LabelWriterType;  

  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
  typedef itk::TransformFileWriter TransformWriterType;

  if( option->GetNumberOfParameters( 0 ) < 1 )
    {
    std::cerr << "image-transform:  Incorrect number of parameters." << std::endl;
    return EXIT_FAILURE; 
    }

  std::string value = option->GetValue( 0 );
  ConvertToLowerCase( value );

  if( strcmp( value.c_str(), "getwrap" ) == 0 ) 
    {
      
    typedef itk::TranslationTransform< double, Dimension> TransformType;
    typedef itk::ConnectedComponentImageFilter< LabelImageType, LabelImageType > ConnectedFilterType;
    typedef itk::RelabelComponentImageFilter< LabelImageType, LabelImageType > RelabelFilterType;
    typedef itk::LabelStatisticsImageFilter< LabelImageType, LabelImageType > LabelStatsFilter;
    typedef itk::WrapNearestNeighborExtrapolateImageFunction<LabelImageType, double> ExtrapolatorType;
    
    std::string imageName = option->GetParameter( 0 );
    unsigned int wrapDim = option->Convert<unsigned int>( option->GetParameter( 1 ) );

    typename LabelReaderType::Pointer reader = LabelReaderType::New();
    reader->SetFileName (imageName);
    reader->Update();
    
    typename LabelImageType::RegionType region = reader->GetOutput()->GetLargestPossibleRegion();
    typename LabelImageType::RegionType::SizeType size = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
    unsigned long origDim = size[wrapDim];
    size[wrapDim] = 2*size[wrapDim];

    typename LabelImageType::Pointer wideImage = LabelImageType::New();
    region.SetSize( size );
    wideImage->SetRegions( region );
    wideImage->SetDirection( reader->GetOutput()->GetDirection() );
    wideImage->SetSpacing( reader->GetOutput()->GetSpacing() );
    wideImage->Allocate();   
    
    typename itk::ImageRegionIteratorWithIndex<LabelImageType> iter( wideImage, wideImage->GetLargestPossibleRegion() );
    while ( !iter.IsAtEnd() )
      {
      typename LabelImageType::IndexType idx = iter.GetIndex();
      typename LabelImageType::IndexType origIdx = idx;
      
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
    
    typename ConnectedFilterType::Pointer filter = ConnectedFilterType::New();
    typename RelabelFilterType::Pointer relabel = RelabelFilterType::New();
    
    filter->SetInput( wideImage );
    filter->SetFullyConnected( 0 );
    relabel->SetInput( filter->GetOutput() );
    relabel->SetMinimumObjectSize( 20 );
    relabel->Update();
    
    typename itk::ImageRegionIteratorWithIndex<LabelImageType> iter2( relabel->GetOutput(), relabel->GetOutput()->GetLargestPossibleRegion() );
    bool foundSubRegions = false;
    
    typename TransformType::Pointer transform = TransformType::New();
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
      typename LabelStatsFilter::Pointer stats = LabelStatsFilter::New();
      stats->SetInput( relabel->GetOutput() );
      stats->SetLabelInput( relabel->GetOutput() );
      stats->Update();
      
      typename LabelStatsFilter::BoundingBoxType box = stats->GetBoundingBox( 1 );
      unsigned int low = box[2*wrapDim];
      unsigned int high = box[2*wrapDim + 1];
      
      unsigned int objSize = high-low+1;
      unsigned int totalmargin = origDim - objSize;
      unsigned int offset = totalmargin - totalmargin/2;
      
      if (objSize > origDim)
        {
        std::cout << "ERROR: Object is bigger than original image space, unable to detect wrapping transform" << std::endl;
        return EXIT_FAILURE;
        }
      else
        {      
        int shift = low - origDim/2;
        
        typename TransformType::OutputVectorType offsetVector;
        offsetVector.Fill( 0 );
        offsetVector[wrapDim] = (shift - offset)*reader->GetOutput()->GetSpacing()[wrapDim];
        typename LabelImageType::DirectionType direction = reader->GetOutput()->GetDirection();
        offsetVector = direction * offsetVector;
        transform->SetOffset( offsetVector );
        }
      }
    
    typename TransformWriterType::Pointer transWriter = TransformWriterType::New();
    transWriter->SetFileName( outputOption->GetValue( 0 ) );
    transWriter->SetInput( transform );
    transWriter->Update();
    
    return EXIT_SUCCESS;
    
    }
  else if( strcmp( value.c_str(), "unwrap" ) == 0 ) 
    {   
    return EXIT_SUCCESS;      
    }

  return EXIT_SUCCESS;
}

int petioleImages( itk::ants::CommandLineParser *parser )
{
  // Get dimensionality
  unsigned int dimension = 3;

  itk::ants::CommandLineParser::OptionType::Pointer dimOption =
    parser->GetOption( "dimension" );
  if( dimOption && dimOption->GetNumberOfValues() > 0 )
    {
    dimension = dimOption->Convert<unsigned int>( dimOption->GetValue() );
    }

  // Get output option
  itk::ants::CommandLineParser::OptionType::Pointer outputOption =
    parser->GetOption( "output" );

  if( !outputOption || outputOption->GetNumberOfValues() == 0 )
    {
    std::cerr << "Warning:  no output option set." << std::endl;
    }

  // Simple unary ops on fibers
  itk::ants::CommandLineParser::OptionType::Pointer imageTransformOption = parser->GetOption( "image-transform" );
  if( imageTransformOption && imageTransformOption->GetNumberOfValues() > 0 )
    {
    switch( dimension )
      {
      case 2:
        {
        ImageTransformOps<2, float>( imageTransformOption, outputOption );
        break;
        }
      case 3:
        {
        ImageTransformOps<3, float>( imageTransformOption, outputOption );
        break;
        }
      default:
        {
        std::cerr << "Unsupported dimension." << std::endl;
        return EXIT_FAILURE;
        break;
        }
      }
    return EXIT_SUCCESS;
    } 

  return EXIT_SUCCESS;
}

void InitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

  {
  std::string description =
    std::string( "operations relating to images and transforms" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "image-transform" );
  option->SetUsageOption( 0, "GetWrap[ image.ext ]" );
  option->SetUsageOption( 1, "UnWrap[ image.exd, transform.ext  ]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }
  
  {
  std::string description =
    std::string( "Ouput dependent on which option is called." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "output" );
  option->SetShortName( 'o' );
  option->SetUsageOption( 0, "output" );
  option->SetDescription( description );
  parser->AddOption( option );
  }
}

int main( int argc, char *argv[] )
{
  itk::ants::CommandLineParser::Pointer parser =
    itk::ants::CommandLineParser::New();
  parser->SetCommand( argv[0] );

  std::string commandDescription =
    std::string( "Collection of common routines for processing of streamlines" );

  parser->SetCommandDescription( commandDescription );
  InitializeCommandLineOptions( parser );

  parser->Parse( argc, argv );

  // Print the entire help menu
  itk::ants::CommandLineParser::OptionType::Pointer longHelpOption =
    parser->GetOption( "help" );
  if( argc == 1 || ( longHelpOption &&
    longHelpOption->Convert<unsigned int>( longHelpOption->GetValue() ) == 1 ) )
    {
    parser->PrintMenu( std::cout, 5, false );
    exit( EXIT_FAILURE );
    }

  itk::ants::CommandLineParser::OptionType::Pointer shortHelpOption =
    parser->GetOption( 'h' );
  if( argc == 1 || ( shortHelpOption &&
    shortHelpOption->Convert<unsigned int>( shortHelpOption->GetValue() ) == 1 ) )
    {
    parser->PrintMenu( std::cout, 5, true );
    exit( EXIT_FAILURE );
    }

  // Print the long help menu for specific items
  if( longHelpOption && longHelpOption->GetNumberOfValues() > 0
    && longHelpOption->Convert<unsigned int>( longHelpOption->GetValue() ) != 0 )
    {
    itk::ants::CommandLineParser::OptionListType options =
      parser->GetOptions();
    for( unsigned int n = 0; n < longHelpOption->GetNumberOfValues(); n++ )
      {
      std::string value = longHelpOption->GetValue( n );
      itk::ants::CommandLineParser::OptionListType::const_iterator it;
      for( it = options.begin(); it != options.end(); ++it )
        {
        const char *longName = ( ( *it )->GetLongName() ).c_str();
        if( strstr( longName, value.c_str() ) == longName  )
          {
          parser->PrintMenu( std::cout, 5, false, *it );
          }
        }
      }
    exit( EXIT_FAILURE );
    }

  // Call main routine
  petioleImages( parser );

  exit( EXIT_SUCCESS );
}

