#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <utility>


#include "itkMesh.h"
#include "itkSize.h"
#include "itkImage.h"
#include "itkVector.h"
#include "antsCommandLineParser.h"
#include "antsCommandLineOption.h"


/*
#include "itkRGBPixel.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkDiffusionTensor3D.h"

#include "itkLinearInterpolateImageFunction.h"
#include "itkChainCodePath.h"

#include "ReadWriteImage.h"
#include "TensorFunctions.h"
#include "itkLabeledPointSetFileReader.h"
#include "antsVtkPolyDataFileReader.h"
#include "antsVtkPolyDataFileWriter.h"

#include "itkGrahamScanConvexHull.h"

#include "itkFiberFileReader.h"
#include "itkFiberFileWriter.h"
#include "itkVectorContainer.h"

#include "itkEllipseFunctions.h"
*/

void ConvertToLowerCase( std::string& str )
{
  std::transform( str.begin(), str.end(), str.begin(), tolower );
// You may need to cast the above line to (int(*)(int))
// tolower - this works as is on VC 7.1 but may not work on
// other compilers
} 


template <unsigned int ImageDimension, class PixelType>
int UnaryFiberOps( itk::ants::CommandLineParser::OptionType *option,
  itk::ants::CommandLineParser::OptionType *outputOption = NULL )
{

  if( option->GetNumberOfParameters( 0 ) < 1 )
    {
    std::cerr << "unary:  Incorrect number of parameters." << std::endl;
    return EXIT_FAILURE; 
    }

  std::string value = option->GetValue( 0 );
  ConvertToLowerCase( value );


  if( strcmp( value.c_str(), "parallelize" ) == 0 ) 
    {
    std::cout << "parallelize not yet implemented" << std::endl;
    }
  else if( strcmp( value.c_str(), "reverse" ) == 0 ) 
    {
    std::cout << "reverse not yet implemented" << std::endl;
    }
  else if( strcmp( value.c_str(), "extract" ) == 0 ) 
    {
    std::cout << "extract not yet implemented" << std::endl;
    }    
  else
    { 
    std::cerr << "unary:  Unrecognized option." << std::endl;
    return EXIT_FAILURE;
    }

  if( outputOption )
    {

    }




  return EXIT_SUCCESS;
}



int petioleFibers( itk::ants::CommandLineParser *parser )
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
  itk::ants::CommandLineParser::OptionType::Pointer imageInfoOption =
    parser->GetOption( "fiber-info" );
  if( !imageInfoOption && (!outputOption || outputOption->GetNumberOfValues() == 0 ))
    {
    std::cerr << "Warning:  no output option set." << std::endl;
    }

  // Simple unary ops on fibers
  itk::ants::CommandLineParser::OptionType::Pointer unaryFiberOption = parser->GetOption( "unary" );
  if( unaryFiberOption && unaryFiberOption->GetNumberOfValues() > 0 )
    {
    switch( dimension )
      {
      case 2:
        {
        UnaryFiberOps<2, float>( unaryFiberOption, outputOption );
        break;
        }
      case 3:
        {
        UnaryFiberOps<3, float>( unaryFiberOption, outputOption );
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
    std::string( "Utility operations on fibers (streamlines)" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "unary" );
  option->SetUsageOption( 0, "parallelize[ fibers ]" );
  option->SetUsageOption( 1, "extract[ fibers, start_index, end_index ]" );
  option->SetUsageOption( 2, "reverse[ fibers  ]" );
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
  petioleFibers( parser );

  exit( EXIT_SUCCESS );
}

