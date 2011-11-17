#include <limits.h>
#include <iomanip>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkDiffusionTensor3D.h"
#include "itkImageToGraphFilter.h"
#include "itkZaleskyDiffusionTensorConnectivityGraphFilter.h"
#include "antsDijkstrasPathGraphFilter.h"
#include "antsDijkstrasGraphTraits.h"
#include "itkMaskImageToGraphFunctor.h"
//#include "itkMaskImageEdgesToGraphFunctor.h"
#include "itkGraph.h"
#include "itkBoykovGraphTraits.h"
#include "antsDijkstrasGraphTraits.h"
#include "itkGraphFileWriter.h"
#include "itkGraphFileReader.h"
#include "antsCommandLineOption.h"
#include "antsCommandLineParser.h"

#include <limits.h>

/* 
 * Forward definitions of functions
 */

// Handle graphs derived from diffusion tensor images
template <unsigned int ImageDimension, class PixelType>
int DiffusionTensorConnectivity( itk::ants::CommandLineParser::OptionType *option,
                                 itk::ants::CommandLineParser::OptionType *outputOption = NULL );

// General handler for graph generation
template <unsigned int ImageDimension, class PixelType>
int MakeGraph( itk::ants::CommandLineParser::OptionType *option,
                                 itk::ants::CommandLineParser::OptionType *outputOption = NULL );

// Handle use of Dijkstras to find "shortest-paths" in graphs
template <unsigned int ImageDimension, class PixelType>
int Dijkstras( itk::ants::CommandLineParser::OptionType *option,
               itk::ants::CommandLineParser::OptionType *outputOption = NULL );


void ConvertToLowerCase( std::string& str )
{
  std::transform( str.begin(), str.end(), str.begin(), tolower );
// You may need to cast the above line to (int(*)(int))
// tolower - this works as is on VC 7.1 but may not work on
// other compilers
} 


template <unsigned int ImageDimension, class PixelType>
int MakeGraph( itk::ants::CommandLineParser::OptionType *option,
  itk::ants::CommandLineParser::OptionType *outputOption = NULL )
{

  typedef itk::Image<unsigned int, ImageDimension> LabelImageType;
  typedef itk::ImageFileReader<LabelImageType> LabelImageReaderType;
  
  typedef itk::ants::DijkstrasGraphTraits<float, ImageDimension> GraphTraitsType;
  typedef itk::Graph<GraphTraitsType>                      GraphType;
  

  typedef itk::GraphFileWriter<GraphType>      GraphWriterType;
  typedef itk::GraphFileReader<GraphType>            GraphReaderType;

  if( option->GetNumberOfParameters( 0 ) < 1 )
    {
    std::cerr << "make-graph:  Incorrect number of parameters." << std::endl;
    return EXIT_FAILURE; 
    }



  std::string value = option->GetValue( 0 );
  ConvertToLowerCase( value );
  
  typename GraphType::Pointer output = NULL;
  bool handleOutput = true;

  if( strcmp( value.c_str(), "binary-image" ) == 0 ) 
    {
    typename LabelImageReaderType::Pointer reader = LabelImageReaderType::New();
    reader->SetFileName( option->GetParameter( 0 ) );
    reader->Update(); 
    
    typedef itk::MaskImageToGraphFunctor<LabelImageType,GraphType> FunctorType;
    typename FunctorType::Pointer functor = FunctorType::New();
    functor->SetExcludeBackground( true );
    functor->ActivateAllNeighbors();
    
    typedef itk::ImageToGraphFilter<LabelImageType,GraphType> GraphSourceType;
    typename GraphSourceType::Pointer makeGraph = GraphSourceType::New();
    makeGraph->SetImageToGraphFunctor( functor );
    makeGraph->SetInput( reader->GetOutput() );
    makeGraph->Update();
    
    output = makeGraph->GetOutput(); 
    
    std::cout << "# Nodes = " << makeGraph->GetOutput()->GetTotalNumberOfNodes() << std::endl;
    std::cout << "# Edges = " << makeGraph->GetOutput()->GetTotalNumberOfEdges() << std::endl;

    }
  else if (strcmp( value.c_str(), "diffusion-tensor-zalesky" ) == 0 )
    {
    DiffusionTensorConnectivity<ImageDimension, PixelType>( option, outputOption );
    handleOutput = false;
    }
  else
    { 
    std::cerr << "--make-graph:  Unrecognized option." << std::endl;
    return EXIT_FAILURE;
    }

  if( outputOption && handleOutput )
    {
    typename GraphWriterType::Pointer writeGraph = GraphWriterType::New();
    writeGraph->SetFileName( outputOption->GetValue( 0 ) ); 
    writeGraph->SetInput( output );
    writeGraph->Update();
    }




  return EXIT_SUCCESS;
}

template <unsigned int ImageDimension, class PixelType>
int Dijkstras( itk::ants::CommandLineParser::OptionType *option,
  itk::ants::CommandLineParser::OptionType *outputOption = NULL )
{
  typedef itk::Image<unsigned int, ImageDimension> LabeledImageType;
  typedef itk::ImageFileReader<LabeledImageType> LabeledImageReader;

  typedef itk::ants::DijkstrasGraphTraits<float, 3> GraphTraitsType;
  typedef itk::Graph<GraphTraitsType>               GraphType;
  typedef GraphType::NodePointerType                NodePointerType;
  typedef itk::GraphFileWriter<GraphType>           GraphWriterType;
  typedef itk::GraphFileReader<GraphType>           GraphReaderType;
  typedef itk::ants::DijkstrasPathGraphFilter<GraphType> FilterType;

  typedef std::list< GraphType::NodeIdentifierType >    ListType;

  typedef itk::Mesh<float, 3> MeshType;

  typedef itk::ants::VtkPolyDataFileWriter<MeshType> MeshWriterType;

  typename GraphReaderType::Pointer reader = GraphReaderType::New();
  reader->SetFileName( option->GetParameter( 0 ) );
  reader->Update();   
  typename GraphType::Pointer graph = reader->GetOutput();

  std::string value = option->GetValue( 0 );
  ConvertToLowerCase( value );

  typename GraphType::Pointer output = NULL;

  if( strcmp( value.c_str(), "dijkstras-shortest-paths" ) == 0 ) 
    {
    if( option->GetNumberOfParameters( 0 ) < 2 )
      {
      std::cerr << "Incorrect number of parameters" << std::endl;
      return EXIT_FAILURE;
      }

    typename LabeledImageReader::Pointer reader2 = LabeledImageReader::New();
    reader2->SetFileName( option->GetParameter(1) );
    reader2->Update();
    
    for (unsigned int i=0; i<graph->GetTotalNumberOfEdges(); i++)
      {
      // FIXME - Replace with using edge-weights via data-name
      // passed as option

      typename LabeledImageType::PointType ptS;
      typename LabeledImageType::PointType ptT;
      
      reader2->GetOutput()->TransformIndexToPhysicalPoint( graph->GetNodePointer(graph->GetEdgePointer(i)->SourceIdentifier)->ImageIndex, ptS);
      reader2->GetOutput()->TransformIndexToPhysicalPoint( graph->GetNodePointer(graph->GetEdgePointer(i)->TargetIdentifier)->ImageIndex, ptT);
      graph->GetEdgePointer(i)->Weight = ptS.EuclideanDistanceTo( ptT );
      
      typename LabeledImageType::PointType::VectorType vector = ptS - ptT;
      graph->GetEdgePointer(i)->Weight = vector.GetNorm();
      }

    typename FilterType::Pointer dijkstras = FilterType::New();
    dijkstras->SetInput( graph );
    dijkstras->SetLabels( reader2->GetOutput() );
    //dijkstras->SetSourceNodes( sources );
    //dijkstras->SetTargetNodes( targets );
    dijkstras->Update();

    graph = dijkstras->GetOutput();

    typename FilterType::NodeListContainerType path = dijkstras->GetTree( );
    
    // FIXME - Extract points that are in the paths for smaller output
    // file
    ListType pointList;
    for (unsigned int i=0; i<path.size(); i++)
      {
      for (unsigned int j=0; j<path[j].size(); j++)
        {
        //GraphType::NodeIdentifierType nodeId = path[j][i]->Identifier;
        //pointList.insert( nodeId );
        pointList.push_back( path[j][i]->Identifier );
        }
      }
    pointList.unique();

    std::cout << "# Paths = " << path.size() << std::endl;
    std::cout << "# Points = " << pointList.size() << std::endl;

    MeshType::Pointer mesh = MeshType::New();
    mesh->GetPoints()->Initialize();
    mesh->GetPoints()->Reserve( pointList.size() );

    ListType::iterator it = pointList.begin();
    while (it != pointList.end() )
      {
      

    /*
    vtkPoints * vtkPoints = vtkPoints::New();
    vtkPoints->Initialize();
    
    vtkCellArray * vtkPath = vtkCellArray::New();
    vtkPath->Initialize();
    
    vtkIdType nPoints = 0;
    
    for (unsigned int j=0; j<path.size(); j++)
      {
      vtkIdType * pathLine = new vtkIdType [ path[j].size() ];
      
      for (unsigned int i=0; i<path[j].size(); i++)
        {
        vtkPoints->InsertNextPoint(path[j][i]->ImageIndex[0], path[j][i]->ImageIndex[1], path[j][i]->ImageIndex[2] );
        pathLine[i] = nPoints;
        ++nPoints;
        }
      vtkPath->InsertNextCell( path[j].size(), pathLine );
    
      delete [] pathLine;
      }
      
    vtkPolyData * polydata = vtkPolyData::New();
    polydata->SetPoints( vtkPoints );
    polydata->SetLines( vtkPath );
    
    vtkPolyDataWriter * writer = vtkPolyDataWriter::New();
    writer->SetInput( polydata );
    writer->SetFileName( outputOption->GetValue( 0 ).c_str() );
    writer->Update();
    */
      }
   }
  else
    { 
    std::cerr << "Dijkstras:  Unrecognized option." << std::endl;
    return EXIT_FAILURE;
    }

  if( outputOption )
    {
    typename GraphWriterType::Pointer writeGraph = GraphWriterType::New();
    writeGraph->SetFileName( outputOption->GetValue( 0 ) ); 
    writeGraph->SetInput( graph );
    //writeGraph->Update();
    
    }

  return EXIT_SUCCESS;
}

template <unsigned int ImageDimension, class PixelType>
int GraphWeights( itk::ants::CommandLineParser::OptionType *option,
  itk::ants::CommandLineParser::OptionType *outputOption = NULL )
{
  typedef itk::Image<float, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ImageReader;

  typedef itk::ants::DijkstrasGraphTraits<float, 3> GraphTraitsType;
  typedef itk::Graph<GraphTraitsType>           GraphType;
  typedef GraphType::NodePointerType                   NodePointerType;
  typedef itk::GraphFileWriter<GraphType>      GraphWriterType;
  typedef itk::GraphFileReader<GraphType>      GraphReaderType;
  typedef itk::ants::DijkstrasPathGraphFilter<GraphType> FilterType;
  
  typedef itk::DiffusionTensor3D<float> TensorType;
  typedef itk::Image<TensorType,ImageDimension> TensorImageType;
  typedef itk::ImageFileReader<TensorImageType> TensorImageReader;

  typename GraphReaderType::Pointer reader = GraphReaderType::New();
  reader->SetFileName( option->GetParameter( 0 ) );
  reader->Update();   
  typename GraphType::Pointer graph = reader->GetOutput();
  
  std::string value = option->GetValue( 0 );
  ConvertToLowerCase( value );

  typename GraphType::Pointer output = NULL;

  if( strcmp( value.c_str(), "edge-euclidean-distance" ) == 0 ) 
    {
    if( option->GetNumberOfParameters( 0 ) < 2 )
      {
      std::cerr << "Incorrect number of parameters" << std::endl;
      return EXIT_FAILURE;
      }

    typename ImageReader::Pointer imgReader = ImageReader::New();
    imgReader->SetFileName( option->GetParameter(1) );
    imgReader->Update();  
      
    for (unsigned int i=0; i<graph->GetTotalNumberOfEdges(); i++)
      {
      typename ImageType::PointType ptS;
      typename ImageType::PointType ptT;
      
      imgReader->GetOutput()->TransformIndexToPhysicalPoint( graph->GetNodePointer(graph->GetEdgePointer(i)->SourceIdentifier)->ImageIndex, ptS);
      imgReader->GetOutput()->TransformIndexToPhysicalPoint( graph->GetNodePointer(graph->GetEdgePointer(i)->TargetIdentifier)->ImageIndex, ptT);
      graph->GetEdgePointer(i)->Weight = ptS.EuclideanDistanceTo( ptT );
      }         
   }
  else if( strcmp( value.c_str(), "edge-tensor-fa" ) == 0 ) 
    {
    if( option->GetNumberOfParameters( 0 ) < 2 )
      {
      std::cerr << "Incorrect number of parameters" << std::endl;
      return EXIT_FAILURE;
      }

    typename TensorImageReader::Pointer imgReader = TensorImageReader::New();
    imgReader->SetFileName( option->GetParameter(1) );
    imgReader->Update(); 

    for (unsigned int i=0; i<graph->GetTotalNumberOfEdges(); i++)
      {
      TensorType tenS = imgReader->GetOutput()->GetPixel( graph->GetNodePointer(graph->GetEdgePointer(i)->SourceIdentifier)->ImageIndex);
      TensorType tenT = imgReader->GetOutput()->GetPixel( graph->GetNodePointer(graph->GetEdgePointer(i)->TargetIdentifier)->ImageIndex);
      graph->GetEdgePointer(i)->Weight = 1.00001 - tenS.GetFractionalAnisotropy() * tenT.GetFractionalAnisotropy();
      }         
    }   
  else if( strcmp( value.c_str(), "edge-tensor-pdd" ) == 0 ) 
    {
    if( option->GetNumberOfParameters( 0 ) < 2 )
      {
      std::cerr << "Incorrect number of parameters" << std::endl;
      return EXIT_FAILURE;
      }

    typename TensorImageReader::Pointer imgReader = TensorImageReader::New();
    imgReader->SetFileName( option->GetParameter(1) );
    imgReader->Update(); 

    for (unsigned int i=0; i<graph->GetTotalNumberOfEdges(); i++)
      {
      TensorType tenS = imgReader->GetOutput()->GetPixel( graph->GetNodePointer(graph->GetEdgePointer(i)->SourceIdentifier)->ImageIndex);
      TensorType tenT = imgReader->GetOutput()->GetPixel( graph->GetNodePointer(graph->GetEdgePointer(i)->TargetIdentifier)->ImageIndex);
      
      typename TensorType::EigenValuesArrayType sourceEVals;
      typename TensorType::EigenValuesArrayType targetEVals;
      
      typename TensorType::EigenVectorsMatrixType sourceEVecs;
      typename TensorType::EigenVectorsMatrixType targetEVecs;
      
      tenS.ComputeEigenAnalysis(sourceEVals,sourceEVecs);
      tenT.ComputeEigenAnalysis(targetEVals,targetEVecs);
        
      typename ImageType::PointType::VectorType sourcePDD;
      typename ImageType::PointType::VectorType targetPDD;

      for (unsigned int j=0; j<ImageDimension; j++)
        {
        sourcePDD[j] = sourceEVecs[2][j];
        targetPDD[j] = targetEVecs[2][j];
        }
      sourcePDD.Normalize();
      targetPDD.Normalize();
        
      typename TensorImageType::PointType ptS;
      typename TensorImageType::PointType ptT;
      
      imgReader->GetOutput()->TransformIndexToPhysicalPoint( graph->GetNodePointer(graph->GetEdgePointer(i)->SourceIdentifier)->ImageIndex, ptS);
      imgReader->GetOutput()->TransformIndexToPhysicalPoint( graph->GetNodePointer(graph->GetEdgePointer(i)->TargetIdentifier)->ImageIndex, ptT);
      graph->GetEdgePointer(i)->Weight = ptS.EuclideanDistanceTo( ptT );
      
      typename TensorImageType::PointType::VectorType vector = ptS - ptT; 
      //float distance = vector.GetNorm();       
      vector.Normalize();
      
      graph->GetEdgePointer(i)->Weight = 1.0 - vcl_fabs(vector*sourcePDD) * vcl_fabs(vector*targetPDD);
      }         
    } 
  else if( strcmp( value.c_str(), "edge-tensor-tsp" ) == 0 ) 
    {
    if( option->GetNumberOfParameters( 0 ) < 2 )
      {
      std::cerr << "Incorrect number of parameters" << std::endl;
      return EXIT_FAILURE;
      }

    typename TensorImageReader::Pointer imgReader = TensorImageReader::New();
    imgReader->SetFileName( option->GetParameter(1) );
    imgReader->Update(); 

    for (unsigned int i=0; i<graph->GetTotalNumberOfEdges(); i++)
      {
      TensorType tenS = imgReader->GetOutput()->GetPixel( graph->GetNodePointer(graph->GetEdgePointer(i)->SourceIdentifier)->ImageIndex);
      TensorType tenT = imgReader->GetOutput()->GetPixel( graph->GetNodePointer(graph->GetEdgePointer(i)->TargetIdentifier)->ImageIndex);
      
      float tsp = tenS[0]*tenT[0] + tenS[3]*tenT[3] + tenS[5]*tenT[5] + 2.0*( tenS[1]*tenT[1] + tenS[2]*tenT[2] + tenS[4]*tenT[4]);
      
      graph->GetEdgePointer(i)->Weight = tsp;
      }         
    }     
  else if( strcmp( value.c_str(), "edge-tensor-projection" ) == 0 ) 
    {
    if( option->GetNumberOfParameters( 0 ) < 2 )
      {
      std::cerr << "Incorrect number of parameters" << std::endl;
      return EXIT_FAILURE;
      }

    typename TensorImageReader::Pointer imgReader = TensorImageReader::New();
    imgReader->SetFileName( option->GetParameter(1) );
    imgReader->Update(); 

    for (unsigned int i=0; i<graph->GetTotalNumberOfEdges(); i++)
      {
      TensorType tenS = imgReader->GetOutput()->GetPixel( graph->GetNodePointer(graph->GetEdgePointer(i)->SourceIdentifier)->ImageIndex);
      TensorType tenT = imgReader->GetOutput()->GetPixel( graph->GetNodePointer(graph->GetEdgePointer(i)->TargetIdentifier)->ImageIndex);
      
      typename TensorType::MatrixType matS;
      typename TensorType::MatrixType matT;
      
      matS(0,0) = tenS[0];
      matS(0,1) = matS(1,0) = tenS[1];
      matS(0,2) = matS(2,0) = tenS[2];
      matS(1,1) = tenS[3];
      matS(1,2) = matS(2,1) = tenS[4];
      matS(2,2) = tenS[5];
      
      matT(0,0) = tenT[0];
      matT(0,1) = matT(1,0) = tenT[1];
      matT(0,2) = matT(2,0) = tenT[2];
      matT(1,1) = tenT[3];
      matT(1,2) = matT(2,1) = tenT[4];
      matT(2,2) = tenT[5];
      
      typename TensorType::EigenValuesArrayType sourceEVals;
      typename TensorType::EigenValuesArrayType targetEVals;
      
      typename TensorType::EigenVectorsMatrixType sourceEVecs;
      typename TensorType::EigenVectorsMatrixType targetEVecs;
      
      tenS.ComputeEigenAnalysis(sourceEVals,sourceEVecs);
      tenT.ComputeEigenAnalysis(targetEVals,targetEVecs);
        
      typename ImageType::PointType::VectorType sourcePDD;
      typename ImageType::PointType::VectorType targetPDD;

      for (unsigned int j=0; j<ImageDimension; j++)
        {
        sourcePDD[j] = sourceEVecs[2][j];
        targetPDD[j] = targetEVecs[2][j];
        }
      sourcePDD.Normalize();
      targetPDD.Normalize();
        
      typename TensorImageType::PointType ptS;
      typename TensorImageType::PointType ptT;
      
      imgReader->GetOutput()->TransformIndexToPhysicalPoint( graph->GetNodePointer(graph->GetEdgePointer(i)->SourceIdentifier)->ImageIndex, ptS);
      imgReader->GetOutput()->TransformIndexToPhysicalPoint( graph->GetNodePointer(graph->GetEdgePointer(i)->TargetIdentifier)->ImageIndex, ptT);
      graph->GetEdgePointer(i)->Weight = ptS.EuclideanDistanceTo( ptT );
      
      typename TensorImageType::PointType::VectorType vector = ptS - ptT; 
      //float distance = vector.GetNorm();       
      vector.Normalize();
      
      graph->GetEdgePointer(i)->Weight = 1.0 - vcl_fabs(vector*sourcePDD) * vcl_fabs(vector*targetPDD);
      }         
    }     
  else
    { 
    std::cerr << "GraphWeight:  Unrecognized option." << std::endl;
    return EXIT_FAILURE;
    }

  if( outputOption )
    {
    typename GraphWriterType::Pointer writeGraph = GraphWriterType::New();
    writeGraph->SetFileName( outputOption->GetValue( 0 ) ); 
    writeGraph->SetInput( graph );
    writeGraph->Update();
    }

  return EXIT_SUCCESS;
}



template <unsigned int ImageDimension, class PixelType>
int DiffusionTensorConnectivity( itk::ants::CommandLineParser::OptionType *option,
  itk::ants::CommandLineParser::OptionType *outputOption = NULL )
{
	
  typedef itk::DiffusionTensor3D<float> DiffusionTensor;
  typedef itk::Image<DiffusionTensor, ImageDimension> ImageType;
  typedef itk::Image<unsigned int, ImageDimension> LabelImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
  
  typedef itk::ants::DijkstrasGraphTraits<unsigned int, 3> GraphTraitsType;
  typedef itk::Graph<GraphTraitsType>           GraphType;
  typedef itk::ZaleskyDiffusionTensorConnectivityGraphFilter<LabelImageType,ImageType,GraphType> GraphSourceType;
  typedef itk::GraphFileWriter<GraphType>      GraphWriterType;
  typedef itk::GraphFileReader<GraphType>      GraphReaderType;

  if( option->GetNumberOfParameters( 0 ) < 1 )
    {
    std::cerr << "DiffusionTensorConnectivity:  Incorrect number of parameters." << std::endl;
    return EXIT_FAILURE; 
    }

  typename ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( option->GetParameter( 0 ) );
  reader1->Update(); 

  typename LabelImageType::Pointer labels = NULL;
  if( option->GetNumberOfParameters( 0 ) > 1) 
    {
    typename LabelReaderType::Pointer reader2 = LabelReaderType::New();
    reader2->SetFileName( option->GetParameter(1) );
    reader2->Update();
    labels = reader2->GetOutput();
    }

  std::string value = option->GetValue( 0 );
  ConvertToLowerCase( value );
  
  typename GraphType::Pointer output = NULL;

  if( strcmp( value.c_str(), "diffusion-tensor-zalesky" ) == 0 ) 
    {

    std::cout << "Computing the DTI-based Zalesky graph" << std::endl;

    typename GraphSourceType::Pointer makeGraph = GraphSourceType::New();
    makeGraph->SetInput( labels ); 
    makeGraph->SetDiffusionTensorImage( reader1->GetOutput() );
    makeGraph->GetImageToGraphFunctor()->SetExcludeBackground( true );
    makeGraph->Update();
    
    output = makeGraph->GetOutput(); 
    
    std::cout << "# Nodes = " << makeGraph->GetOutput()->GetTotalNumberOfNodes() << std::endl;
    std::cout << "# Edges = " << makeGraph->GetOutput()->GetTotalNumberOfEdges() << std::endl;

    }
  else
    { 
    std::cerr << "DiffusionTensorConnectivity:  Unrecognized option." << std::endl;
    return EXIT_FAILURE;
    }

  if( outputOption )
    {
    typename GraphWriterType::Pointer writeGraph = GraphWriterType::New();
    writeGraph->SetFileName( outputOption->GetValue( 0 ) ); 
    writeGraph->SetInput( output );
    writeGraph->Update();
    }




  return EXIT_SUCCESS;
}



template <unsigned int ImageDimension, class PixelType>
int TimeConnectivity( itk::ants::CommandLineParser::OptionType *option,
  itk::ants::CommandLineParser::OptionType *outputOption = NULL )
{
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typename ImageType::Pointer outputImage = NULL;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( option->GetParameter( 0 ) );
  reader->Update();

  //PixelType constant = option->Convert<PixelType>( option->GetParameter( 1 ) );

  std::string value = option->GetValue( 0 );
  ConvertToLowerCase( value );

  if( strcmp( value.c_str(), "pearson-correlation" ) == 0 )
    {
		std::cout << "Pearson Correlation not yet implemented" << std::endl;
    }
  else if( strcmp( value.c_str(), "spearman-correlation" ) == 0 )
    {
		std::cout << "Spearman Correlation not yet implemented" << std::endl;
    }
  else
    {
    std::cerr << "TimeConnectivity:  Unrecognized option." << std::endl;
    return EXIT_FAILURE;
    }

  if( outputOption )
    {
		std::cout << "Writer not yet implemented" << std::endl;
    }

  return EXIT_SUCCESS;
}

int Connectivity( itk::ants::CommandLineParser *parser )
{
  // Get dimensionality
  unsigned int dimension = 3;

  itk::ants::CommandLineParser::OptionType::Pointer dimOption =
    parser->GetOption( "image-dimensionality" );
  if( dimOption && dimOption->GetNumberOfValues() > 0 )
    {
    dimension = dimOption->Convert<unsigned int>( dimOption->GetValue() );
    }

  // Get output option
  itk::ants::CommandLineParser::OptionType::Pointer outputOption =
    parser->GetOption( "output" );
  itk::ants::CommandLineParser::OptionType::Pointer imageInfoOption =
    parser->GetOption( "image-info" );
  if( !imageInfoOption && (!outputOption || outputOption->GetNumberOfValues() == 0 ))
    {
    std::cerr << "Warning:  no output option set." << std::endl;
    }

  // Create a graph
  itk::ants::CommandLineParser::OptionType::Pointer makeGraphOption = parser->GetOption( "make-graph" );
  if( makeGraphOption && makeGraphOption->GetNumberOfValues() > 0 )
    {
    switch( dimension )
      {
      case 2:
        {
        MakeGraph<2, float>( makeGraphOption, outputOption );
        break;
        }
      case 3:
        {
        MakeGraph<3, float>( makeGraphOption, outputOption );
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

  // Voxel-wise binary operations
  itk::ants::CommandLineParser::OptionType::Pointer timeOption =
    parser->GetOption( "time" );
  if( timeOption && timeOption->GetNumberOfValues() > 0 )
    {
    switch( dimension )
      {
      case 2:
        {
        //TimeConnectivity<2, float>( timeOption, outputOption );
        break;
        }
      case 3:
        {
        //TimeConnectivity<3, float>( timeOption, outputOption );
        break;
        }
      case 4:
        {
        //TimeConnectivity<4, float>( timeOption, outputOption );
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
    
  itk::ants::CommandLineParser::OptionType::Pointer dijkstrasOption =
    parser->GetOption( "graph" );
  if( dijkstrasOption && dijkstrasOption->GetNumberOfValues() > 0 )
    {
    switch( dimension )
      {
      case 2:
        {
        //Dijkstras<2, float>( dijkstrasOption, outputOption );
        break;
        }
      case 3:
        {
        Dijkstras<3, float>( dijkstrasOption, outputOption );
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
    
    
    // Voxel-wise binary operations
  itk::ants::CommandLineParser::OptionType::Pointer graphWeightsOption =
    parser->GetOption( "graph-weights" );
  if( graphWeightsOption && graphWeightsOption->GetNumberOfValues() > 0 )
    {
    switch( dimension )
      {
      case 2:
        {
        //GraphWeights<2, float>( graphWeightsOption, outputOption );
        break;
        }
      case 3:
        {
        GraphWeights<3, float>( graphWeightsOption, outputOption );
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
    std::string( "Algorithms for creating a graph" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "make-graph" );
  option->SetUsageOption( 0, "binary-image[image]" );
  option->SetUsageOption( 1, "diffusion-tensor-zalesky[ image, labels ]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Connectivity as measured with time-signals such as fMRI" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "time" );
  option->SetUsageOption( 0, "pearson-correlation[image1,<labelimage>]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Graph methods" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "graph" );
  option->SetUsageOption( 0, "dijkstras-shortest-paths[graph,labelimage]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }
  
  {
  std::string description =
    std::string( "Graph edge and node weighting methods" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "graph-weights" );
  option->SetUsageOption( 0, "edge-euclidean-distance[graph,referenceimage]");
  option->SetUsageOption( 1, "edge-tensor-fa[graph,tensorimage]");
  option->SetUsageOption( 2, "edge-tensor-ppd[ graph, tensorimage ]");
  option->SetUsageOption( 3, "edge-tensor-tsp[ graph, tensorimage ]");
  option->SetUsageOption( 4, "edge-tensor-projection[ graph, tensorimage ]");
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Ouput dependent on which option is called." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "output" );
  option->SetShortName( 'o' );
  option->SetUsageOption( 0, "outputImage" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Print the help menu (short version)." );

  OptionType::Pointer option = OptionType::New();
  option->SetShortName( 'h' );
  option->SetDescription( description );
  option->AddValue( std::string( "0" ) );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Print the help menu (long version)." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "help" );
  option->SetDescription( description );
  option->AddValue( std::string( "0" ) );
  parser->AddOption( option );
  }
}

/* 
 * Main progam
 */
int main( int argc, char *argv[] )
{
  itk::ants::CommandLineParser::Pointer parser =
    itk::ants::CommandLineParser::New();
  parser->SetCommand( argv[0] );

  std::string commandDescription =
    std::string( "Collection of common routines for estimating connectivity" );

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
  Connectivity( parser );

  exit( EXIT_SUCCESS );
} 
