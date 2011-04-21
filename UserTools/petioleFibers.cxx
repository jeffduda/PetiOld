#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <utility>


#include "itkMesh.h"
#include "itkSize.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkVector.h"

#include "antsCommandLineParser.h"
#include "antsCommandLineOption.h"

#include "antsVtkPolyDataFileReader.h"
#include "antsVtkPolyDataFileWriter.h"

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkBSplineControlPointImageFunction.h"
#include "itkBSplineControlPointImageFilter.h"

/*
#include "itkRGBPixel.h"

#include "itkDiffusionTensor3D.h"

#include "itkLinearInterpolateImageFunction.h"
#include "itkChainCodePath.h"

#include "ReadWriteImage.h"
#include "TensorFunctions.h"
#include "itkLabeledPointSetFileReader.h"


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

template<typename VectorType>
float VectorDistance( VectorType a, VectorType b)
{
  float sumsq = 0;
  for (unsigned int i=0; i<VectorType::Length; i++)
  {
    sumsq += ( (a[i]-b[i])*(a[i]-b[i]) );
  }
  
  return vcl_sqrt(sumsq); 
}


/* Find constant arc lenght parameterization of BSpline
 * This assumes that the spline spacing is finer than the desired parameterization 
 */
template<class FunctionType>
typename itk::VectorContainer<unsigned int, float>::Pointer ParameterizeBSplineByArcLength( FunctionType * function, float stepsize )
{

  typedef typename FunctionType::InputImageType::PixelType PointType;
  typedef typename FunctionType::PointType ArcType;
  typename itk::VectorContainer<unsigned int, float>::Pointer points = itk::VectorContainer<unsigned int, float>::New();
  points->Initialize();
  
  float epsilon = 1.e-5; // lowering this by an order of mag drastically slows the algorithm on tested data
  
  float uStep = 1.0 / (function->GetSize()[0] - 1.0);

  ArcType u1;
  u1[0] = 0.0;
  ArcType u2;
  u2[0] = uStep;
  
  bool finished = false;
  
  //std::cout << u1 << " " << std::endl;
  
  unsigned long failCount = 0;
  while ( (!finished) && (failCount < 10000))
  {
    ++failCount;
    PointType u1Pt;
    u1Pt = function->EvaluateAtParametricPoint(u1);
    PointType u2Pt;
    u2Pt = function->EvaluateAtParametricPoint(u2);
    PointType endPt;
    ArcType endU;
    endU[0] = 1.0;
    endPt = function->EvaluateAtParametricPoint(endU);
    
    float distance = VectorDistance<PointType>(u1Pt,u2Pt);
    
    float endPoint = 1-function->GetSpacing()[0];
    
    while( (distance < 2.0*stepsize) && (u2[0] < endPoint) ) 
    {
      u2[0] += uStep;
      if (u2[0] > endPoint) u2[0] = endPoint;
      u2Pt = function->EvaluateAtParametricPoint(u2);
      distance = VectorDistance<PointType>(u1Pt,u2Pt);
    }
       
    ArcType uMid;
    uMid[0] = u2[0];
    PointType uMidPt = u2Pt;
    
    ArcType tempU1 = u1;
    
    unsigned int failCountLcl = 0;
    while( (fabs(distance-stepsize) > epsilon) && (failCountLcl < 1000000) )
    {
      
      //std::cout << " - " << std::flush;
      
      uMid[0] = (tempU1[0]+u2[0])/2.0;
      uMidPt = function->EvaluateAtParametricPoint(uMid);
      
      distance = VectorDistance<PointType>(u1Pt,uMidPt);
      
      //std::cout << "  (" << tempU1 << "->" << u2 << "=" << fabs(distance-stepsize) << ")" << std::endl;
      
      
      if (distance > stepsize)
      {
        u2[0] = uMid[0];
        u2Pt = uMidPt;
      }
      else if (distance < stepsize)
      {
        tempU1[0] = uMid[0];
      }
      
      if (fabs(u2[0] - u1[0]) < epsilon)
        distance = 0;
        
      ++failCountLcl;
            
    }
    
    //std::cout << uMid << " " << std::endl;
    points->InsertElement(points->Size(), uMid[0]);
    float distFromEnd = VectorDistance<PointType>(endPt,uMidPt);
    //std::cout << "error = " << fabs(distance - stepsize) << "(" << failCountLcl << ")" << std::endl;
    
    if (distFromEnd < stepsize)
    {
      finished = true;
      //std::cout << "done." << std::endl;
    }
    else {
      float lastStep = uMid[0] - u1[0];
      u1[0] = uMid[0];
      u2[0] = u1[0] + lastStep;
      if (u2[0] > endPoint) u2[0] = endPoint;
    }
    
    if ( fabs(u1[0] - u2[0]) < epsilon )
    {
      finished = true;
      //std::cout << "done." << std::endl;
    }
  }
 
  return points;
}



template <unsigned int Dimension, class PixelType>
int FiberOps( itk::ants::CommandLineParser::OptionType *option,
  itk::ants::CommandLineParser::OptionType *outputOption = NULL )
{

  if( option->GetNumberOfParameters( 0 ) < 1 )
    {
    std::cerr << "fibers:  Incorrect number of parameters." << std::endl;
    return EXIT_FAILURE; 
    }

  std::string value = option->GetValue( 0 );
  ConvertToLowerCase( value );

  typedef itk::Mesh<float,Dimension> MeshType;
  typedef itk::ants::VtkPolyDataFileReader<MeshType> MeshReaderType;
  typedef itk::ants::VtkPolyDataFileWriter<MeshType> MeshWriterType;
    
  typename MeshReaderType::Pointer meshReader = MeshReaderType::New();



  if( strcmp( value.c_str(), "parallelize" ) == 0 ) 
    {
    
    meshReader->SetFileName( option->GetParameter( 0 ) );
    meshReader->Update(); 
    
    std::cout << "output " << outputOption->GetValue(0) << std::endl;
    
    
    typename MeshReaderType::LineSetType::Pointer lines = meshReader->GetLines();
    typename MeshType::Pointer points = meshReader->GetOutput();
    
    typename MeshReaderType::LineType refLine;
    typename MeshType::PointsContainer::Pointer refPoints = MeshType::PointsContainer::New();
    refPoints->Initialize();
    
    refLine = lines->GetElement(0);
    refPoints->Reserve( lines->GetElement(0).Size() );
    for (unsigned int i=0; i<lines->GetElement(0).Size(); i++ )
      {
      refPoints->InsertElement(i, points->GetPoints()->GetElement( refLine[i] ) );
      }
    
    int rCount = 0;  
    typename MeshReaderType::LineSetType::Pointer outLines = MeshReaderType::LineSetType::New();
    outLines->Initialize();
    
    
    for (unsigned int i=0; i<lines->Size(); i++)
      {

      typename MeshReaderType::LineType line = lines->GetElement(i);
      
      float fDist = 0;
      float rDist = 0;

      
      for (unsigned int j=0; j<refLine.Size(); j++)
        {
        typename MeshType::PointType refPoint = refPoints->GetElement( refLine[j] );
        
        float cIndex = (float)j/refLine.Size();
        int fIndex = floor( cIndex * line.Size() );
        int rIndex = line.Size() - 1 - fIndex;
        
        typename MeshType::PointType fPoint = points->GetPoints()->GetElement( line[fIndex] );
        typename MeshType::PointType rPoint = points->GetPoints()->GetElement( line[rIndex] );
        
        fDist += refPoint.EuclideanDistanceTo( fPoint )/refLine.Size();
        rDist += refPoint.EuclideanDistanceTo( rPoint )/refLine.Size();
        
        }
      
      // Need to reverse points in line
      if (rDist < fDist)
        {
        typename MeshReaderType::LineType rLine;
        rLine.SetSize( line.Size() );
        
        unsigned int * ptBuffer = new unsigned int [ line.Size() ];
        for (unsigned int j=0; j<line.Size(); j++)
          {
          ptBuffer[j] = line[j];
          }
        for (unsigned int j=0; j<line.Size(); j++)
          {
          rLine[j] = line[line.Size() - 1 - j];
          //lines->GetElement(i)[j] = ptBuffer[line.Size() - 1 - j]; 
          }
        lines->InsertElement(i,rLine);
        
        delete [] ptBuffer;
        ++rCount;
        }
       
      }
    
    std::cout << "Reversed " << rCount << "/" << lines->Size() << " fibers " << std::endl;
    
    typename MeshWriterType::Pointer meshWriter = MeshWriterType::New();
    meshWriter->SetFileName( outputOption->GetValue( 0 ) );
    meshWriter->SetInput( points );
    meshWriter->SetLines( lines );
    meshWriter->Update();
    
    return EXIT_SUCCESS;
    
    }
  else if( strcmp( value.c_str(), "reverse" ) == 0 ) 
    {
   
    meshReader->SetFileName( option->GetParameter( 0 ) );
    meshReader->Update();   
    
    typename MeshReaderType::LineSetType::Pointer lines = meshReader->GetLines();
    typename MeshType::Pointer points = meshReader->GetOutput();
    typename MeshReaderType::LineSetType::Pointer outLines = MeshReaderType::LineSetType::New();
    outLines->Initialize();
    
    
    for (unsigned int i=0; i<lines->Size(); i++)
      {

      typename MeshReaderType::LineType line = lines->GetElement(i);

      typename MeshReaderType::LineType rLine;
      rLine.SetSize( line.Size() );
        
      unsigned int * ptBuffer = new unsigned int [ line.Size() ];
      for (unsigned int j=0; j<line.Size(); j++)
        {
        ptBuffer[j] = line[j];
        }
      for (unsigned int j=0; j<line.Size(); j++)
        {
        rLine[j] = line[line.Size() - 1 - j]; 
        }
      lines->InsertElement(i,rLine);
      
      delete [] ptBuffer;
      }
    
    typename MeshWriterType::Pointer meshWriter = MeshWriterType::New();
    meshWriter->SetFileName( outputOption->GetValue( 0 ) );
    meshWriter->SetInput( points );
    meshWriter->SetLines( lines );
    meshWriter->Update();
    
    return EXIT_SUCCESS;    
    
    
    }
  else if( strcmp( value.c_str(), "trim" ) == 0 ) 
    {
    
    meshReader->SetFileName( option->GetParameter( 0 ) );
    meshReader->Update(); 
    typename MeshReaderType::LineSetType::Pointer lines = meshReader->GetLines();
    typename MeshType::Pointer points = meshReader->GetOutput();
    
    typename MeshReaderType::Pointer refReader = MeshReaderType::New();
    refReader->SetFileName( option->GetParameter( 1 ) );
    refReader->Update();
         
    typename MeshType::PointsContainer::Pointer refPoints = refReader->GetOutput()->GetPoints();
    typename MeshReaderType::LineType refLine = refReader->GetLines()->GetElement(0);

    typename MeshType::PointType startPt = refPoints->GetElement( refLine[0] );
    typename MeshType::PointType endPt = refPoints->GetElement( refLine[refLine.Size()-1] );    
    
    unsigned long * starts = new unsigned long [ lines->Size() ];
    unsigned long * ends = new unsigned long [ lines->Size() ];
    unsigned long nOutPoints = 0;  
  
    for (unsigned int i=0; i<lines->Size(); i++)
      {

      typename MeshReaderType::LineType line = lines->GetElement(i);
      
      float sDist = 0;
      float eDist = 0;
      starts[i] = 0;
      ends[i] = 0;
      
      for (unsigned int j=0; j<line.Size(); j++)
        {
        
        typename MeshType::PointType refPoint = points->GetPoints()->GetElement( line[j] );
        
        float sRefDist = refPoint.EuclideanDistanceTo( startPt );
        float eRefDist = refPoint.EuclideanDistanceTo( endPt );
        
        if (j==0)
          {
          sDist = sRefDist;
          eDist = eRefDist;
          }
        else
          {
          if (sRefDist < sDist)
            {
            sDist = sRefDist;
            starts[i] = j;
            }
          if (eRefDist < eDist)
            {
            eDist = eRefDist;
            ends[i] = j;
            }
          }
        }
      
      if (ends[i] < starts[i])
        {
        unsigned long tmp = starts[i];
        starts[i] = ends[i];
        ends[i] = tmp;
        }
    
      //std::cout << starts[i] << " - " << ends[i] << std::endl;
        
      nOutPoints += (ends[i] - starts[i] + 1);
      }
    
    
    typename MeshType::Pointer outMesh = MeshType::New();
    outMesh->GetPoints()->Reserve(nOutPoints);
    
    typename MeshType::PointDataContainer::Pointer pointData = MeshType::PointDataContainer::New();
    pointData->Initialize();
    outMesh->SetPointData( pointData );

    typename MeshWriterType::LineSetType::Pointer outLines = MeshWriterType::LineSetType::New();
    outLines->Initialize();
    
    unsigned long ptIdx = 0;
    
    for (unsigned int i=0; i<lines->Size(); i++)
      {
      typename MeshWriterType::LineType line = lines->GetElement(i);
      
      typename MeshWriterType::LineType subline;
      subline.SetSize(ends[i] - starts[i] + 1);
      unsigned long subIdx = 0;
      
      for (unsigned long j=starts[i]; j<=ends[i]; j++)
        {
        outMesh->SetPoint(ptIdx, points->GetPoints()->GetElement( line[j] ));
        subline[subIdx] = ptIdx;
        ++subIdx;
        ++ptIdx;
        }
      outLines->InsertElement(outLines->Size(), subline);
      }   

    
    typename MeshWriterType::Pointer meshWriter = MeshWriterType::New();
    meshWriter->SetFileName( outputOption->GetValue( 0 ) );
    meshWriter->SetInput( outMesh );
    meshWriter->SetLines( outLines );
    meshWriter->Update();
    
    
    std::cout << "extract not yet implemented" << std::endl;
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


template <unsigned int Dimension, class PixelType>
int PointData( itk::ants::CommandLineParser::OptionType *option,
  itk::ants::CommandLineParser::OptionType *outputOption, itk::ants::CommandLineParser::OptionListType optional )
{

  if( option->GetNumberOfParameters( 0 ) < 1 )
    {
    std::cerr << "fibers:  Incorrect number of parameters." << std::endl;
    return EXIT_FAILURE; 
    }

  std::string value = option->GetValue( 0 );
  ConvertToLowerCase( value );

  typedef itk::Mesh<float,Dimension> MeshType;
  typedef itk::ants::VtkPolyDataFileReader<MeshType> MeshReaderType;
  typedef itk::ants::VtkPolyDataFileWriter<MeshType> MeshWriterType;
    
  typename MeshReaderType::Pointer meshReader = MeshReaderType::New();



  if( strcmp( value.c_str(), "arc-length" ) == 0 ) 
    {
    
    meshReader->SetFileName( option->GetParameter( 0 ) );
    meshReader->Update(); 
    
    std::string dataName = "arc-length";
    if ( option->GetNumberOfParameters( 0 ) > 1 )
    {
      dataName = option->GetParameter( 1 );
    }

    bool normalize = false;
    itk::ants::CommandLineParser::OptionListType::iterator it = optional.begin();
    while ( !(it == optional.end()) )
      {
      
      if ( (*it)->GetLongName() == "normalize" )
        {
        normalize = true;
        }
      ++it;      
      }
    
    typename MeshReaderType::LineSetType::Pointer lines = meshReader->GetLines();
    typename MeshType::Pointer points = meshReader->GetOutput();
    
    typename MeshWriterType::MultiComponentScalarSetType::Pointer arclengths 
      = MeshWriterType::MultiComponentScalarSetType::New(); 
    arclengths->Reserve( points->GetNumberOfPoints() );
    
    float maxS = 0.0;
    
    for (unsigned int i=0; i<lines->Size(); i++)
      {

      typename MeshReaderType::LineType line = lines->GetElement(i);

      float s = 0.0;
      typename MeshWriterType::MultiComponentScalarType sValue;
      sValue.SetSize(1);
      sValue[0] = s;      
      arclengths->InsertElement( line[0], sValue );
      
      for (unsigned int j=1; j<line.Size(); j++)
        {
        typename MeshType::PointType prevPoint = points->GetPoint(j-1);
        typename MeshType::PointType point = points->GetPoint(j);
        
        typename MeshWriterType::MultiComponentScalarType sValue;
        sValue.SetSize(1);
        s = s + prevPoint.EuclideanDistanceTo( point );
        
        if (s > maxS)
          {
          maxS = s;
          }
          
        sValue[0] = s;
        arclengths->InsertElement( line[j], sValue );        
        }
      }      
     
    if ( normalize )  
      {
      for (unsigned int i=0; i<arclengths->Size(); i++)
        {
        typename MeshWriterType::MultiComponentScalarType sValue = arclengths->GetElement(i);
        sValue[0] = sValue[0] / maxS;
        arclengths->InsertElement(i, sValue);
        }
      }


    // FIXME - copy existing data from reader to retain in final output
    typename MeshWriterType::MultiComponentScalarMultiSetType::Pointer dataSets 
      = MeshWriterType::MultiComponentScalarMultiSetType::New();
    dataSets->Reserve(1);
    dataSets->SetElement(0, arclengths);
    
    typename MeshWriterType::MultiComponentScalarSetNamesType::Pointer dataNames
      =  MeshWriterType::MultiComponentScalarSetNamesType::New();
    dataNames->Reserve(1);
    dataNames->SetElement(0, dataName);

    
    typename MeshWriterType::Pointer meshWriter = MeshWriterType::New();
    meshWriter->SetFileName( outputOption->GetValue( 0 ) );
    meshWriter->SetInput( points );
    meshWriter->SetLines( lines );
    meshWriter->SetMultiComponentScalarSets( dataSets );
    meshWriter->SetMultiComponentScalarSetNames( dataNames );
    meshWriter->Update();
    
    return EXIT_SUCCESS;
    
    }
  else if( strcmp( value.c_str(), "reverse" ) == 0 ) 
    {
   
    meshReader->SetFileName( option->GetParameter( 0 ) );
    meshReader->Update();   
    
    typename MeshReaderType::LineSetType::Pointer lines = meshReader->GetLines();
    typename MeshType::Pointer points = meshReader->GetOutput();
    typename MeshReaderType::LineSetType::Pointer outLines = MeshReaderType::LineSetType::New();
    outLines->Initialize();
    
    
    for (unsigned int i=0; i<lines->Size(); i++)
      {

      typename MeshReaderType::LineType line = lines->GetElement(i);

      typename MeshReaderType::LineType rLine;
      rLine.SetSize( line.Size() );
        
      unsigned int * ptBuffer = new unsigned int [ line.Size() ];
      for (unsigned int j=0; j<line.Size(); j++)
        {
        ptBuffer[j] = line[j];
        }
      for (unsigned int j=0; j<line.Size(); j++)
        {
        rLine[j] = line[line.Size() - 1 - j]; 
        }
      lines->InsertElement(i,rLine);
      
      delete [] ptBuffer;
      }
    
    typename MeshWriterType::Pointer meshWriter = MeshWriterType::New();
    meshWriter->SetFileName( outputOption->GetValue( 0 ) );
    meshWriter->SetInput( points );
    meshWriter->SetLines( lines );
    meshWriter->Update();
    
    return EXIT_SUCCESS;    
    
    
    }
  else if( strcmp( value.c_str(), "trim" ) == 0 ) 
    {
    
    meshReader->SetFileName( option->GetParameter( 0 ) );
    meshReader->Update(); 
    typename MeshReaderType::LineSetType::Pointer lines = meshReader->GetLines();
    typename MeshType::Pointer points = meshReader->GetOutput();
    
    typename MeshReaderType::Pointer refReader = MeshReaderType::New();
    refReader->SetFileName( option->GetParameter( 1 ) );
    refReader->Update();
         
    typename MeshType::PointsContainer::Pointer refPoints = refReader->GetOutput()->GetPoints();
    typename MeshReaderType::LineType refLine = refReader->GetLines()->GetElement(0);

    typename MeshType::PointType startPt = refPoints->GetElement( refLine[0] );
    typename MeshType::PointType endPt = refPoints->GetElement( refLine[refLine.Size()-1] );    
    
    unsigned long * starts = new unsigned long [ lines->Size() ];
    unsigned long * ends = new unsigned long [ lines->Size() ];
    unsigned long nOutPoints = 0;  
  
    for (unsigned int i=0; i<lines->Size(); i++)
      {

      typename MeshReaderType::LineType line = lines->GetElement(i);
      
      float sDist = 0;
      float eDist = 0;
      starts[i] = 0;
      ends[i] = 0;
      
      for (unsigned int j=0; j<line.Size(); j++)
        {
        
        typename MeshType::PointType refPoint = points->GetPoints()->GetElement( line[j] );
        
        float sRefDist = refPoint.EuclideanDistanceTo( startPt );
        float eRefDist = refPoint.EuclideanDistanceTo( endPt );
        
        if (j==0)
          {
          sDist = sRefDist;
          eDist = eRefDist;
          }
        else
          {
          if (sRefDist < sDist)
            {
            sDist = sRefDist;
            starts[i] = j;
            }
          if (eRefDist < eDist)
            {
            eDist = eRefDist;
            ends[i] = j;
            }
          }
        }
      
      if (ends[i] < starts[i])
        {
        unsigned long tmp = starts[i];
        starts[i] = ends[i];
        ends[i] = tmp;
        }
    
      //std::cout << starts[i] << " - " << ends[i] << std::endl;
        
      nOutPoints += (ends[i] - starts[i] + 1);
      }
    
    
    typename MeshType::Pointer outMesh = MeshType::New();
    outMesh->GetPoints()->Reserve(nOutPoints);
    
    typename MeshType::PointDataContainer::Pointer pointData = MeshType::PointDataContainer::New();
    pointData->Initialize();
    outMesh->SetPointData( pointData );

    typename MeshWriterType::LineSetType::Pointer outLines = MeshWriterType::LineSetType::New();
    outLines->Initialize();
    
    unsigned long ptIdx = 0;
    
    for (unsigned int i=0; i<lines->Size(); i++)
      {
      typename MeshWriterType::LineType line = lines->GetElement(i);
      
      typename MeshWriterType::LineType subline;
      subline.SetSize(ends[i] - starts[i] + 1);
      unsigned long subIdx = 0;
      
      for (unsigned long j=starts[i]; j<=ends[i]; j++)
        {
        outMesh->SetPoint(ptIdx, points->GetPoints()->GetElement( line[j] ));
        subline[subIdx] = ptIdx;
        ++subIdx;
        ++ptIdx;
        }
      outLines->InsertElement(outLines->Size(), subline);
      }   

    
    typename MeshWriterType::Pointer meshWriter = MeshWriterType::New();
    meshWriter->SetFileName( outputOption->GetValue( 0 ) );
    meshWriter->SetInput( outMesh );
    meshWriter->SetLines( outLines );
    meshWriter->Update();
    
    
    std::cout << "extract not yet implemented" << std::endl;
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





template <unsigned int Dimension, class PixelType>
int BSplineOps( itk::ants::CommandLineParser::OptionType *option,
  itk::ants::CommandLineParser::OptionType *outputOption = NULL )
{

  std::cout << "b-spline operations" << std::endl;

  if( option->GetNumberOfParameters( 0 ) < 1 )
    {
    std::cerr << "b-spline:  Incorrect number of parameters." << std::endl;
    return EXIT_FAILURE; 
    }

  std::string value = option->GetValue( 0 );
  ConvertToLowerCase( value );

  typedef itk::Mesh<float,Dimension> MeshType;
  typedef itk::ants::VtkPolyDataFileReader<MeshType> MeshReaderType;
  typedef itk::ants::VtkPolyDataFileWriter<MeshType> MeshWriterType;
    




  if( strcmp( value.c_str(), "fit-curve" ) == 0 ) 
    {
    if( option->GetNumberOfParameters( 0 ) < 1 )
      {
      std::cerr << "b-spline:  Incorrect number of parameters." << std::endl;
      return EXIT_FAILURE; 
      }
    
    typename MeshReaderType::Pointer meshReader = MeshReaderType::New();
    meshReader->SetFileName( option->GetParameter( 0 ) );
    meshReader->Update();   
    
    std::cout << "Fit BSpline Curve" << std::endl;
       
    // Hard code for now as these parameters tend to work well for this application (i.e. fiber tractography)
    unsigned int SplineOrder = 8;
    unsigned int NumberOfControlPoints = SplineOrder+2;
    unsigned int NumberOfLevels = 4;
    unsigned int nPoints = 100;
    unsigned int nIterations = 1;
    double convergence_thresh = 0.001;
   
    if ( option->GetNumberOfParameters( 0 ) > 1 )
      {
      SplineOrder = option->Convert<unsigned int>( option->GetValue(1) );
      }  

    typedef typename MeshType::PointType           PointType;   
    typedef typename itk::Vector<float,Dimension>  ParametricPixelType;
    typedef itk::Mesh<ParametricPixelType, 1>      ParametricMeshType;
    typedef typename ParametricMeshType::PointType ParametricPointType;  
    typedef itk::Image<ParametricPixelType, 1>     ParametricImageType;
    typedef itk::BSplineScatteredDataPointSetToImageFilter<ParametricMeshType, ParametricImageType > BSplineFilterType;
    typedef itk::BSplineControlPointImageFilter< ParametricImageType > ControlPointFilterType;
    typedef itk::ImageFileWriter<ParametricImageType> BSplineImageWriterType;

    typedef std::pair<ParametricPointType,ParametricPixelType> MeshPointPair;
    std::vector<MeshPointPair> parametricPoints;
    
    typename ParametricImageType::SizeType psize;
    typename ParametricImageType::PointType porigin;
    typename ParametricImageType::SpacingType pspacing;
    typename BSplineFilterType::ArrayType ncps;
    
    psize.Fill( nPoints );
    porigin.Fill( 0.0 );
    pspacing.Fill( 1.0 / (nPoints-1) );
    ncps[0] = NumberOfControlPoints;
    
    typename BSplineFilterType::Pointer smoother = BSplineFilterType::New();
    smoother->SetSize( psize );
    smoother->SetSpacing( pspacing );
    smoother->SetOrigin( porigin );
    smoother->SetSplineOrder( SplineOrder );
    smoother->SetNumberOfControlPoints(ncps);
    smoother->SetNumberOfLevels(NumberOfLevels);

    
    std::cout << "Examining " << meshReader->GetLines()->Size() << " lines" << std::endl;
    
    unsigned long count = 0;
    
    typename ParametricMeshType::Pointer parametricMesh = ParametricMeshType::New();
    typename ParametricMeshType::PointDataContainer::Pointer dataContainer = ParametricMeshType::PointDataContainer::New();
    parametricMesh->SetPointData( dataContainer );
    parametricMesh->GetPoints()->Initialize();
    //parametricMesh->GetPoints()->Reserve( parametricPoints.size() );
    parametricMesh->GetPointData()->Initialize();
    //parametricMesh->GetPointData()->Reserve( parametricPoints.size() );
    
    for (unsigned int i=0; i<meshReader->GetLines()->Size(); i++)
    {

      typename MeshReaderType::LineType line = meshReader->GetLines()->GetElement(i);
      count += line.Size();
      
      // Use all points in each line
      for (unsigned int j=0; j < line.Size(); j++)
      {
          
        typename ParametricMeshType::PointType sPoint;
        typename ParametricMeshType::PixelType sPix;
        sPoint[0] = porigin[0] + j/static_cast<float>(line.Size()-1.0);
        for (unsigned int k=0; k<Dimension; k++)
        {
          sPix[k] = meshReader->GetOutput()->GetPoints()->GetElement( line[j] )[k];
        }
        
        unsigned long meshPointId = parametricMesh->GetNumberOfPoints();
        parametricMesh->SetPoint( meshPointId, sPoint );
        parametricMesh->SetPointData( meshPointId, sPix );
              
      }
    }
    
    smoother->SetInput( parametricMesh );
    smoother->Update();
    
    if (nIterations > 1)
      std::cout << "1 - 0.0 (" << parametricMesh->GetNumberOfPoints() << " points)" << std::endl;
    
    typename ParametricImageType::Pointer img = ParametricImageType::New();
    img->SetSpacing( smoother->GetOutput()->GetSpacing() );
    img->SetOrigin( smoother->GetOutput()->GetOrigin() );
    img->SetRegions( smoother->GetOutput()->GetLargestPossibleRegion() );
    img->Allocate();
    
    typename itk::ImageRegionIteratorWithIndex<ParametricImageType> imgIt(smoother->GetOutput(),
      smoother->GetOutput()->GetLargestPossibleRegion() );
    while(!imgIt.IsAtEnd())
    {
      img->SetPixel(imgIt.GetIndex(), imgIt.Get() );
      ++imgIt;  
    }
       
    
    float lastDiff = 1e12;
    float slope = convergence_thresh+1;
    unsigned int n=1;
    while ( (n < nIterations) && (slope > convergence_thresh) )
    {
      
      // reset parmetric mesh
      parametricMesh->Initialize();
      parametricMesh->GetPoints()->Initialize();
      dataContainer = ParametricMeshType::PointDataContainer::New();
      dataContainer->Initialize();
      parametricMesh->SetPointData( dataContainer );
      
      typename ParametricImageType::IndexType idxStart;
      idxStart[0] = 0;
      typename ParametricImageType::IndexType idxEnd;
      idxEnd[0] = nPoints-1;
      
      ParametricPixelType refLineStartPix = smoother->GetOutput()->GetPixel(idxStart);
      ParametricPixelType refLineEndPix  = smoother->GetOutput()->GetPixel(idxEnd);
          
      PointType refLineStartPt, refLineEndPt;
      for (unsigned int i=0; i<Dimension; i++)
      {
        refLineStartPt[i] = refLineStartPix[i];
        refLineEndPt[i] = refLineEndPix[i];
      }
      
      typename ControlPointFilterType::Pointer splineFunctions = ControlPointFilterType::New();
      splineFunctions->SetInput( smoother->GetOutput() );
      splineFunctions->SetSize( psize );
      splineFunctions->SetSpacing( pspacing );
      splineFunctions->SetOrigin( porigin );
      splineFunctions->SetSplineOrder( SplineOrder );
          
      for (unsigned int i=0; i<meshReader->GetLines()->Size(); i++)
      {

        //std::cout << i << " " << std::flush;

        typename MeshReaderType::LineType line = meshReader->GetLines()->GetElement(i);
        unsigned int lineStart = 0;
        unsigned int lineEnd = line.Size()-1;
        
        float lineStartDist = VectorDistance<PointType>(refLineStartPt, meshReader->GetOutput()->GetPoints()->GetElement( line[lineStart] ));
        float lineEndDist = VectorDistance<PointType>(refLineStartPt, meshReader->GetOutput()->GetPoints()->GetElement( line[lineEnd] ));
          
        // Find end points on each fiber
        for (unsigned int j=0; j < line.Size(); j++)
        {
          PointType currentPoint = meshReader->GetOutput()->GetPoints()->GetElement( line[j] );
          
          float sDist = VectorDistance<PointType>(refLineStartPt, currentPoint);
          if ( sDist < lineStartDist )
          {
            lineStartDist = sDist;
            lineStart = j;
          }
          float eDist = VectorDistance<PointType>(refLineEndPt, currentPoint);
          if ( eDist < lineEndDist )
          {
            lineEndDist = eDist;
            lineEnd = j;
          }
        }
        
        lineStart = 0;
        lineEnd = line.Size()-1;
        
        // Find closest points on spline and add to parametric mesh
        for (unsigned int j=lineStart; j<=lineEnd; j++)
        {
          PointType point = meshReader->GetOutput()->GetPoints()->GetElement( line[j] );
          typename ParametricImageType::PixelType pix;
          for (unsigned int d=0; d<Dimension; d++)
            pix[d] = point[d];
            
          typename ParametricImageType::PointType u; 
          //splineFunctions->CalculateParametersClosestToDataPoint(pix, u);
          //u[0] = static_cast<float>(j-lineStart)/(lineEnd - lineStart + 1.0);
          
          typename ParametricImageType::IndexType imageIndex;
          imageIndex[0] = 0;
          u[0] = 0;
          float minDist = VectorDistance<ParametricPixelType>(img->GetPixel(imageIndex), pix);
          for (unsigned int k=1; k<nPoints; k++)
          {
            imageIndex[0] = k;
            float currentDist = VectorDistance<ParametricPixelType>(img->GetPixel(imageIndex), pix);
            if (currentDist < minDist)
            {
              minDist = currentDist;
              u[0] = static_cast<float>(k)/(nPoints-1);
            }
          }
          
          unsigned long meshPointId = parametricMesh->GetNumberOfPoints();
          parametricMesh->SetPoint(meshPointId, u);
          parametricMesh->SetPointData(meshPointId, pix);                     
        }
        
      }
      
      // Update fit spline
      smoother->SetInput( parametricMesh );
      smoother->Update();
      
      // Calculate difference from last iteration
      float diff = 0;
      for (unsigned int i=0; i<nPoints; i++)
      {
        typename ParametricImageType::IndexType idx;
        idx[0] = i;
        
        diff += VectorDistance<ParametricPixelType>(img->GetPixel(idx), smoother->GetOutput()->GetPixel(idx));
        img->SetPixel(idx, smoother->GetOutput()->GetPixel(idx));      
        
      }
      
      slope = lastDiff - diff;
      lastDiff = diff;
      
      std::cout << (n+1) << " - " << diff/nPoints << " (slope=" << slope << ")" << std::endl;
      
      ++n;
    }

    typename BSplineImageWriterType::Pointer imgWriter = BSplineImageWriterType::New();
    imgWriter->SetFileName( outputOption->GetValue( 0 ) );
    imgWriter->SetInput( smoother->GetPhiLattice() );
    imgWriter->Update();
      
    return 1;



    }
  else if  (strcmp( value.c_str(), "resample-curve" ) == 0 ) 
    {
    
    
    // Hard code for now as these parameters tend to work well for this application (i.e. fiber tractography)
    unsigned int SplineOrder = 8;
    unsigned int NumberOfControlPoints = SplineOrder+2;
    //unsigned int NumberOfLevels = 4;
    unsigned int size = 100;
    float spacing = 1;
    
    typedef itk::Vector<PixelType, Dimension> VectorType;
    typedef itk::Image<VectorType, 1> ImageType;
    typedef itk::ImageFileReader<ImageType> ImageReaderType;
    typedef itk::BSplineControlPointImageFunction<ImageType> FunctionType;
    typedef itk::Mesh<VectorType,1>                   ParametricMeshType;
    
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName ( option->GetParameter( 0 ) );
    reader->Update();
    
    typename FunctionType::Pointer function = FunctionType::New();

    typename ImageType::SpacingType imgSpacing;
    imgSpacing.Fill( 1.0 / (size - 1.0) );
    function->SetSpacing( imgSpacing );
    typename ImageType::PointType origin;
    origin[0] = 0;
    function->SetOrigin( origin );
    typename ImageType::SizeType imgSize;
    imgSize[0] = size;
    function->SetSize( imgSize );

    typename FunctionType::ArrayType ncps;
    ncps[0] = NumberOfControlPoints;
    
    function->SetSplineOrder( SplineOrder ); 
    function->SetInputImage( reader->GetOutput() );
    
    //std::cout << "Arc length parameterize" << std::endl;
    typedef typename itk::VectorContainer<unsigned int, float> PointsHolderType;
    typename PointsHolderType::Pointer arcPoints = ParameterizeBSplineByArcLength<FunctionType>(function, spacing);
    
    typename MeshWriterType::LineType line;
    line.SetSize(arcPoints->Size());
    
    typename MeshType::Pointer outMesh = MeshType::New();
    outMesh->GetPoints()->Initialize();
    
    for (unsigned int i=0; i<arcPoints->Size(); i++)
    {
      line[i]= i;    
      
      typename MeshType::PointType outputPoint;
      typename ParametricMeshType::PointType sPoint;
      typename ParametricMeshType::PixelType sPix;
      sPoint[0] = arcPoints->GetElement(i); 
      
      sPix = function->EvaluateAtParametricPoint(sPoint);
      
      for (unsigned int j=0; j<Dimension; j++)
      {
        outputPoint[j] = sPix[j];
      }

      outMesh->SetPoint(i, outputPoint);
    }
    
    typename MeshWriterType::LineSetType::Pointer lines = MeshWriterType::LineSetType::New();
    lines->Initialize();
    lines->InsertElement(0,line);
    
    typename MeshWriterType::Pointer meshWriter = MeshWriterType::New();
    meshWriter->SetFileName( outputOption->GetValue( 0 ) );
    meshWriter->SetInput( outMesh );
    meshWriter->SetLines( lines );
    meshWriter->Update();
    
    return EXIT_SUCCESS;
    
    
    }
  else
    { 
    std::cerr << "b-spline:  Unrecognized option." << std::endl;
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
  itk::ants::CommandLineParser::OptionType::Pointer fibersOption = parser->GetOption( "fibers" );
  if( fibersOption && fibersOption->GetNumberOfValues() > 0 )
    {
    switch( dimension )
      {
      case 2:
        {
        FiberOps<2, float>( fibersOption, outputOption );
        break;
        }
      case 3:
        {
        FiberOps<3, float>( fibersOption, outputOption );
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
    
  // routines for manipulating the point-data of a fiber bundle
  itk::ants::CommandLineParser::OptionType::Pointer pointdataOption = parser->GetOption( "point-data" );
  if( pointdataOption && pointdataOption->GetNumberOfValues() > 0 )
    {
    
    
    switch( dimension )
      {
      case 2:
        {
        PointData<2, float>( pointdataOption, outputOption, parser->GetUnknownOptions() );
        break;
        }
      case 3:
        {
        PointData<3, float>( pointdataOption, outputOption, parser->GetUnknownOptions() );
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
    


  // BSpline related operations
  itk::ants::CommandLineParser::OptionType::Pointer bsplineOption = parser->GetOption( "b-spline" );
  if( bsplineOption && bsplineOption->GetNumberOfValues() > 0 )
    {
    switch( dimension )
      {
      case 2:
        {
        BSplineOps<2, float>( bsplineOption, outputOption );
        break;
        }
      case 3:
        {
        BSplineOps<3, float>( bsplineOption, outputOption );
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
    std::string( "fibers operations on fibers (streamlines)" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "fibers" );
  option->SetUsageOption( 0, "parallelize[ fibers ]" );
  option->SetUsageOption( 1, "extract[ fibers, start_index, end_index ]" );
  option->SetUsageOption( 2, "reverse[ fibers ]" );
  option->SetUsageOption( 3, "trim[ fibers, reference ]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Routines related to point-data" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "point-data" );
  option->SetUsageOption( 0, "arc-length[ fibers, data-name ]" );
  option->SetUsageOption( 1, "normalize[ fibers, data-name ]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }
  
  {
  std::string description =
    std::string( "BSpline related options" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "b-spline" );
  option->SetUsageOption( 0, "fit-curve[ fibers ]" );
  option->SetUsageOption( 1, "resample-curve[ fibers ]" );
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
  petioleFibers( parser );

  exit( EXIT_SUCCESS );
}

