/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: antsVtkPolyDataFileReader.txx,v $
  Language:  C++
  Date:      $Date: 2009/03/13 19:48:16 $
  Version:   $Revision: 1.23 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __antsVtkPolyDataFileReader_txx
#define __antsVtkPolyDataFileReader_txx

#include "antsVtkPolyDataFileReader.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLabelContourImageFilter.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkByteSwapper.h"

#include <fstream>
#include <stdio.h>
#include <string>

namespace itk {
namespace ants {

//
// Constructor
//
template<class TOutputMesh>
VtkPolyDataFileReader<TOutputMesh>
::VtkPolyDataFileReader()
{
  this->m_RandomPercentage = 1.0;
  this->m_ExtractBoundaryPoints = false;

  this->m_MultiComponentScalarSets = NULL;
  this->m_CellMultiComponentScalarSets = NULL;
  this->m_Lines = NULL;
  //
  // Create the output
  //
  typename TOutputMesh::Pointer output = TOutputMesh::New();
  this->ProcessObject::SetNumberOfRequiredOutputs( 1 );
  this->ProcessObject::SetNthOutput( 0, output.GetPointer() );
}

template<class TOutputMesh>
void
VtkPolyDataFileReader<TOutputMesh>
::GenerateData()
{
  if( this->m_FileName == "" )
    {
    itkExceptionMacro( "No input FileName" );
    return;
    }

  //
  // Read input file
  //
  std::ifstream inputFile( m_FileName.c_str() );

  if( !inputFile.is_open() )
    {
    itkExceptionMacro("Unable to open file\n"
        "inputFilename= " << m_FileName );
    return;
    }
  else
    {
    inputFile.close();
    }

  /**
   * Get filename extension
   */

  std::string::size_type pos = this->m_FileName.rfind( "." );
  std::string extension( this->m_FileName, pos+1, this->m_FileName.length()-1 );

  if( extension == "txt" )
    {
    
    }
  else if( extension == "vtk" )
    {
    this->ReadVTKFile();
    }
  else // try reading the file as an image
    {
    
    }

  if( this->m_RandomPercentage < 1.0 )
    {
    typedef Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
    typename GeneratorType::Pointer generator = GeneratorType::New();
    generator->SetSeed();

    typename OutputMeshType::Pointer output = OutputMeshType::New();
    output->Initialize();

    if( this->GetOutput()->GetNumberOfPoints() > 0 )
      {
      typename OutputMeshType::PointsContainerIterator It =
        this->GetOutput()->GetPoints()->Begin();

      unsigned long count = 0;
      while( It != this->GetOutput()->GetPoints()->End() )
        {
        if( generator->GetVariateWithClosedRange() <= this->m_RandomPercentage )
          {
          output->SetPoint( count, It.Value() );
          PixelType label;
          bool elementExists = this->GetOutput()->GetPointData()
            ->GetElementIfIndexExists( It.Index(), &label );
          if( elementExists )
            {
            output->SetPointData( count, label );
            }
          count++;
          }
        ++It;
        }
      }
    this->GraftOutput( output );
    }
  this->m_LabelSet.clear();

  /**
   * If the number of points does not match the number of 
   * point data, fill the point data with zeros
   */
  if( !this->GetOutput()->GetPointData() || 
    this->GetOutput()->GetPointData()->Size() != 
    this->GetOutput()->GetPoints()->Size() )
    {
    //itkWarningMacro( "Number of points does not match number of labels. "
    //  << "Filling point data with label zero." ); 
     typename OutputMeshType::PointsContainerIterator It =
       this->GetOutput()->GetPoints()->Begin();

     while( It != this->GetOutput()->GetPoints()->End() )
       {
       this->GetOutput()->SetPointData( It.Index(), 
         NumericTraits<PixelType>::Zero );
       ++It;
       }
    } 

  if( this->GetOutput()->GetNumberOfPoints() > 0 )
    {
    typename OutputMeshType::PointDataContainerIterator ItD =
      this->GetOutput()->GetPointData()->Begin();
    while( ItD != this->GetOutput()->GetPointData()->End() )
      {
      if( find( this->m_LabelSet.begin(), 
        this->m_LabelSet.end(), ItD.Value() ) == this->m_LabelSet.end() )
        {
        this->m_LabelSet.push_back( ItD.Value() );
        }
      ++ItD;
      }
    }
}



template<class TOutputMesh>
void
VtkPolyDataFileReader<TOutputMesh>
::ReadVTKFile()
{
  this->ReadVTKPolyDataFile();
}

template<class TOutputMesh>
void
VtkPolyDataFileReader<TOutputMesh>
::ReadPointsFromVTKFile()
{
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();

  std::ifstream inputFile( this->m_FileName.c_str() );

  std::string line;
  
  bool isBinary = false;

  while( !inputFile.eof() )
    {
    std::getline( inputFile, line );

    if (line.find( "BINARY" ) != std::string::npos )
      {
      isBinary = true;
      }

    if( line.find( "POINTS" ) != std::string::npos )
      {
      break;
      }
    }

  itkDebugMacro( "POINTS line" << line );

  std::string pointLine( line, strlen( "POINTS " ), line.length() );
  itkDebugMacro( "pointLine " << pointLine );

  int numberOfPoints = -1;

  if( sscanf( pointLine.c_str(),"%d",&numberOfPoints ) != 1 )
    {
    itkExceptionMacro( "ERROR: Failed to read numberOfPoints\n"
        "       pointLine = " << pointLine );
    return;
    }

  itkDebugMacro( "numberOfPoints = " << numberOfPoints );

  if( numberOfPoints < 1 )
    {
    itkExceptionMacro( "numberOfPoints < 1"
        << "       numberOfPoints = " << numberOfPoints );
    return;
    }

  outputMesh->GetPoints()->Reserve( numberOfPoints );

  //
  // Load the point coordinates into the itk::Mesh
  //
  PointType point;

  if (isBinary)
    {
    itkDebugMacro( "Data is binary" );
          
    float p;
    float * ptData = new float [ numberOfPoints*3 ];
    inputFile.read( reinterpret_cast< char * >( ptData ), 3 * numberOfPoints * sizeof(p) );
    ByteSwapper<float>::SwapRangeFromSystemToBigEndian(ptData,numberOfPoints*3); 
    
    for (long i = 0; i < numberOfPoints; i++ )
      {
      for (long j = 0; j < Dimension; j++ )
        {
        point[j] = ptData[i*3+j];
        }
      outputMesh->SetPoint( i, point );
      }
      
    delete [] ptData;
    }
  else 
    {
    for( long i = 0; i < numberOfPoints; i++ )
      {
      if( Dimension == 2 )
        {
        float trash;
        inputFile >> point >> trash;
        }
      else  // Dimension = 3
        {
        inputFile >> point;
        }
      outputMesh->SetPoint( i, point );
      }
    }

  inputFile.close();
}


template<class TOutputMesh>
typename VtkPolyDataFileReader<TOutputMesh>::MultiComponentScalarSetType *
VtkPolyDataFileReader<TOutputMesh>
::GetMultiComponentScalarSetByName( const char * name )
{

  for (unsigned int i=0; i<this->m_MultiComponentScalarSetNames->Size(); i++)
  {
    if (!(this->m_MultiComponentScalarSetNames->GetElement(i).compare(name)))
    {
      return this->m_MultiComponentScalarSets->GetElement(i).GetPointer();
    }   
  }
  
  std::cout << "WARNING: No point data found with name: " << name << std::endl;
  
  return NULL;

}




template<class TOutputMesh>
bool 
VtkPolyDataFileReader<TOutputMesh>
::ReadVTKPolyDataFile( )
{
  
  this->m_InputFile.open( this->m_FileName.c_str(), std::ifstream::in );
  if (!this->m_InputFile )
  {
    std::cerr << "Unable to open file: " << this->m_FileName << std::endl;
    return false;
  }
  
  //typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  
  std::string line;
  std::getline(this->m_InputFile, line);
  
  if (line.find( "# vtk DataFile Version" ) == std::string::npos)
  {
    std::cerr << "Invalid header: " << line << std::endl;
    return false;
  }
  
  // Get header info line (ignore for now)
  std::getline(this->m_InputFile, line);
  
  // Get ASCII || BINARY flag
  std::getline(this->m_InputFile, line);
  if (line.find("BINARY") != std::string::npos)
  {
    this->m_BinaryData = true;
  }
  else if (line.find("ASCII") != std::string::npos)
  {
    this->m_BinaryData = false;
  }
  else
  {
    std::cerr << "Unknown data format: " << line << std::endl;
    return false;
  }
  
  // Get Dataset type
  std::getline(this->m_InputFile, line);
  if (line.find("DATASET POLYDATA") != std::string::npos)
  {
    return this->ReadVTKPolyData( );
  }
  else
  {
    std::cerr << "Only POLYDATA currently supported" << std::endl;
    return false;
  }
  
  return false;

  
}

template<class TOutputMesh>
bool 
VtkPolyDataFileReader<TOutputMesh>
::ReadVTKPolyData( )
{
  std::string line;
  
  
  while (!this->m_InputFile.eof())
  {
    std::getline(this->m_InputFile, line);
  
    if (line.find("POINTS") !=  std::string::npos)
    {
      
      std::string::size_type pos = line.rfind( " " );
      std::string dataType = std::string( line, pos+1, line.length()-1 );
      
      pos = line.find( " " );
      std::string::size_type pos2 = line.find( " ", pos+1 );
      std::string temp = std::string( line, pos+1, (pos2-pos) );
      unsigned long nPoints = std::atoi( temp.c_str() );
      
      if (!this->ReadVTKPoints(nPoints, dataType))
      {
        return false;
      }
    }
    else if (line.find("LINES") != std::string::npos)
    {
      std::string::size_type sp1 = line.find( " " );
      std::string::size_type sp2 = line.find( " ", sp1+1);
      
      unsigned long nLines = std::atoi( std::string( line, sp1+1, (sp2-sp1) ).c_str() );
      unsigned long nValues = std::atoi( std::string( line, sp2+1, line.length()-sp2 ).c_str()  );
      
      if (!this->ReadVTKLines( nLines, nValues ))
      {
        return false;
      } 
    }
    else if (line.find("POLYGONS") != std::string::npos)
    {
      
      std::string::size_type sp1 = line.find( " " );
      std::string::size_type sp2 = line.find( " ", sp1+1);
      
      unsigned long nPolys = std::atoi( std::string( line, sp1+1, (sp2-sp1) ).c_str() );
      unsigned long nValues = std::atoi( std::string( line, sp2+1, line.length()-sp2 ).c_str()  );
      
      if (!this->ReadVTKPolygons( nPolys, nValues ))
      {
        return false;
      }  
      
      //std::cerr << "POLYGONS not yet supported" << std::endl;
    }
    else if (line.find("VERTICES") != std::string::npos)
    {
      std::cerr << "VERTICES not yet supported" << std::endl;
    }
    else if (line.find("TRIANGLE_STRIPS") != std::string::npos)
    {
      std::cerr << "TRIANGE_STRIPS not yet supported" << std::endl; 
    }
    else if (line.find("POINT_DATA") != std::string::npos)
    {
      std::string::size_type sp1 = line.find( " " );
      
      unsigned long nPoints = std::atoi( std::string( line, sp1, line.length()-1 ).c_str() );
      
      if (!this->ReadVTKPointData( nPoints ))
      {
        return false;
      } 
    }
    else if (line.find("CELL_DATA") != std::string::npos)
    {
      std::string::size_type sp1 = line.find( " " );
      
      unsigned long nCells = std::atoi( std::string( line, sp1, line.length()-1 ).c_str() );
      
      if (!this->ReadVTKCellData( nCells ))
      {
        return false;
      }  

    }
    
  }
  
  this->m_InputFile.close();
  return true;
  
}


template<class TOutputMesh>
bool 
VtkPolyDataFileReader<TOutputMesh>
::ReadVTKPoints( unsigned long nPoints, std::string dataType )
{

  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  outputMesh->GetPoints()->Initialize();
  outputMesh->GetPoints()->Reserve(nPoints);
  
  typename OutputMeshType::PointDataContainer::Pointer pointData = OutputMeshType::PointDataContainer::New();
  pointData->Initialize();
  outputMesh->SetPointData( pointData );
  
  //std::cout << "Reading " << nPoints << " points of type " << dataType << std::endl;
  
  if (!this->m_BinaryData)
  {
    for (unsigned long i=0; i<nPoints; i++)
    {
      float x,y,z;
      this->m_InputFile >> x >> y >> z;
      typename OutputMeshType::PointType point;
      point[0] = x;
      if (point.Size() > 1)
        point[1] = y;
      if (point.Size() > 2)
        point[2] = z;
    
      outputMesh->GetPoints()->InsertElement(i,point);
        
    }
  }
  else
  {
    itkDebugMacro( "Data is binary" );
          
    float p;
    float * ptData = new float [ nPoints*3 ];
    this->m_InputFile.read( reinterpret_cast< char * >( ptData ), 3 * nPoints * sizeof(p) );
    ByteSwapper<float>::SwapRangeFromSystemToBigEndian(ptData,nPoints*3); 
    
    for (unsigned long i = 0; i < nPoints; i++ )
    {
      typename OutputMeshType::PointType point;
      for (long j = 0; j < 3; j++ )
      {
        if (j < OutputMeshType::PointType::Dimension) 
          point[j] = ptData[i*3+j];
      }

      outputMesh->SetPoint( i, point );
    }
      
    delete [] ptData;    
    
    
  }
  
  return true;
}

template<class TOutputMesh>
bool 
VtkPolyDataFileReader<TOutputMesh>
::ReadVTKLines( unsigned long nLines, unsigned long nValues )
{

  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  this->m_Lines = LineSetType::New();
  
  
  //std::cout << "Reading " << nLines << " lines with " << nValues << " values" << std::endl;
  
  
  if (!this->m_BinaryData)
  {
    for (unsigned long i=0; i<nLines; i++)
    {
      
     
      unsigned long lineSize;
      this->m_InputFile >> lineSize;
      
      LineType line(lineSize);
           
      for (unsigned long j=0; j<lineSize; j++)
      {
        unsigned long index;
        this->m_InputFile >> index;
        line[j] = index;
      }
      
      this->m_Lines->InsertElement(i,line);  
        
    }
  }
  else
  {
    
    int p;
    int * lineData = new int [ nValues ];
    this->m_InputFile.read( reinterpret_cast< char * >( lineData ), nValues * sizeof(p) );
    ByteSwapper<int>::SwapRangeFromSystemToBigEndian(lineData,nValues); 
    
    unsigned long valueId = 0;
    unsigned long lineId = 0;
    while (valueId < nValues)
    {
      unsigned int lineLength = lineData[valueId];
      ++valueId;
      
      LineType polyLine;
      polyLine.SetSize( lineLength );
        
      for (unsigned long i = 0; i < lineLength; i++)
        {
        polyLine[i] = lineData[valueId];
        ++valueId;
        }
        
      this->m_Lines->InsertElement( lineId, polyLine );
      ++lineId;
    }  
    
    delete [] lineData;    
    
    
  }
  
  return true;
}



template<class TOutputMesh>
bool 
VtkPolyDataFileReader<TOutputMesh>
::ReadVTKPolygons( unsigned long nPolygons, unsigned long nValues )
{

  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  this->m_Polygons = LineSetType::New();
  
  
  //std::cout << "Reading " << nPolygons << " polygons with " << nValues << " values" << std::endl;
  
  
  if (!this->m_BinaryData)
  {
    for (unsigned long i=0; i<nPolygons; i++)
    {
      
     
      unsigned long lineSize;
      this->m_InputFile >> lineSize;
      
      LineType line(lineSize);
           
      for (unsigned long j=0; j<lineSize; j++)
      {
        unsigned long index;
        this->m_InputFile >> index;
        line[j] = index;
      }
      
      this->m_Polygons->InsertElement(i,line);  
        
    }
  }
  else
  {
    
    int p;
    int * lineData = new int [ nValues ];
    this->m_InputFile.read( reinterpret_cast< char * >( lineData ), nValues * sizeof(p) );
    ByteSwapper<int>::SwapRangeFromSystemToBigEndian(lineData,nValues); 
    
    unsigned long valueId = 0;
    unsigned long lineId = 0;
    while (valueId < nValues)
    {
      unsigned int lineLength = lineData[valueId];
      ++valueId;
      
      LineType polyLine;
      polyLine.SetSize( lineLength );
        
      for (unsigned long i = 0; i < lineLength; i++)
        {
        polyLine[i] = lineData[valueId];
        ++valueId;
        }
        
      this->m_Polygons->InsertElement( lineId, polyLine );
      ++lineId;
    }  
    
    delete [] lineData;    
    
    
  }
  
  return true;
}



template<class TOutputMesh>
bool 
VtkPolyDataFileReader<TOutputMesh>
::ReadVTKPointData( unsigned long nPoints )
{

  std::string line;
  std::getline(this->m_InputFile, line);
  
  this->m_MultiComponentScalarSets = MultiComponentScalarMultiSetType::New();
  this->m_MultiComponentScalarSets->Initialize();
  
  this->m_MultiComponentScalarSetNames = MultiComponentScalarSetNamesType::New();
  this->m_MultiComponentScalarSetNames->Initialize();
  
  while (line.length() > 0)
  {
    if (line.find("SCALARS") != std::string::npos)
    {
      
      unsigned long nComponents = 1;
      std::string::size_type sp1 = line.find( " " );
      std::string::size_type sp2 = line.find( " ", sp1+1);
      std::string::size_type sp3 = line.find( " ", sp2+1);
      if (sp3 == std::string::npos)
      {
        sp3 = line.length()-1;
      }
      else
      {
        nComponents = std::atoi( std::string( line, sp3+1, line.length()-sp2-1 ).c_str()  );
      }
      
      std::string dataName = std::string( line, sp1+1, (sp2-sp1)-1 );
      std::string dataType = std::string( line, sp2+1, (sp3-sp2)-1 );
      
      this->ReadVTKScalars(dataName, dataType, nPoints, nComponents, true);
      
    }
    else if (line.find("COLOR_SCALARS") != std::string::npos)
    {
      std::cerr << "COLOR_SCALARS not yet supported" << std::endl;
    }
    else if (line.find("LOOKUP_TABLE") != std::string::npos)
    {
      std::cerr << "LOOKUP_TABLE not yet supported" << std::endl;
    }
    else if (line.find("VECTORS") != std::string::npos)
    {
      std::cerr << "VECTORS not yet supported" << std::endl;
    }
    else if (line.find("NORMALS") != std::string::npos)
    {
      std::cerr << "NORMALS not yet supported" << std::endl;
    }
    else if (line.find("TEXTURE_COORDINATES") != std::string::npos)
    {
      std::cerr << "TEXTURE_COORDINATES not yet supported" << std::endl;
    }
    else if (line.find("TENSORS") != std::string::npos)
    {
      std::cerr << "TENSORS not yet supported" << std::endl;
    }
    else if (line.find("FIELD") != std::string::npos)
    {
      std::cerr << "FIELD not yet supported" << std::endl;
    }

    std::getline( this->m_InputFile, line );
  }
  
  return true;
}

template<class TOutputMesh>
bool 
VtkPolyDataFileReader<TOutputMesh>
::ReadVTKCellData( unsigned long nCells )
{
  
  //std::cout << "Reading cell data for " << nCells << " cells" << std::endl;

  std::string line;
  std::getline(this->m_InputFile, line);
  
  this->m_CellMultiComponentScalarSets = MultiComponentScalarMultiSetType::New();
  this->m_CellMultiComponentScalarSets->Initialize();
  
  this->m_CellMultiComponentScalarSetNames = MultiComponentScalarSetNamesType::New();
  this->m_CellMultiComponentScalarSetNames->Initialize();
  
  while (line.length() > 0)
  {
    if (line.find("SCALARS") != std::string::npos)
    {
      
      unsigned long nComponents = 1;
      std::string::size_type sp1 = line.find( " " );
      std::string::size_type sp2 = line.find( " ", sp1+1);
      std::string::size_type sp3 = line.find( " ", sp2+1);
      if (sp3 == std::string::npos)
      {
        sp3 = line.length()-1;
      }
      else
      {
        nComponents = std::atoi( std::string( line, sp3+1, line.length()-sp2-1 ).c_str()  );
      }
      
      std::string dataName = std::string( line, sp1+1, (sp2-sp1)-1 );
      std::string dataType = std::string( line, sp2+1, (sp3-sp2)-1 );
      
      this->ReadVTKScalars(dataName, dataType, nCells, nComponents, false);
      
    }
    else if (line.find("COLOR_SCALARS") != std::string::npos)
    {
      std::cerr << "COLOR_SCALARS not yet supported" << std::endl;
    }
    else if (line.find("LOOKUP_TABLE") != std::string::npos)
    {
      std::cerr << "LOOKUP_TABLE not yet supported" << std::endl;
    }
    else if (line.find("VECTORS") != std::string::npos)
    {
      std::cerr << "VECTORS not yet supported" << std::endl;
    }
    else if (line.find("NORMALS") != std::string::npos)
    {
      std::cerr << "NORMALS not yet supported" << std::endl;
    }
    else if (line.find("TEXTURE_COORDINATES") != std::string::npos)
    {
      std::cerr << "TEXTURE_COORDINATES not yet supported" << std::endl;
    }
    else if (line.find("TENSORS") != std::string::npos)
    {
      std::cerr << "TENSORS not yet supported" << std::endl;
    }
    else if (line.find("FIELD") != std::string::npos)
    {
      std::cerr << "FIELD not yet supported" << std::endl;
    }

    std::getline( this->m_InputFile, line );
  }
  
  return true;
}


template<class TOutputMesh>
bool 
VtkPolyDataFileReader<TOutputMesh>
::ReadVTKScalars( std::string dataName, std::string dataType, unsigned long nPoints, unsigned long nComponents, bool isPointData=true )
{
 
  //std::cout << "Reading " << nPoints << " scalar data points with " << nComponents << " components of type " << dataType << " named " << dataName << std::endl;
  
  // get LOOKUP_TABLE name
  std::string line;
  std::getline(this->m_InputFile, line);
  
  bool labels = false;
  if (dataName.find("pointLabels") != std::string::npos)
  {
    labels = true;
    
  }
  
  typename MultiComponentScalarSetType::Pointer set = MultiComponentScalarSetType::New();
  set->Initialize();
  set->Reserve( nPoints );
  
  if (!this->m_BinaryData)
  {
    for (unsigned long i=0; i<nPoints; i++)
    {
      MultiComponentScalarType value(nComponents);
      for (unsigned long j=0; j<nComponents; j++)
      {
        float component;
        this->m_InputFile >> component;
        value[j] = component;
      }
      if (labels)
      {
        typename OutputMeshType::PixelType label = static_cast<typename OutputMeshType::PixelType>( value[0] );
        this->GetOutput()->GetPointData()->InsertElement(i, label);
      }
      else
      {
        set->InsertElement(i, value);
      }
    }
  }
  else
  {
        
  }
  
  if (!labels)
  {
    if (isPointData)
      {
      long id = this->m_MultiComponentScalarSets->Size();
      this->m_MultiComponentScalarSets->InsertElement( id, set );
      this->m_MultiComponentScalarSetNames->InsertElement( id, dataName );
    }
    else
      {
      long id = this->m_CellMultiComponentScalarSets->Size();
      this->m_CellMultiComponentScalarSets->InsertElement( id, set );
      this->m_CellMultiComponentScalarSetNames->InsertElement( id, dataName );  
      }
  }
  return true;
}



template<class TOutputMesh>
void
VtkPolyDataFileReader<TOutputMesh>
::ReadScalarsFromVTKFile( )
{
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();

  std::ifstream inputFile( this->m_FileName.c_str() );

  std::string line;

  bool isBinary = false;

  //
  // Find the labels associated with each pixel
  //
  while( !inputFile.eof() )
    {
    std::getline( inputFile, line );

   if (line.find( "BINARY" ) != std::string::npos )
      {
      isBinary = true;
      }

    if( line.find( "SCALARS" ) != std::string::npos )
      {
      break;
      }
    }

  if( inputFile.eof() )
    {
    inputFile.close();
    return;
    }


  std::string::size_type pos = line.rfind( " " );
  std::string temp = std::string( line, pos+1, line.length()-1 );
  unsigned int numberOfComponents = std::atoi( temp.c_str() );

  pos = line.find( " " );
  std::string::size_type pos2 = line.find( " ", pos+1 );
  temp = std::string( line, pos+1, (pos2-pos) );
  //std::cout << "Reading data named: " << temp << std::endl;

  if (isBinary)
    {
 
  
    int numberOfValues = outputMesh->GetNumberOfPoints()*numberOfComponents;
    float p;
    int * scalarData = new int [ numberOfValues ];
    inputFile.read( reinterpret_cast< char * >( scalarData ), numberOfComponents * sizeof(p) );
    ByteSwapper<int>::SwapRangeFromSystemToBigEndian(scalarData,numberOfValues); 
  
    if( numberOfComponents == 1 )
      {
      //PixelType label;
      for( unsigned long i = 0; i < outputMesh->GetNumberOfPoints(); i++ )
        {
        outputMesh->SetPointData( i, scalarData[i] );
        }
  //    itkExceptionMacro( "Only single label components are readable" );
      }
    else
      {
        
      this->m_MultiComponentScalarSets = MultiComponentScalarMultiSetType::New(); 
      this->m_MultiComponentScalarSets->Initialize();
      
      typename MultiComponentScalarSetType::Pointer set = MultiComponentScalarSetType::New();
      set->Initialize();
      
      
      for( unsigned long i = 0; i < outputMesh->GetNumberOfPoints(); i++ )
        {
        MultiComponentScalarType scalar;
        scalar.SetSize( numberOfComponents );
        for( unsigned int d = 0; d < numberOfComponents; d++ )
          {
          scalar[d] = scalarData[i*numberOfComponents + d];
          }
        set->InsertElement( i, scalar );
        }
      } 
    
      delete [] scalarData;
      //this->m_MultiComponentScalarSets->InsertElement(0,set);
       
    }
  else
    {
      
      

    while (numberOfComponents > 0)
      {

      if( numberOfComponents == 1 )
        {
        PixelType label;
        for( unsigned long i = 0; i < outputMesh->GetNumberOfPoints(); i++ )
          {
          inputFile >> label;
          outputMesh->SetPointData( i, label );
          }
    //    itkExceptionMacro( "Only single label components are readable" );
        }
      else
        {
        this->m_MultiComponentScalarSets = MultiComponentScalarMultiSetType::New(); 
        this->m_MultiComponentScalarSets->Initialize();
        
        for( unsigned long i = 0; i < outputMesh->GetNumberOfPoints(); i++ )
          {
          MultiComponentScalarType scalar;
          scalar.SetSize( numberOfComponents );
          for( unsigned int d = 0; d < numberOfComponents; d++ )
            {
            inputFile >> scalar[d];
            }
          //this->m_MultiComponentScalars->InsertElement( i, scalar );
          }
        }
        
        numberOfComponents = 0;
        
      }
      
    }

  inputFile.close();
}

template<class TOutputMesh>
void
VtkPolyDataFileReader<TOutputMesh>
::ReadLinesFromVTKFile()
{
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();

  std::ifstream inputFile( this->m_FileName.c_str() );

  std::string line;
  
  bool isBinary = false;

  //
  // Find the labels associated with each pixel
  //
  while( !inputFile.eof() )
    {
    std::getline( inputFile, line );
    
    if (line.find( "BINARY" ) != std::string::npos )
      {
      isBinary = true;
      }

    if( line.find( "LINES" ) != std::string::npos )
      {
      break;
      }
    }

  if( inputFile.eof() )
    {
    inputFile.close();
    return;
    }

  std::string::size_type pos = line.rfind( " " );

  std::string temp = std::string( line, 6, pos-1 );
  unsigned int numberOfLines = std::atoi( temp.c_str() );

  temp = std::string(line, pos, line.length()-1 );
  unsigned int numberOfValues = std::atoi( temp.c_str() );

  this->m_Lines = LineSetType::New(); 
  this->m_Lines->Initialize();
  
  if (isBinary)
    {
    int p;
    int * lineData = new int [ numberOfValues ];
    inputFile.read( reinterpret_cast< char * >( lineData ), numberOfValues * sizeof(p) );
    ByteSwapper<int>::SwapRangeFromSystemToBigEndian(lineData,numberOfValues); 
    
    long valueId = 0;
    long lineId = 0;
    while (valueId < numberOfValues)
    {
      int lineLength = lineData[valueId];
      ++valueId;
      
      LineType polyLine;
      polyLine.SetSize( lineLength );
        
      for (long i = 0; i < lineLength; i++)
        {
        polyLine[i] = lineData[valueId];
        ++valueId;
        }
      this->m_Lines->InsertElement( lineId, polyLine );
      ++lineId;
    }  
    
    delete [] lineData;
  }
  else 
  {
    
    
    for( unsigned int i = 0; i < numberOfLines; i++ )
      {
      LineType line;
      unsigned int numberOfPoints;
      inputFile >> numberOfPoints;
      line.SetSize( numberOfPoints );
      for( unsigned int d = 0; d < numberOfPoints; d++ )
        {
        inputFile >> line[d];
        }
      this->m_Lines->InsertElement( i, line );
      }
    }

  inputFile.close();
}


template<class TOutputMesh>
void
VtkPolyDataFileReader<TOutputMesh>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "FileName: "
     << this->m_FileName << std::endl;
  os << indent << "RandomPercentage: "
     << this->m_RandomPercentage << std::endl;

}

}
} //end of namespace itk


#endif
