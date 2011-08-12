/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit (ITK)
  Module:    $RCSfile: itkGeodesicPaths.hxx,v $
  Language:  C++
  Date:      $Date: 2009/10/20 20:16:43 $
  Version:   $Revision: 1.1 $

=========================================================================*/
#ifndef _itkGeodesicPaths_cxx_
#define _itkGeodesicPaths_cxx_

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "itkGeodesicPaths.h"

namespace itk {

/*
 * compute the local cost
 */
template<class TGraphSearchNode, typename TCostImageType > 
typename GeodesicPaths<TGraphSearchNode,TCostImageType>::
PixelType GeodesicPaths<TGraphSearchNode,TCostImageType>::LocalCost() 
{
  if (!m_LocalCostFP) {
    NodeLocationType loc1=m_CurrentNode->GetLocation();
    NodeLocationType loc2=m_NeighborNode->GetLocation();
//    NodeLocationType stepsize=loc2-loc1;
    float cost=0;
    for (int i=0; i<GraphDimension; i++) 
      cost+=(loc2[i]-loc1[i])*(loc2[i]-loc1[i]);
    cost=sqrt(cost);
    return cost;
//m_CostIndex[i]=loc2[i];
//    return stepsize.magnitude();
    //*(1. + (PixelType) m_GraphCostImage->GetPixel(m_CostIndex));
  }
  else return m_LocalCostFP(this);
};



template<class TGraphSearchNode, typename TCostImageType >
void GeodesicPaths<TGraphSearchNode,TCostImageType>::SearchEdgeSet() 
// loops over the neighbors in the graph 
{
  int i=0,j=0;  
  GraphCostImagePixelType pix=0;
  NodeLocationType loc,oloc,dif;	

  GraphNeighborhoodIteratorType GHood(m_Radius, m_Graph,m_Graph->GetLargestPossibleRegion());
  GraphNeighborhoodIndexType	GNI;

  bool nearedge=false;
  for (i=0; i < GraphDimension; i++)
  {   
  	GNI[i]=(long int) m_CurrentNode->GetLocation()[i];
    oloc[i]=(float)GNI[i];
    if (GNI[i] > m_GraphSize[i]-1 || GNI[i] < 1) nearedge=true; 
  } 

  m_MeanGeodesic+=m_CurrentCost;

  if (nearedge) return;

  GHood.SetLocation(GNI);
  for (i = 0; i < GHood.Size(); i++)
  {     
  	m_NeighborNode=GHood.GetPixel(i);
    m_GraphIndex = GHood.GetIndex(i); 
    pix=m_GraphCostImage->GetPixel(m_GraphIndex);

    if (m_NeighborNode != m_CurrentNode && pix == m_SurfaceLabel)
    { 
      if (!m_NeighborNode)
      {   

//   std::cout << " new ind " << loc << " orig " << oloc << std::endl;
      	typename GraphSearchNode<PixelType,CoordRep,GraphDimension>::Pointer G=
  	  	GraphSearchNode<PixelType,CoordRep,GraphDimension>::New();
    	  G->SetUnVisited();
  	    G->SetTotalCost(m_MaxCost);
        for (int k=0; k<GraphDimension; k++) loc[k]=m_GraphIndex[k];
    	  G->SetLocation(loc);
  	    G->SetPredecessor(NULL);
  	    m_Graph->SetPixel(m_GraphIndex,G);

        m_NeighborNode=G;

        //dif=loc-oloc;

      }
      TerminationCondition();
    
	    if (!m_SearchFinished && m_CurrentNode != m_NeighborNode &&
		    !m_NeighborNode->GetDelivered()  )
//        && dif.magnitude() <= 1.0)
  	  {
          m_NewCost = m_CurrentCost + LocalCost();
//        std::cout << "  entering " << (int) pix << " cost "<< m_NewCost << std::endl;
	        CheckNodeStatus();
	    }
    }
  }
}


template<class TGraphSearchNode, typename TCostImageType > 
void 
GeodesicPaths<TGraphSearchNode,TCostImageType>
::PathVariation() 
{
  float maxunitvecdiff=0.0;
  float mag=0.0;
  GraphIteratorType Iter( m_Graph, m_GraphRegion );
  Iter.GoToBegin();
  
  for( ; !Iter.IsAtEnd(); ++Iter  )
  {
    m_GraphIndex = Iter.GetIndex(); 
    m_CurrentNode = Iter.Get();
    unsigned char curpix=
    m_GraphCostImage->GetPixel(m_GraphIndex);
    if (m_CurrentNode && curpix == m_SurfaceLabel)
    {                 
      int i=0,j=0;  
      GraphCostImagePixelType pix=0;

      // keeps the current unit vector defined by the predecessor
      NodeLocationType curunitvec;	
      // keeps the neighbor unit vector defined by the neighbor predecessor
      NodeLocationType nunitvec;	
      NodeLocationType unitvecdif;	

      GraphNeighborhoodIteratorType GHood(m_Radius, m_Graph,m_Graph->GetRequestedRegion());
      GraphNeighborhoodIndexType	GNI;

      for (i=0; i < GraphDimension; i++)
      {   
  	    GNI[i]=m_CurrentNode->GetLocation()[i];
        curunitvec[i]=m_CurrentNode->GetLocation()[i]
          -m_CurrentNode->GetPredecessor()->GetLocation()[i];
      } 

      mag=curunitvec.magnitude();
      if (mag != 0) curunitvec/=mag;

      // keeps track of the difference between cur and neigh unit vec
      maxunitvecdiff=0.0;
      mag=0.0;

      GHood.SetLocation(GNI);
      for (i = 0; i < GHood.Size(); i++)
      {     
  	    m_NeighborNode=GHood.GetPixel(i);
        m_GraphIndex = GHood.GetIndex(i); 
        pix=m_GraphCostImage->GetPixel(m_GraphIndex);

        if (m_NeighborNode && m_NeighborNode != m_CurrentNode && pix == m_SurfaceLabel)
        { 
          for (j=0; j < GraphDimension; j++)
          {   
            nunitvec[j]=m_NeighborNode->GetLocation()[j]
               -m_NeighborNode->GetPredecessor()->GetLocation()[j];
          }
          //mag=nunitvec.magnitude();
          //if (mag !=0.0) nunitvec/=mag;

          //unitvecdif=nunitvec-curunitvec;
          //mag=unitvecdif.magnitude();
          mag += fabs(m_CurrentNode->GetTotalCost()-m_NeighborNode->GetTotalCost())/100.;
          //if (mag > maxunitvecdiff) 
          maxunitvecdiff=mag;
        }
      }
//      if (fabs(m_CurrentNode->GetTotalCost() - 200.0) < 2.5) 
//        m_CurrentNode->SetTotalCost(255);
//      if (fabs(m_CurrentNode->GetTotalCost() - 150.0) < 2.5) 
//        m_CurrentNode->SetTotalCost(195);
//      else 
        m_CurrentNode->SetTotalCost(maxunitvecdiff);
    }
  }

}



template<class TGraphSearchNode, typename TCostImageType > 
void 
GeodesicPaths<TGraphSearchNode,TCostImageType>
::CopyTotalCostToCostImage() 
{
  typedef typename GraphCostImageType::PixelType PixelType;
  GraphIteratorType Iter( m_Graph, m_GraphRegion );
  Iter.GoToBegin();
  
  for( ; !Iter.IsAtEnd(); ++Iter  )
  {
    m_GraphIndex = Iter.GetIndex(); 
    m_CurrentNode = Iter.Get();
    TCostImageType curpix=
      m_GraphCostImage->GetPixel(m_GraphIndex);
    if (m_CurrentNode && curpix == m_SurfaceLabel)
    {   
      m_GraphCostImage->SetPixel(m_GraphIndex,
       (PixelType) m_CurrentNode->GetTotalCost() );
    }
  }

}


template<class TGraphSearchNode, typename TCostImageType > 
void 
GeodesicPaths<TGraphSearchNode,TCostImageType>
::MapEndCostToSource() 
{
  typedef typename GraphCostImageType::PixelType PixelType;
  GraphIteratorType Iter( m_Graph, m_GraphRegion );

  Iter.GoToBegin();
  for( ; !Iter.IsAtEnd(); ++Iter  )
  {
    m_GraphIndex = Iter.GetIndex(); 
    m_CurrentNode = Iter.Get();
    TCostImageType curpix=
      m_GraphCostImage->GetPixel(m_GraphIndex);
    if (m_CurrentNode )
    {   
      this->BackTrack(m_CurrentNode);
       m_GraphCostImage->SetPixel(m_GraphIndex,10);
    }
  }
  


  Iter.GoToBegin();
  for( ; !Iter.IsAtEnd(); ++Iter  )
  {
    m_GraphIndex = Iter.GetIndex(); 
    m_CurrentNode = Iter.Get();
     TCostImageType curpix=
      m_GraphCostImage->GetPixel(m_GraphIndex);
    if (m_CurrentNode) 
    //&& curpix==m_SourceLabel)
    {   
    
      this->ForwardTrack(m_CurrentNode);

      m_GraphCostImage->SetPixel(m_GraphIndex,m_CurrentNode->GetValue()+1.0);
//        (TCostImageType) m_CurrentNode->GetValue() );
    }
//    else m_GraphCostImage->SetPixel(m_GraphIndex,0 );
  }

}



template<class TGraphSearchNode, typename TCostImageType > 
typename GeodesicPaths<TGraphSearchNode,TCostImageType>::FieldTypePointer
GeodesicPaths<TGraphSearchNode,TCostImageType>
::MapToVectorField() 
{

  typename FieldType::Pointer    field = FieldType::New();
  field->SetLargestPossibleRegion(  m_GraphCostImage->GetLargestPossibleRegion() ); 
  field->SetRequestedRegion(   m_GraphCostImage->GetLargestPossibleRegion() );  
  field->SetBufferedRegion(    m_GraphCostImage->GetLargestPossibleRegion());
  field->Allocate(); 

  itk::ImageRegionIteratorWithIndex<FieldType> 
    fit( field, field->GetLargestPossibleRegion() );
  
  VectorType zero; 
  zero.Fill(0.0);
  fit.GoToBegin();
  while(!fit.IsAtEnd()) {fit.Set(zero); ++fit;}
 

  typedef typename GraphCostImageType::PixelType PixelType;
  GraphIteratorType Iter( m_Graph, m_GraphRegion );

  Iter.GoToBegin();
  for( ; !Iter.IsAtEnd(); ++Iter  )
  {
    m_GraphIndex = Iter.GetIndex(); 
    m_CurrentNode = Iter.Get();


    if (m_CurrentNode )
    {   

    typename TGraphSearchNode::Pointer G=m_CurrentNode;
    typename TGraphSearchNode::Pointer P=m_CurrentNode->GetPredecessor();
    
    float highcost=G->GetTotalCost();
   
    VectorType vec;
    for (int ii=0; ii<GraphDimension;ii++) vec[ii]=(float)0;
    if (highcost > G->GetValue()) 
    {
    field->SetPixel(m_GraphIndex,vec);
    G->SetValue(highcost);
    }
    typename FieldType::IndexType nextind;
    
    while(P && G != P)
    {
    
      for (int ii=0; ii<GraphDimension;ii++) 
       {
       nextind[ii]=(long)(P->GetLocation()[ii]+0.5);
       vec[ii]=(float)(m_GraphIndex[ii]-nextind[ii]);
       }
//       std::cout << " Vec " << vec << std::endl;
      if (highcost > P->GetValue()) 
      {
      field->SetPixel(nextind,vec);
      P->SetValue(highcost);
      m_GraphCostImage->SetPixel(nextind,highcost);
      }

	    G=P;
	    P=G->GetPredecessor();
      
    }

    }
  }

  return field;
  

}




template<class TGraphSearchNode, typename TCostImageType > 
void 
GeodesicPaths<TGraphSearchNode,TCostImageType>
::ConnectedSurface(float frac) 
{


// collect all the indices of points that are in the surface

  typedef std::vector<CostImageIndexType>      indexContainerType;
  indexContainerType  surfaceindexlist;

  m_GraphSize=m_GraphCostImage->GetLargestPossibleRegion().GetSize();

  this->InitializeGraph();

  unsigned long numpts=0;
  unsigned long ii=0,mct=0;

  typedef ImageRegionIteratorWithIndex<GraphCostImageType> CostImageIteratorType; 
  GraphCostImagePixelType pix=0;
  CostImageIteratorType CostIterator( m_GraphCostImage, m_GraphRegion );
  CostIterator.GoToBegin();
  while(  !CostIterator.IsAtEnd()  )
  {
    pix=CostIterator.Get();
    
    if (pix == m_SurfaceLabel) 
    {
      surfaceindexlist.insert(surfaceindexlist.end(),CostIterator.GetIndex());
      numpts++;
    }
    ++CostIterator; 
  }

  if (numpts == 0)
  { 
    std::cout << " No surface points ";
    return;
  }
// store the number of these points - use e.g. 10% for randomized alg
  srand ( time(NULL) );

  bool done = false;
  float fracvisited = 0;
  while (!done )
  {

  // first initialize the graph
    this->InitializeGraph();

    unsigned long randind = rand() % numpts;

  // select one at random
    NodeLocationType pos1;
    GraphNeighborhoodIndexType	index=surfaceindexlist[randind];

    for (unsigned int t=0; t<GraphDimension; t++)
      pos1[t]=index[t];

    typename SearchNode::Pointer G1 = SearchNode::New();
    G1->SetLocation(pos1);
    this->SetSource(G1);	
    this->InitializeQueue();
    this->SetRadius(1);  
    this->FindPath(); 
    
    fracvisited = (float)m_QS->GetTimer()/(float)numpts;
    if ( fracvisited > frac || mct > 100)
    {  
  	  GraphIteratorType Iter( m_Graph, m_GraphRegion );
      Iter.GoToBegin();
  
      for( ; !Iter.IsAtEnd(); ++Iter  )
	    {
  	    m_GraphIndex = Iter.GetIndex(); 
        m_CurrentNode = Iter.Get();
        unsigned char curpix=
          m_GraphCostImage->GetPixel(m_GraphIndex);
        if (!m_CurrentNode && curpix == m_SurfaceLabel)
        {                 
          m_GraphCostImage->SetPixel(m_GraphIndex,1);
        }
      }
      done = true;
    } 
    std::cout << " frac visited " << fracvisited << std::endl;

    this->EmptyQ();
    mct++;

  }

}




template<class TGraphSearchNode, typename TCostImageType > 
void 
GeodesicPaths<TGraphSearchNode,TCostImageType>
::LabelConnectedSurface(float frac) 
{


// collect all the indices of points that are in the surface

  typedef std::vector<CostImageIndexType>      indexContainerType;
  indexContainerType  surfaceindexlist;

  m_GraphSize=m_GraphCostImage->GetLargestPossibleRegion().GetSize();

  this->InitializeGraph();

  unsigned long numpts=0;
  unsigned long ii=0,mct=0;

  typedef ImageRegionIteratorWithIndex<GraphCostImageType> CostImageIteratorType; 
  GraphCostImagePixelType pix=0;
  CostImageIteratorType CostIterator( m_GraphCostImage, m_GraphRegion );
  CostIterator.GoToBegin();
  while(  !CostIterator.IsAtEnd()  )
  {
    pix=CostIterator.Get();
    
    if (pix == m_SurfaceLabel) 
    {
      surfaceindexlist.insert(surfaceindexlist.end(),CostIterator.GetIndex());
      numpts++;
    }
    ++CostIterator; 
  }

  if (numpts == 0)
  { 
    std::cout << " No surface points ";
    return;
  }
// store the number of these points - use e.g. 10% for randomized alg
  srand ( time(NULL) );

  bool done = false,innerdone=false;
  float fracvisited = 0, totalfrac=0.0;
  unsigned int label=m_SurfaceLabel+1;  // surface is labeled as 2
  

  CostIterator.GoToBegin();

  while (!done )
  {

    GraphNeighborhoodIndexType	index;
    innerdone=false;
    while(  !CostIterator.IsAtEnd() && !innerdone )
    {
      pix=CostIterator.Get();
    
      if (pix == m_SurfaceLabel) 
      { 
        index=CostIterator.GetIndex();
        innerdone=true;
      }
      ++CostIterator; 
    }

    // first initialize the graph
    this->InitializeGraph();


  // select one at random
    NodeLocationType pos1;

    for (unsigned int t=0; t<GraphDimension; t++)
      pos1[t]=index[t];

    typename SearchNode::Pointer G1 = SearchNode::New();
    G1->SetLocation(pos1);
    this->SetSource(G1);	
    this->InitializeQueue();
    this->SetRadius(1);  
    this->FindPath(); 
    
    fracvisited = (float)m_QS->GetTimer()/(float)numpts;
    std::cout << " index of start point " << index << std::endl;
    if ( m_QS->GetTimer() > 30 )
    {
      // label the visited positions with the label
  	  GraphIteratorType Iter( m_Graph, m_GraphRegion );
      Iter.GoToBegin();

      unsigned int uselabel= (unsigned int) m_QS->GetTimer();//rand();
      if (uselabel < 3 ) uselabel=uselabel+10;
  
      for( ; !Iter.IsAtEnd(); ++Iter  )
	    {
  	    m_GraphIndex = Iter.GetIndex(); 
        m_CurrentNode = Iter.Get();
        unsigned char curpix=
          m_GraphCostImage->GetPixel(m_GraphIndex);
        if (m_CurrentNode && curpix == m_SurfaceLabel)
        {                 
          m_GraphCostImage->SetPixel(m_GraphIndex,uselabel);
        }
      } 

      totalfrac+=fracvisited;
      label++;
      if (totalfrac > frac) done = true;
    } 
    else
    {  
      // label the visited positions with zero (delete them as they're too small)
  	  GraphIteratorType Iter( m_Graph, m_GraphRegion );
      Iter.GoToBegin();
  
      for( ; !Iter.IsAtEnd(); ++Iter  )
	    {
  	    m_GraphIndex = Iter.GetIndex(); 
        m_CurrentNode = Iter.Get();
        unsigned char curpix=
          m_GraphCostImage->GetPixel(m_GraphIndex);
        if (m_CurrentNode && curpix == m_SurfaceLabel)
        {                 
          m_GraphCostImage->SetPixel(m_GraphIndex,0);
        }
      } 
      totalfrac+=fracvisited;
      if (totalfrac > frac) done = true;
    } 
    
    std::cout << " frac visited " << totalfrac << " label " << label << std::endl;

    this->EmptyQ();
    mct++;
  }

}



template<class TGraphSearchNode, typename TCostImageType > 
void 
GeodesicPaths<TGraphSearchNode,TCostImageType>
::MeanGeodesic(unsigned long num ,float frac) 
{

// initialize the list of mean geodesic distances 
  typedef std::vector<float>      geodesicdistanceContainerType;
  geodesicdistanceContainerType   meandistancelist;



// collect all the indices of points that are in the surface

  typedef std::vector<CostImageIndexType>      indexContainerType;
  indexContainerType  surfaceindexlist;

  m_GraphSize=m_GraphCostImage->GetLargestPossibleRegion().GetSize();

  this->InitializeGraph();

  typedef ImageRegionIteratorWithIndex<GraphCostImageType> CostImageIteratorType; 
  GraphCostImagePixelType pix=0;
  CostImageIteratorType CostIterator( m_GraphCostImage, m_GraphRegion );
  CostIterator.GoToBegin();
  while(  !CostIterator.IsAtEnd()  )
  {
    pix=CostIterator.Get();
    
    if (pix == m_SurfaceLabel) 
    {
     surfaceindexlist.insert(surfaceindexlist.end(),CostIterator.GetIndex());
     meandistancelist.insert(meandistancelist.end(),0.0);
    }
    ++CostIterator; 
   }

// store the number of these points - use e.g. 10% for randomized alg

   unsigned long numpts=surfaceindexlist.size();

// estimate radius of equivalent sphere

   float rad = sqrt((float)numpts/(4.0*3.14159));
   std::cout << " estimated radius " << rad << std::endl; 

   float sphereG=3.14159/2.0*(float) rad;
   std::cout << " sphere G " <<  sphereG << std::endl;

   if (frac != 0) 
   {
    m_SphereGeodesic = sphereG;
    this->GenerateEquivalentSphere(numpts,rad);
   } 
   else sphereG = m_SphereGeodesic;
  
   unsigned long trynum = 0;
   if (num == 0) trynum =(unsigned long) ((float) numpts * frac);
   else trynum = num;

   if (trynum > numpts) trynum=numpts;

   unsigned long step = numpts/trynum;

   if (step == 0) step=1;

   std::cout << " random attempts " << trynum << " step " << step  << 
    " num " << num << " numpts " << numpts << std::endl;
   srand ( time(NULL) );


   float aggmeandist=0,minmeandist=1.e9;
   
   unsigned long ii=0,mct=0;

   for (ii=0; ii< numpts; ii=ii+step)
   {

  // first initialize the graph
    this->InitializeGraph();

//    float randmax=32767.0;
    unsigned long randind = ii; //rand() % numpts;

  // select one at random
    NodeLocationType pos1;
    GraphNeighborhoodIndexType	index=surfaceindexlist[randind];

    for (unsigned int t=0; t<GraphDimension; t++)
      pos1[t]=index[t];

    std::cout << " index " << index << std::endl;

    typename SearchNode::Pointer G1 = SearchNode::New();
    G1->SetLocation(pos1);
    this->SetSource(G1);	
    this->InitializeQueue();
    this->SetRadius(1);  
   
    m_MeanGeodesic=0;
    this->FindPath(); 

    float meandist=m_MeanGeodesic/(double)m_QS->GetTimer();

    if (m_QS->GetTimer() > (0.8*numpts) )
    {
      meandistancelist[randind]=meandist;      
      aggmeandist+=meandist;
      mct++;
      if (meandist < minmeandist) minmeandist=meandist;
    }
    std::cout << " mean dist " << meandist << " at try " << mct << " of " << trynum << std::endl;
    this->EmptyQ();

    float avgmd=0.0;
    if (mct > 0) avgmd = aggmeandist/(float)mct;
    std::cout << " aggregate mean dist " << avgmd 
    << " compactness " << avgmd/m_SphereGeodesic << " at try " << ii << std::endl;    
    
    m_GraphCostImage->SetPixel(index,m_SurfaceLabel);
   }

   aggmeandist/=(float)mct;
   std::cout << " sphere G " << sphereG << 
    " our G " << aggmeandist << std::endl; 

   m_MeanGeodesic = aggmeandist;

   for (ii=0; ii<numpts; ii++)
   {
     float c=meandistancelist[ii];
     if ( c == 0.0) c = minmeandist;
     typename SearchNode::Pointer sn= m_Graph->GetPixel(surfaceindexlist[ii]);
     if (sn) sn->SetTotalCost(c);
   }

}


template<class TGraphSearchNode, typename TCostImageType > 
void 
GeodesicPaths<TGraphSearchNode,TCostImageType>
::GenerateEquivalentSphere(unsigned int numPts, float rad) 
{
  typedef GraphCostImageType sphereType;
  typedef typename sphereType::RegionType  sphereRegionType;
  typename sphereType::Pointer     sphere = sphereType::New();
  sphereRegionType        sphereRegion;

  typename sphereType::SizeType size;

  unsigned int j=0;
  for (j=0; j< GraphDimension; j++) 
    size[j]=(unsigned long) (2.*(rad+10.));
 
  sphereRegion.SetSize( size );
  sphere->SetLargestPossibleRegion( sphereRegion ); 
  sphere->SetRequestedRegion( sphereRegion );  
  sphere->SetBufferedRegion( sphereRegion );
  sphere->Allocate(); 




  unsigned int surfcount=0;

  while (surfcount < numPts)
  {
    surfcount=0;
    typedef typename sphereType::IndexType IndexType;
    IndexType cind;
  
    for (j=0; j<GraphDimension; j++) cind[j]=size[j]/2;

    ImageRegionIteratorWithIndex<sphereType> 
      ti( sphere, sphere->GetLargestPossibleRegion() );

    ti.GoToBegin(); 
    while(!ti.IsAtEnd()  )
    {
      IndexType index=ti.GetIndex();

      float d=0;
      for (unsigned int j=0; j<GraphDimension; j++)
        d+=(cind[j]-index[j])*(cind[j]-index[j]);

      d=sqrt(d);
  
      if (fabs(d-rad) <= 0.5) 
      {
        ti.Set( (GraphCostImagePixelType) m_SurfaceLabel);
        surfcount++; 
      }
      else ti.Set(0);

    	++ti;
    }
    std::cout << " surf count " << surfcount << " goal " << numPts << std::endl;
    rad+=1.0;
  } 

  typename GraphCostImageType::Pointer tempimage=m_GraphCostImage;

  m_GraphCostImage=sphere;

  this->MeanGeodesic(1,0);
  m_SphereGeodesic = m_MeanGeodesic;

  m_GraphCostImage = tempimage; 
  m_GraphSize=m_GraphCostImage->GetLargestPossibleRegion().GetSize();

}


template<class TGraphSearchNode, typename TCostImageType > 
void 
GeodesicPaths<TGraphSearchNode,TCostImageType>::writeimage(const char* fn) 
{

  std::cout << " writing image " << std::endl;
  FILE *fbin; 
  fbin=fopen(fn,"wb"); 
	GraphIteratorType Iter( m_Graph, m_GraphRegion );
  Iter.GoToBegin();
  
  for( ; !Iter.IsAtEnd(); ++Iter  )
  {
    // write binary image
    m_GraphIndex = Iter.GetIndex(); 
    m_CurrentNode = Iter.Get();    float tf=0;
    if (m_CurrentNode)
    {
      tf = (float ) m_CurrentNode->GetTotalCost();
      if (fabs(tf) > 1.e5) tf=0;
    }
    fwrite(&tf,sizeof(tf),1,fbin); 
  }
  fclose(fbin);	
}

} // end namespace itk

#endif
