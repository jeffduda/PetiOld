/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit (ITK)
  Module:    $RCSfile: itkDijkstrasAlgorithm.h,v $ Language:  C++
  Date:      $Date: 2009/10/20 20:16:43 $
  Version:   $Revision: 1.1 $

Copyright (c) 2001 Insight Consortium
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

 * The name of the Insight Consortium, nor the names of any consortium members,
   nor of any contributors, may be used to endorse or promote products derived
   from this software without specific prior written permission.

  * Modified source versions must be plainly marked as such, and must not be
    misrepresented as being the original software.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
#ifndef _antsDijkstrasAlgorithm_h_
#define _antsDijkstrasAlgorithm_h_

#include <string>
#include <iostream>
#include <stack>
#include <vector>
#include <list>
#include <queue>
#include <map>

#include "vnl/vnl_math.h"
//#include "vnl/vnl_matrix_fixed.h"
//#include "vnl/vnl_vector_fixed.h"
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNeighborhoodIterator.h"
#include "itkVector.h"
#include "antsGraphSearchNode.h"

using namespace std;

namespace itk {
namespace ants {


// Forward declaration of DijkstrasAlgorithm so it can be declared a friend
template<class TGraphSearchNode> class DijkstrasAlgorithm;

template<class TGraphSearchNode>
class DijkstrasAlgorithmQueue : public itk::LightObject
/** \class DijkstrasAlgorithmQueue
the class containing the priority queue and associated data. 
*/
{    
private:
  
  template<class GraphSearchNode>
  class GraphSearchNodePriority /* defines the comparison operator for the prioritiy queue */
  {
    public:  
    bool operator() ( typename GraphSearchNode::Pointer Node1, 
                      typename GraphSearchNode::Pointer Node2) 
	  {  
      return Node1->GetTotalCost() > Node2->GetTotalCost();
	  }
  };

public: /* Standard typedefs.*/
  typedef DijkstrasAlgorithmQueue Self;                   
  typedef LightObject Superclass;         
  typedef SmartPointer<Self> Pointer;   
  typedef SmartPointer<const Self> ConstPointer;  
  itkTypeMacro(DijkstrasAlgorithmQueue,LightObject);
  itkNewMacro(Self);  /** Method for creation through the object factory.   */

  typedef typename TGraphSearchNode::Pointer TGraphSearchNodePointer;
  typedef typename TGraphSearchNode::PixelType PixelType; /** pixel type for the cost */
  typedef typename TGraphSearchNode::CoordRep CoordRep;   /** type for coordinates */
  typedef typename std::priority_queue< typename TGraphSearchNode::Pointer,std::vector< typename TGraphSearchNode::Pointer>, 
	  GraphSearchNodePriority<TGraphSearchNode> >QType; /** the queue we are using */
  typedef vector < typename TGraphSearchNode::Pointer >  NodeListType;
  
  inline QType GetQ() { return m_Q; }

  void   AddToPath(TGraphSearchNodePointer G) { this->m_Path.push_back(G); }

  inline NodeListType GetPath() { return m_Path; }

  void   EmptyPath() { m_Path.clear(); }

  inline NodeListType GetSourceNodes() { return m_SourceNodes; }

  inline NodeListType GetSinkNodes() { return m_SinkNodes; }

  inline void IncrementTimer() { m_timer++;}

  inline long GetTimer() { return m_timer; }

  inline void EmptyQ() { while (!m_Q.empty()) m_Q.pop(); m_timer=0; m_SourceNodes.clear(); m_SinkNodes.clear();  }
  
  NodeListType m_SinkNodes; 
  NodeListType m_SourceNodes;       
  QType m_Q;
  NodeListType m_Path; 

protected:  
  friend class DijkstrasAlgorithm<TGraphSearchNode>; // so it can access this data easily
  
  DijkstrasAlgorithmQueue() { m_timer=0; }
  ~DijkstrasAlgorithmQueue() {}
  
private: 
  unsigned long m_timer;
  DijkstrasAlgorithmQueue(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented  
};


/**
 * \class DijkstrasAlgorithm
 * \brief General shortest path / greedy dynamic programming solver.
 *
 *  This class uses Dijkstra's algorithm to solve the shortest path problem.
 *  We use the stl priority queue which is not optimal for this problem, but which
 *  works o.k.  It's implemented to be used for general regularly connected
 *  graphs with fixed costs, or for the variational solution of an integral 
 *  curve matching energy.
 *  Note: we assume all edge weights are positive.  
 *  The class is implemented as a an abstract base class, with virtual functions for 
 *  LocalCost, SearchEdgeSet and FindPath.  LocalCost must be implemented by derived classes. 
 *  The connectivity of the graph defined
 *  here is always regular and is controlled by a set of neighborhood indices.
 *  the default is a radius 1 neighborhood with all entries used.  However, the
 *  user may also supply her own regular connectivity by selecting the size of
 *  the neighborhood and a subset of the indices which define the edges.  If
 *  the GraphSearchNode contains its edges, they may be used with minor modification
 *  to the function SearchEdgeSet.
 *  Another improvement would be to make the LocalCost function a pointer 
 *  to a function which could be set.   
 */ 
template<class TGraphSearchNode>
class DijkstrasAlgorithm : public itk::LightObject
{
public: 
  typedef DijkstrasAlgorithm       Self;                   
  typedef LightObject              Superclass;         
  typedef SmartPointer<Self>       Pointer;   
  typedef SmartPointer<const Self> ConstPointer;  
  itkTypeMacro(DijkstrasAlgorithm,LightObject);
  itkNewMacro(Self);  

  enum{GraphDimension=TGraphSearchNode::GraphDimension};  /** dimension of the graph */

// Computation Data
  typedef TGraphSearchNode                                                  SearchNode; /** dimension of the graph */
  typedef typename SearchNode::Pointer                                      SearchNodePointer;
  typedef typename SearchNode::PixelType                                    PixelType; /**  pixel type for the cost */
  typedef typename SearchNode::CoordRep                                     CoordRep;  /** coordinate type */
  typedef Image<SearchNodePointer,GraphDimension>                           GraphType;
  typedef typename GraphType::SizeType                                      GraphSizeType;
  typedef ImageRegionIteratorWithIndex<GraphType>                           GraphIteratorType; 
  typedef typename GraphType::RegionType                                    GraphRegionType;
  typedef typename DijkstrasAlgorithmQueue<TGraphSearchNode>::Pointer       QType;
  typedef typename DijkstrasAlgorithmQueue<TGraphSearchNode>::NodeListType  NodeListType;
  typedef itk::NeighborhoodIterator<GraphType>                              GraphNeighborhoodIteratorType; 
  typedef typename GraphNeighborhoodIteratorType::IndexType                 GraphNeighborhoodIndexType;
  typedef typename GraphNeighborhoodIteratorType::RadiusType                RadiusType;
  typedef typename TGraphSearchNode::NodeLocationType                       NodeLocationType;
  typedef  typename GraphType::IndexType                                    IndexType;
  
  
// FUNCTIONS
  void InitializeGraph();  /** initializes all graph values appropriately */
  void InitializeQueue();  /** initializes all queue values appropriately 
								call AFTER source and sink are set*/
  void InitializeEdgeTemplate(); /** helper function initializes edge set appropriately */
  void InitializeEdgeTemplate(vector < unsigned int >,unsigned int); /** user supplied edge template */
  void SetGraphSize(typename GraphType::SizeType Sz); /** the rectangular size of the graph */

  inline void EmptyQ() { m_QS->EmptyQ(); this->m_TotalCost=0; }

  /* adds a source to the source set */
  void SetSource(typename TGraphSearchNode::Pointer G)
  {
    m_QS->m_SourceNodes.push_back(G); 
    for (int i=0; i<GraphDimension; i++)
      {
	m_GraphIndex[i]=(long int)(G->GetLocation()[i]+0.5);
//      std::cout << " mgi " << m_GraphIndex[i];
      }
      m_Graph->SetPixel(m_GraphIndex,G);
  };  

  typename TGraphSearchNode::Pointer GetGraphNode(   IndexType index)
  {
    //    std::cout << " get node "  << index << std::endl;
    return m_Graph->GetPixel(index);
  };  

  // adds a sink to the sink set 
  void SetSink(typename TGraphSearchNode::Pointer G)
  {
    m_QS->m_SinkNodes.push_back(G); 
  }

  // Backtracks from the given node to its source node; 
  void BackTrack(typename TGraphSearchNode::Pointer SinkNode)
  {

    m_QS->m_Path.clear();

    typename TGraphSearchNode::Pointer G=SinkNode;
    typename TGraphSearchNode::Pointer P=SinkNode->GetPredecessor();
    if (!P || !G) return;
    float highcost=G->GetValue();
    if (G->GetTotalCost() > P->GetValue())
      {
      P->SetAncestor(G);
      P->SetValue(G->GetTotalCost());
      highcost=G->GetTotalCost();
      }
   

    while(P && G != P)
    {
      m_QS->m_Path.push_back(G);
      G=P;
      P=G->GetPredecessor();
      if (P->GetValue() < highcost)
      {
        P->SetValue(highcost);
        P->SetAncestor(G);
      }
    }

    if (!P) cout << " null pred "; //else cout << " pred == self \n";
    return;
  }


  // Inverse of backtrack - from the given node to its sink node; 
  void ForwardTrack(typename TGraphSearchNode::Pointer SinkNode)
  {
    typename TGraphSearchNode::Pointer G=SinkNode;
    typename TGraphSearchNode::Pointer P=SinkNode->GetAncestor();
    while(P && G != P && G)
    {
      if (P->GetValue() > G->GetValue()) G->SetValue(P->GetValue());
      if (G->GetValue() > P->GetValue()) P->SetValue(G->GetValue());
      G=P;
      P=G->GetAncestor();
    }
    return;
  }


  virtual  bool TerminationCondition();  /** decides when the algorithm stops */

  virtual void SearchEdgeSet();  /** loops over the neighbors in the graph */
  
  void CheckNodeStatus();  /** checks if the node has been explored already, its cost, etc. */

  virtual PixelType LocalCost();      /* computes the local cost */
   /* alternatively, we could pass the function as a template parameter 
      or set a function pointer.  the latter method is used in dijkstrasegment. */

  virtual void FindPath();  /* runs the algorithm */

  inline unsigned int GetPathSize()
  {
	return m_QS->m_Path.size();
  }

  inline typename TGraphSearchNode::Pointer GetPathAtIndex(unsigned int i)
  {
	return m_QS->m_Path[i];
  }

  inline typename TGraphSearchNode::Pointer GetNeighborNode()
  {
    return m_NeighborNode;
  }
  inline typename TGraphSearchNode::Pointer GetCurrentNode()
  {
    return m_CurrentNode;
  }

  void SetMaxCost(PixelType m){ m_MaxCost=m;}
  double GetTotalCost() { return m_TotalCost; }
  void SetSearchFinished(bool m){ m_SearchFinished=m;} 
  /** sets the boolean that indicates if the algorithm is done */

protected:


  QType m_QS;

  /** defines neighborhood connectivity */
  vector < unsigned int >   m_EdgeTemplate;

  /** used by the neighborhood iterator */
  RadiusType                m_Radius; 

  /** holds the predecessor node */
  typename TGraphSearchNode::Pointer       m_PredecessorNode;  

  /** holds the current node */
  typename TGraphSearchNode::Pointer       m_CurrentNode;  

  /** holds the current neighbor node */
  typename TGraphSearchNode::Pointer       m_NeighborNode; 
  
  /** holds all the graph information */
  typename GraphType::Pointer        m_Graph;  

  GraphRegionType           m_GraphRegion;

  GraphSizeType             m_GraphSize;  /** rectangular size of graph */ 

  typename GraphType::IndexType      m_GraphIndex;

  bool                      m_SearchFinished;

  PixelType                 m_NewCost;

  PixelType                 m_CurrentCost;

  PixelType                 m_MaxCost;  // This creates an insurmountable barrier unless all costs are max
  
  double m_TotalCost;

  unsigned long             m_NumberSearched;
  DijkstrasAlgorithm(); 
  ~DijkstrasAlgorithm(){};
 
private:
  
 
  DijkstrasAlgorithm(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented  

};

}
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "antsDijkstrasAlgorithm.hxx"
#endif

#endif
