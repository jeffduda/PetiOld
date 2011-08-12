/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit (ITK)
  Module:    $RCSfile: antsGraphSearchNode.h,v $ Language:  C++
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
#ifndef _antsGraphSearchNode_h_
#define _antsGraphSearchNode_h_

#include <vector>

#include "itkVector.h"
//#include "antsDijkstrasAlgorithm.h"

namespace itk {
namespace ants {


/** The GraphSearchNode class defines a general shortest path graph node.
*   The algorithm requires  
*   the node to have a pointer to itself and entry for the cumulative cost.
*   We also define an index to its location and a couple of booleans
*   for keeping track of the graph node's state.
*   We assume the connectivity between nodes is defined externally.  
*   The class will also be useful for minimum spanning trees and
*   other graph search algorithms.   Connectivity is defined externally.
*   May be worthwhile to implement arbitrary connectivity e.g. for random graphs.
*   One way to do this is to include a list of pointers which define
*   the neighbors of this node, similar to how the predecessor is defined.
*/
template <class TPixelType, typename TCoordRep=unsigned int, 
          unsigned int NGraphDimension=2>
class GraphSearchNode : public itk::LightObject
{
public:

  /* Standard typedefs.*/
  typedef GraphSearchNode Self;                   
  typedef LightObject Superclass;         
  typedef SmartPointer<Self> Pointer;   
  typedef SmartPointer<const Self> ConstPointer;  
  itkTypeMacro(GraphSearchNode,LightObject);
  itkNewMacro(Self);  /** Method for creation through the object factory.   */
  
  enum StateType {UnVisitedState, VisitedState, DeliveredState, UnVisitableState };
  enum{GraphDimension=NGraphDimension};
  typedef TPixelType PixelType;   /** defines the cost data type */
  typedef TCoordRep  CoordRep;    /** defines the location data type */
  typedef itk::Vector<CoordRep,GraphDimension>  NodeLocationType; 
 
//  typedef typename itk::Image<float,GraphDimension>::IndexType  NodeLocationType; 

 typedef std::vector < Pointer >  NodeListType;
 
//  typedef itk::Image<CoordRep,GraphDimension>::IndexType  NodeLocationType;

  inline void SetLocation(NodeLocationType loc)
  {
  	m_Location=loc;
  }
  inline  NodeLocationType GetLocation()
  {
	  return m_Location;
  }

  inline void SetTotalCost(TPixelType cost)
  {
  	m_TotalCost=cost;
  }
  inline void SetValue(TPixelType cost,int which=0)
  {
	  if (which <= 0) m_Value1=cost;
	  if (which == 1) m_Value2=cost;
	  if (which == 2) m_Value3=cost;
	  if (which >= 3) m_Value4=cost;
  }
  inline TPixelType GetValue(int which=0)
  {
	  if (which <= 0) return m_Value1;
	  if (which == 1) return m_Value2;
	  if (which == 2) return m_Value3;
	  if (which >= 3) return m_Value4;
	  return m_Value1;
  }
  inline void SetUnVisited()
  {
	  m_State=UnVisitedState;
  }
  inline void SetUnVisitable()
  {
	  m_State=UnVisitableState;
  }
  inline void SetVisited()
  {
	  m_State=VisitedState;
  }
  inline void SetDelivered()
  {
	  m_State=DeliveredState;
  }

  inline TPixelType GetTotalCost()
  {
	  return m_TotalCost;
  }

  inline void SetPredecessor(Pointer address)
  {
	  m_PredecessorAddress=address;
  }
  inline Pointer GetPredecessor()
  {
	  return m_PredecessorAddress;
  }

  inline void SetAncestor(Pointer address)
  {
	  m_AncestorAddress=address;
  }
  inline Pointer GetAncestor()
  {
	  return m_AncestorAddress;
  }

  inline bool GetUnVisited()
  {
  	if (m_State==UnVisitedState)
	  return true;
	  else return false;
  }


  inline bool GetUnVisitable()
  {
  	if (m_State==UnVisitableState)
	  return true;
	  else return false;
  }

  inline bool GetVisited()
  {
  	if (m_State==VisitedState)
	  return true;
	  else return false;
  }

  inline bool GetDelivered()
  {
	  if (m_State==DeliveredState)
	  return true;
	  else return false;
  }

  inline void SetState(StateType S)
  {
	  m_State=S;
  }

  inline StateType GetState()
  {
	  return m_State;
  }

  inline void SetIdentity(unsigned int i) {  m_Identity=i; }

  inline unsigned int GetIdentity() { return m_Identity; }

  inline int GetNumberOfNeighbors() { return m_Neighbors.size(); }
  inline Pointer GetNeighbor(int i) { return m_Neighbors[i]; }

  void SetNeighborSize(int i){m_Neighbors.resize(i);}

  NodeListType m_Neighbors;
  unsigned short          m_NumberOfNeighbors;
  unsigned int            m_Identity;

protected:

  GraphSearchNode()
  {
	  m_TotalCost=0.0;
      m_Value1=0.0;
      m_Value2=0.0;
      m_Value3=0.0;
      m_Value4=0.0;
	  m_PredecessorAddress=NULL;
	  m_AncestorAddress=NULL;
	  m_State=UnVisitedState;
	  m_NumberOfNeighbors=0;
	  m_Identity=0;
  }

  ~GraphSearchNode(){}

private:
  TPixelType m_TotalCost; /** keeps track of the minimum accumulated cost. */  
  TPixelType m_Value1; /** keeps track of the minimum accumulated cost. */
  TPixelType m_Value2; /** keeps track of the minimum accumulated cost. */
  TPixelType m_Value3; /** keeps track of the minimum accumulated cost. */
  TPixelType m_Value4; /** keeps track of the minimum accumulated cost. */

  StateType m_State;
  NodeLocationType m_Location; /** provides the location in the graph. */
  Pointer m_PredecessorAddress;/** provides the best predecessor address */
  Pointer m_AncestorAddress;/** provides the best predecessor address */ 
  

  GraphSearchNode(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented  

};


// Forward declaration of DijkstrasAlgorithm so it can be declared a friend
// template<class TGraphSearchNode> class DijkstrasAlgorithm;

}
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
//#include "antsGraphSearchNode.hxx"
#endif

#endif
