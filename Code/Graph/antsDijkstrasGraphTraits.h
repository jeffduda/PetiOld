/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: antsDijkstrasGraphTraits.h,v $
  Language:  C++
  Date:
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __antsDijkstrasGraphTraits_h
#define __antsDijkstrasGraphTraits_h

#include "itkImageGraphTraits.h"

namespace itk {
namespace ants {

/**
 * Graph traits class for use with the DijkstraPathGraphFilter class.
 */

template <typename TWeight= float, unsigned int VImageDimension = 3>
class DijkstrasGraphTraits : public ImageGraphTraits<TWeight, VImageDimension>
{
public:
  typedef DijkstrasGraphTraits                          Self;
  typedef ImageGraphTraits<TWeight, VImageDimension> Superclass;

  typedef TWeight                                 NodeWeightType;
  typedef TWeight                                 EdgeWeightType;
  typedef typename Superclass::NodeIdentifierType NodeIdentifierType;
  typedef typename Superclass::EdgeIdentifierType EdgeIdentifierType;
  typedef typename Superclass::EdgeIdentifierContainerType
                                                  EdgeIdentifierContainerType;
  typedef typename Superclass::IndexType          IndexType;
  typedef typename Superclass::EdgeType           EdgeType;
  typedef typename Superclass::EdgePointerType    EdgePointerType;
  typedef typename std::vector<NodeIdentifierType>
                                                  NodeIdentifierContainerType;
  struct  NodeType;
  typedef NodeType* NodePointerType;
  
  itkStaticConstMacro(ImageDimension, unsigned int, VImageDimension);

  struct NodeType
    {
    NodeIdentifierType Identifier;
    EdgeIdentifierContainerType IncomingEdges;
    EdgeIdentifierContainerType OutgoingEdges;
    NodeIdentifierContainerType PreviousNodes;
    
    int TimeStamp;
    bool IsTarget;
    bool IsSource;
    bool Visited;
    
    NodeWeightType Weight;
    IndexType ImageIndex;
    };
    
  //static bool EdgeSortPredicate(const EdgeType& lhs, const EdgeType& rhs)
  //  {
  //  return lhs.Weight < rhs.Weight;
  //  }  
};

}
} // end namespace itk

#endif
