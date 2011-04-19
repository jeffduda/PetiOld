/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMaskImageEdgesToGraphFunctor.txx,v $
  Language:  C++
  Date:      $Date: 2008/11/11 11:51:34 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMaskImageEdgesToGraphFunctor_txx
#define __itkMaskImageEdgesToGraphFunctor_txx

#include "itkMaskImageEdgesToGraphFunctor.h"
#include "itkNeighborhood.h"

namespace itk
{

template<typename TInputImage, typename TOutputGraph>
bool 
MaskImageEdgesToGraphFunctor<TInputImage, TOutputGraph>
::IsPixelANode( IndexType idx )
{

  if ( this->GetInput( idx)->GetPixel(idx) != this->m_BackgroundValue )
    {
    itk::Neighborhood<typename InputImageType::PixelType, InputImageType::Dimension> neighborhood;
    }
  
  for (unsigned int i=0; i<InputImageType::Dimension; i++)
    {
     
   
  
  return false;
  
  


} // end namespace itk

}

#endif
