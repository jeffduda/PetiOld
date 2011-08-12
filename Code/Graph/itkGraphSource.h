/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGraphSource.h,v $
  Language:  C++
  Date:      $Date: 2008/11/11 11:51:34 $
  Version:   $Revision: 1.1 $

  Copyright ( c ) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGraphSource_h
#define __itkGraphSource_h

#if defined( _MSC_VER )
#pragma warning ( disable : 4786 )
#endif
#include "itkProcessObject.h"

namespace itk
{

/** \class GraphSource
 *  \brief Base class for all process objects that output graph data.
 *
 * GraphSource is the base class for all process objects that output
 * graph data. Specifically, this class defines the GetOutput() method
 * that returns a pointer to the output graph. The class also defines
 * some internal private data members that are used to manage streaming
 * of data.
 *
 * \ingroup DataSources
 */
template <class TOutputGraph>
class ITK_EXPORT GraphSource : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef GraphSource               Self;
  typedef ProcessObject             Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( GraphSource,ProcessObject );

  /** Some convenient typedefs. */
  typedef DataObject::Pointer               DataObjectPointer;
  typedef TOutputGraph                      OutputGraphType;
  typedef typename OutputGraphType::Pointer OutputGraphPointer;

  /** Get the graph output of this process object.  */
  OutputGraphType * GetOutput( void );
  OutputGraphType * GetOutput( unsigned int idx );

  /** Set the graph output of this process object. This call is slated
   * to be removed from ITK. You should GraftOutput() and possible
   * DataObject::DisconnectPipeline() to properly change the output. */
  void SetOutput( TOutputGraph *output );

  /** Graft the specified DataObject onto this ProcessObject's output.
   * This method grabs a handle to the specified DataObject's bulk
   * data to used as its output's own bulk data. It also copies the
   * region ivars and meta-data from the
   * specified data object into this filter's output data object. Most
   * importantly, however, it leaves the Source ivar untouched so the
   * original pipeline routing is intact. This method is used when a
   * process object is implemented using a mini-pipeline which is
   * defined in its GenerateData() method.  The usage is:
   *
   * \code
   *    // setup the mini-pipeline to process the input to this filter
   *    firstFilterInMiniPipeline->SetInput(  this->GetInput()  );
   *
   *    // setup the mini-pipeline to calculate the correct regions
   *    // and write to the appropriate bulk data block
   *    lastFilterInMiniPipeline->GraftOutput(  this->GetOutput()  );
   *
   *    // execute the mini-pipeline
   *    lastFilterInMiniPipeline->Update();
   *
   *    // graft the mini-pipeline output back onto this filter's output.
   *    // this is needed to get the appropriate regions passed back.
   *    this->GraftOutput(  lastFilterInMiniPipeline->GetOutput()  );
   * \endcode
   *
   * For proper pipeline execution, a filter using a mini-pipeline
   * must implement the GenerateInputRequestedRegion(),
   * GenerateOutputRequestedRegion(), GenerateOutputInformation() and
   * EnlargeOutputRequestedRegion() methods as necessary to reflect
   * how the mini-pipeline will execute ( in other words, the outer
   * filter's pipeline mechanism must be consistent with what the
   * mini-pipeline will do ). */
  virtual void GraftOutput( OutputGraphType *output );

  /** Make a DataObject of the correct type to used as the specified
   * output.  Every ProcessObject subclass must be able to create a
   * DataObject that can be used as a specified output. This method
   * is automatically called when DataObject::DisconnectPipeline() is
   * called.  DataObject::DisconnectPipeline, disconnects a data object
   * from being an output of its current source.  When the data object
   * is disconnected, the ProcessObject needs to construct a replacement
   * output data object so that the ProcessObject is in a valid state.
   * So DataObject::DisconnectPipeline eventually calls
   * ProcessObject::MakeOutput. Note that MakeOutput always returns a
   * SmartPointer to a DataObject. If a subclass of GraphSource has
   * multiple outputs of different types, then that class must provide
   * an implementation of MakeOutput(). */
  virtual DataObjectPointer MakeOutput( unsigned int idx );

protected:
  GraphSource();
  virtual ~GraphSource() {}
  void PrintSelf( std::ostream& os, Indent indent ) const;

  /** Requested region of Graph is specified as i of N unstructured regions.
   * Since all DataObjects should be able to set the requested region in
   * unstructured form, just copy output->RequestedRegion all inputs. */
  void GenerateInputRequestedRegion();

private:
  GraphSource( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  /** Used by streaming: The requested region of the output being processed
   * by the execute method. Set in the GenerateInputRequestedRegion method. */
  int m_GenerateDataRegion;
  int m_GenerateDataNumberOfRegions;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGraphSource.hxx"
#endif

#endif
