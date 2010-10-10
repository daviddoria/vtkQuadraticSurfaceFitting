// .NAME vtkQuadricSurfaceFitting - Fit a quadric surface to a set of points
// .SECTION Description
// vtkQuadricSurfaceFitting fits a quadric surface to a set of points.

#ifndef __vtkQuadricSurfaceFitting_h
#define __vtkQuadricSurfaceFitting_h

#include "vtkPolyDataAlgorithm.h"

class vtkQuadricSurfaceFitting : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkQuadricSurfaceFitting,vtkPolyDataAlgorithm);
  static vtkQuadricSurfaceFitting *New();
  vtkQuadricSurfaceFitting(){}

protected:
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
  vtkQuadricSurfaceFitting(const vtkQuadricSurfaceFitting&);  // Not implemented.
  void operator=(const vtkQuadricSurfaceFitting&);  // Not implemented.

};

#endif
