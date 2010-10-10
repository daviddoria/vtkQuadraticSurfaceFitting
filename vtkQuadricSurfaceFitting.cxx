#include "vtkQuadricSurfaceFitting.h"
#include "vtkBestFitPlane.h"

#include "vtkObjectFactory.h"
#include "vtkLandmarkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkFieldData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkPolyData.h"
#include "vtkDoubleArray.h"
#include "vtkPlane.h"
#include "vtkImageData.h"
#include "vtkSmartPointer.h"
#include "vtkMath.h"
#include "vtkQuadric.h"
#include "vtkContourFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkSampleFunction.h"

vtkStandardNewMacro(vtkQuadricSurfaceFitting);

namespace //anonymous
{
  /* allocate memory for an nrow x ncol matrix */
  template<class TReal>
      TReal **create_matrix ( long nrow, long ncol )
  {
    typedef TReal* TRealPointer;
    TReal **m = new TRealPointer[nrow];

    TReal* block = ( TReal* ) calloc ( nrow*ncol, sizeof ( TReal ) );
    m[0] = block;
    for ( int row = 1; row < nrow; ++row )
    {
      m[ row ] = &block[ row * ncol ];
    }
    return m;
  }

  /* free a TReal matrix allocated with matrix() */
  template<class TReal>
      void free_matrix ( TReal **m )
  {
    free ( m[0] );
    delete[] m;
  }
} //end anonymous namespace

int vtkQuadricSurfaceFitting::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // Get the input and ouptut
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkPolyData *input = vtkPolyData::SafeDownCast(
      inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkPolyData *output = vtkPolyData::SafeDownCast(
		  outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Compute the best fit plane.
  vtkSmartPointer<vtkBestFitPlane> bestFitFilter =
    vtkSmartPointer<vtkBestFitPlane>::New();
  bestFitFilter->SetInputConnection(input->GetProducerPort());
  bestFitFilter->Update();

  double bestFitNormal[3];
  vtkDoubleArray* normalArray = vtkDoubleArray::SafeDownCast(bestFitFilter->
                      GetOutput()->GetFieldData()->GetArray("BestFitPlaneNormal"));
  normalArray->GetTupleValue(0,bestFitNormal);

  double bestFitOrigin[3];
  vtkDoubleArray* originArray = vtkDoubleArray::SafeDownCast(bestFitFilter->
                      GetOutput()->GetFieldData()->GetArray("BestFitPlaneOrigin"));
  originArray->GetTupleValue(0,bestFitOrigin);

  // Find a coordinate system on the plane
  double v1[3] = {bestFitNormal[0], bestFitNormal[1], bestFitNormal[2]};
  double v2[3];
  double v3[3];
  vtkMath::Perpendiculars(v1, v2, v3, 0);

  // Project points onto the plane
  vtkSmartPointer<vtkPoints> transformedPoints =
    vtkSmartPointer<vtkPoints>::New();

  for(vtkIdType i = 0; i < input->GetNumberOfPoints(); i++)
    {
    double p[3];
    input->GetPoint(i,p);
    double projected[3];
    vtkPlane::ProjectPoint(p, bestFitOrigin, bestFitNormal, projected);
    // Find the vector from the plane origin to the projected point in the new coordinate system
    // by projecting the vector from the plane origin to the projected point in the normal coordinate
    // system onto v2 and v2. This vector projected onto v1 is the height needed for the
    // least squares problem.
    double oldV[3];
    vtkMath::Subtract(projected, bestFitOrigin, oldV);
    double new1[3];
    double new2[3];
    double new3[3];
    vtkMath::ProjectVector(oldV, v1, new1);
    vtkMath::ProjectVector(oldV, v2, new2);
    vtkMath::ProjectVector(oldV, v3, new3);

    transformedPoints->InsertNextPoint(vtkMath::Norm(new2), vtkMath::Norm(new3), vtkMath::Norm(new1));
    }

  // Fit a quadric surface
  int numberOfSamples = input->GetNumberOfPoints();
  int numberOfVariables = 6;
  double **x = create_matrix<double> (numberOfSamples, numberOfVariables);

  double **y = create_matrix<double> ( numberOfSamples, 1 );

  for(unsigned int i = 0; i < numberOfSamples; i++)
  {
    double p[3];
    input->GetPoint(i,p);

    x[i][0] = pow(p[0],2); //x^2
    x[i][1] = pow(p[1],2); //y^2
    x[i][2] = p[0] * p[1]; //x*y
    x[i][3] = p[0]; // x
    x[i][4] = p[1]; //y
    x[i][5] = 1; // constant
    //std::cout << X[i][0] << " " << X[i][1] << " " << X[i][2] << " " << X[i][3] << " " << X[i][4] << " " << X[i][5] << std::endl;
    y[i][0] = p[2]; //z
  }

  double **m = create_matrix<double> ( numberOfVariables, 1 );

  vtkMath::SolveLeastSquares(numberOfSamples, x, numberOfVariables, y, 1, m);

  // Create the quadric function definition
  vtkSmartPointer<vtkQuadric> quadric =
    vtkSmartPointer<vtkQuadric>::New();
  //quadric->SetCoefficients(x^2,   y^2      0  xy       0  0  x        y        -1  z);
  quadric->SetCoefficients(m[0][0], m[1][0], 0, m[2][0], 0, 0, m[3][0], m[4][0], -1, m[5][0]);

  // Sample the quadric function
  vtkSmartPointer<vtkSampleFunction> sample =
    vtkSmartPointer<vtkSampleFunction>::New();
  sample->SetSampleDimensions(50,50,50);
  sample->SetImplicitFunction(quadric);

  double inputBounds[6];
  input->GetBounds(inputBounds);
  double border = .2;
  inputBounds[0] -= border;
  inputBounds[1] += border;
  inputBounds[2] -= border;
  inputBounds[3] += border;
  inputBounds[4] -= border;
  inputBounds[5] += border;
  sample->SetModelBounds(inputBounds);

  vtkSmartPointer<vtkContourFilter> contours =
    vtkSmartPointer<vtkContourFilter>::New();
  contours->SetInput(sample->GetOutput());
  contours->GenerateValues(1, 0.0, 0.0);
  contours->Update();

  free_matrix(x);
  free_matrix(m);


  double sourceOrigin[3];
  for(unsigned int i = 0; i < 3; i++)
    {
    sourceOrigin[i] = bestFitOrigin[i];
    }

  double sourceP1[3];
  vtkMath::Add(sourceOrigin, v1, sourceP1);

  double sourceP2[3];
  vtkMath::Add(sourceOrigin, v2, sourceP2);

  double sourceP3[3];
  vtkMath::Add(sourceOrigin, v3, sourceP3);

  vtkSmartPointer<vtkPoints> sourcePoints =
    vtkSmartPointer<vtkPoints>::New();
  sourcePoints->InsertNextPoint(sourceOrigin);
  sourcePoints->InsertNextPoint(sourceP1);
  sourcePoints->InsertNextPoint(sourceP2);
  sourcePoints->InsertNextPoint(sourceP3);

  vtkSmartPointer<vtkPoints> targetPoints =
    vtkSmartPointer<vtkPoints>::New();
  targetPoints->InsertNextPoint(0,0,0);
  targetPoints->InsertNextPoint(1,0,0);
  targetPoints->InsertNextPoint(0,1,0);
  targetPoints->InsertNextPoint(0,0,1);

  // Find the coordinates of the projected points in the new coordinate system
  vtkSmartPointer<vtkLandmarkTransform> landmarkTransform =
    vtkSmartPointer<vtkLandmarkTransform>::New();
  landmarkTransform->SetSourceLandmarks(sourcePoints);
  landmarkTransform->SetTargetLandmarks(targetPoints);
  landmarkTransform->SetModeToRigidBody();
  landmarkTransform->Update();

  vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter =
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  transformFilter->SetInputConnection(contours->GetOutputPort());
  transformFilter->SetTransform(landmarkTransform);

  output->ShallowCopy(transformFilter->GetOutput());

  return 1;
}
