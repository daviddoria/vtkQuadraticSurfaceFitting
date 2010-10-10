#include "vtkQuadricSurfaceFitting.h"

#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"

#include "vtkPolyData.h"
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

  std::cout << "Solution is: " << std::endl;
  for(unsigned int i = 0; i < numberOfVariables; i++)
    {
    std::cout << m[i][0] << " ";
    }

  std::cout << std::endl;

  // Create the quadric function definition
  vtkSmartPointer<vtkQuadric> quadric =
    vtkSmartPointer<vtkQuadric>::New();
  //quadric->SetCoefficients(x^2,   y^2      0  xy       0  0  x        y        -1  z);
  quadric->SetCoefficients(m[0][0], m[1][0], 0, m[2][0], 0, 0, m[3][0], m[4][0], -1, m[5][0]);

  // Sample the quadric function
  vtkSmartPointer<vtkSampleFunction> sample = vtkSmartPointer<vtkSampleFunction>::New();
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

  vtkSmartPointer<vtkContourFilter> contours = vtkSmartPointer<vtkContourFilter>::New();
  contours->SetInput(sample->GetOutput());
  contours->GenerateValues(1, 0.0, 0.0);
  contours->Update();

  free_matrix(x);
  free_matrix(m);

  output->ShallowCopy(contours->GetOutput());

  return 1;
}
