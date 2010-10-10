#include <vtkSmartPointer.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkProperty.h>
#include <vtkSphereSource.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

#include "vtkQuadricSurfaceFitting.h"

//new method:
// 1) find the best fit plane
// 2) pick an orthogonal basis using the normal
// 3) generate a new set of points using the distance of the points from the plane as the z value and the (x,y) in the new coordinate system
// 4) fit the surface
// 5) transform the surface using the transformation between the world coordinate system and the new coordinate system

void GenerateTopOfSphere(vtkPolyData*);
void GenerateSideOfSphere(vtkPolyData*);

int main(int, char *[])
{
  vtkSmartPointer<vtkPolyData> polydata =
    vtkSmartPointer<vtkPolyData>::New();
  //GenerateTopOfSphere(polydata);
  GenerateSideOfSphere(polydata);

  vtkSmartPointer<vtkPolyDataMapper> pointsMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  pointsMapper->SetInputConnection(polydata->GetProducerPort());

  vtkSmartPointer<vtkActor> pointsActor =
    vtkSmartPointer<vtkActor>::New();
  pointsActor->SetMapper(pointsMapper);
  pointsActor->GetProperty()->SetPointSize(4);

  vtkSmartPointer<vtkQuadricSurfaceFitting> quadricFitting =
    vtkSmartPointer<vtkQuadricSurfaceFitting>::New();
  quadricFitting->SetInputConnection(polydata->GetProducerPort());
  quadricFitting->Update();

  vtkSmartPointer<vtkPolyDataMapper> quadraticMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  quadraticMapper->SetInputConnection(quadricFitting->GetOutputPort());

  vtkSmartPointer<vtkActor> quadraticActor =
    vtkSmartPointer<vtkActor>::New();
  quadraticActor->SetMapper(quadraticMapper);

  // Create a renderer, render window, and interactor
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  // Add the actor to the scene
  renderer->AddActor(pointsActor);
  renderer->AddActor(quadraticActor);
  renderer->SetBackground(.3, .6, .3); // Background color green

  // Render and interact
  renderWindow->Render();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}

void GenerateTopOfSphere(vtkPolyData* polydata)
{
  vtkSmartPointer<vtkSphereSource> sphereSource =
    vtkSmartPointer<vtkSphereSource>::New();
  sphereSource->Update();

  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();

  for(vtkIdType i = 0; i < sphereSource->GetOutput()->GetNumberOfPoints(); i++)
    {
    double p[3];
    sphereSource->GetOutput()->GetPoint(i,p);
    if(p[2] > 0)
      {
      points->InsertNextPoint(p);
      }
    }

  vtkSmartPointer<vtkPolyData> temp =
    vtkSmartPointer<vtkPolyData>::New();

  temp->SetPoints(points);

  vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  glyphFilter->SetInputConnection(temp->GetProducerPort());
  glyphFilter->Update();

  polydata->ShallowCopy(glyphFilter->GetOutput());
}


void GenerateSideOfSphere(vtkPolyData* polydata)
{
  vtkSmartPointer<vtkSphereSource> sphereSource =
    vtkSmartPointer<vtkSphereSource>::New();
  sphereSource->Update();

  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();

  for(vtkIdType i = 0; i < sphereSource->GetOutput()->GetNumberOfPoints(); i++)
    {
    double p[3];
    sphereSource->GetOutput()->GetPoint(i,p);
    if(p[1] > 0)
      {
      points->InsertNextPoint(p);
      }
    }

  vtkSmartPointer<vtkPolyData> temp =
    vtkSmartPointer<vtkPolyData>::New();

  temp->SetPoints(points);

  vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  glyphFilter->SetInputConnection(temp->GetProducerPort());
  glyphFilter->Update();

  polydata->ShallowCopy(glyphFilter->GetOutput());
}