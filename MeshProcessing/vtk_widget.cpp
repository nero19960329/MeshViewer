#include "vtk_widget.h"

#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>

VTKWidget::VTKWidget(QWidget * parent) :
	QVTKWidget(parent) {
	this->renderer = vtkSmartPointer<vtkRenderer>::New();
	this->renderer->GradientBackgroundOn();
	this->renderer->SetBackground2(.1, .1, .4);
	this->renderer->SetBackground(.9, .9, .9);
	this->GetRenderWindow()->AddRenderer(this->renderer);
	this->update();
}

VTKWidget::~VTKWidget() {}

vtkSmartPointer<vtkActor> VTKWidget::addActor(vtkSmartPointer<vtkPolyData> mesh) {
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(mesh);

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->GetProperty()->SetAmbient(0.5);
	actor->GetProperty()->SetDiffuse(0.6);
	actor->GetProperty()->SetSpecular(0.3);
	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(.7, .7, .7);

	this->renderer->AddActor(actor);
	this->update();

	return actor;
}

void VTKWidget::resetCamera() {
	this->renderer->ResetCamera();
}