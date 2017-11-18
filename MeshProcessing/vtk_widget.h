#pragma once

#include <vtkActor.h>
#include <vtkAutoInit.h>
#include <vtkPolyData.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

#include <QVTKWidget.h>

VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);

class VTKWidget : public QVTKWidget {
private:
	vtkSmartPointer<vtkRenderer> renderer;

public:
	VTKWidget(QWidget * parent = nullptr);
	~VTKWidget();

public:
	vtkSmartPointer<vtkActor> addActor(vtkSmartPointer<vtkPolyData> mesh);
	void resetCamera();
};