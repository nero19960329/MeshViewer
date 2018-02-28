#pragma once

#include <vtkActor.h>
#include <vtkAutoInit.h>
#include <vtkPolyData.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkTextActor.h>

#include <QVTKWidget.h>

#include "icp_algorithm.h"
#include "mesh_processing_data_model.h"
#include "mesh_processing_interactor_style.h"

VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);

class VTKWidget : public QVTKWidget {
public:
	MeshProcessingDataModel * mesh_processing_data_model_;

	vtkSmartPointer<vtkRenderer> renderer;

	vtkSmartPointer<vtkTextActor> bottomTextActor;
	vtkSmartPointer<vtkTextActor> topTextActor;

public:
	vtkSmartPointer<MeshProcessingInteractorStyle> style;

public:
	VTKWidget(QWidget * parent = nullptr);
	~VTKWidget();

public:
	vtkSmartPointer<vtkActor> addActor(vtkSmartPointer<vtkPolyData> mesh);
	vtkSmartPointer<vtkActor> addWireFrameActor(vtkSmartPointer<vtkPolyData> mesh);
	void removeActor(vtkSmartPointer<vtkActor> actor);

	void highlightMesh(vtkSmartPointer<vtkActor> actor);
	void unhighlightMesh(vtkSmartPointer<vtkActor> actor);

	vtkSmartPointer<vtkActor> highlightVertex(vtkSmartPointer<vtkPolyData> mesh, vtkIdType id);
	vtkSmartPointer<vtkActor> highlightFace(vtkSmartPointer<vtkPolyData> mesh, vtkIdType id);

	vtkSmartPointer<vtkActor> addLine(double * p1, double * p2);

	void updateTopText();
	void updateBottomText(int number_of_points, int number_of_faces, int number_of_edges);

	void resetCamera();

private:
	void initTextActor();
};