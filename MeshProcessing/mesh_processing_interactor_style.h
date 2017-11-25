#pragma once

#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkSetGet.h>

#include <QtWidgets/QWidget>

#include "mesh_processing_data_model.h"

class MeshProcessingInteractorStyle : public QWidget, public vtkInteractorStyleTrackballCamera{
	Q_OBJECT

private:
	MeshProcessingDataModel * mesh_processing_data_model_;

public:
	static MeshProcessingInteractorStyle * New();
	vtkTypeMacro(MeshProcessingInteractorStyle, vtkInteractorStyleTrackballCamera);

	virtual void OnRightButtonDown();

signals:
	void selectVertex(vtkIdType id);
	void selectFace(vtkIdType id);

private:
	void pickVertex();
	void pickFace();
};