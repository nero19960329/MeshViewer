#pragma once

#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkSetGet.h>
#include <vtkUnsignedCharArray.h>
#include <vtkVector.h>

#include <QtWidgets/QWidget>

#include <vector>

#include "mesh_processing_data_model.h"

class MeshProcessingInteractorStyle : public QWidget, public vtkInteractorStyleTrackballCamera {
	Q_OBJECT

private:
	MeshProcessingDataModel * mesh_processing_data_model_;

	bool is_right_down;
	bool moving;
	int start_pos[2], end_pos[2];
	vtkUnsignedCharArray * pixel_array;

public:
	static MeshProcessingInteractorStyle * New();
	vtkTypeMacro(MeshProcessingInteractorStyle, vtkInteractorStyleTrackballCamera);

	void OnMouseMove() override;
	void OnRightButtonDown() override;
	void OnRightButtonUp() override;

protected:
	MeshProcessingInteractorStyle();

signals:
	void selectVertex(vtkIdType id);
	void selectFace(vtkIdType id);
	void selectMultiVertex(const std::vector<vtkIdType> & ids);

private:
	void pickVertex();
	void pickFace();
	void startPickingMultiVertex();
	void movePickingMultiVertex();
	void endPickingMultiVertex();

	void drawPolygon();
	void drawPixels(const vtkVector2i & startPos, const vtkVector2i & endPos, unsigned char * pixels, int * size);
};