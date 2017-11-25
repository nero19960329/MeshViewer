#include "mesh_processing_interactor_style.h"

#include <vtkCellLocator.h>
#include <vtkCellPicker.h>
#include <vtkObjectFactory.h>
#include <vtkPointLocator.h>
#include <vtkProperty.h>
#include <vtkPropPicker.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>

vtkStandardNewMacro(MeshProcessingInteractorStyle);

void MeshProcessingInteractorStyle::OnRightButtonDown() {
	this->mesh_processing_data_model_ = MeshProcessingDataModel::getInstance();

	if (this->mesh_processing_data_model_->pick_mode_ == MeshProcessingDataModel::VERTEX)
		pickVertex();
	else if (this->mesh_processing_data_model_->pick_mode_ == MeshProcessingDataModel::FACE)
		pickFace();
}

void MeshProcessingInteractorStyle::pickVertex() {
	int * clickPos = this->GetInteractor()->GetEventPosition();

	vtkSmartPointer<vtkPropPicker> picker =
		vtkSmartPointer<vtkPropPicker>::New();
	picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());

	double * pos = picker->GetPickPosition();
	if (pos[0] == 0 && pos[1] == 0 && pos[2] == 0) return;

	vtkSmartPointer<vtkPointLocator> pointLocator =
		vtkSmartPointer<vtkPointLocator>::New();
	pointLocator->SetDataSet(this->mesh_processing_data_model_->combined_mesh_);
	pointLocator->BuildLocator();

	vtkIdType vertex_id = pointLocator->FindClosestPoint(pos);
	emit(selectVertex(vertex_id));
}

void MeshProcessingInteractorStyle::pickFace() {
	int * clickPos = this->GetInteractor()->GetEventPosition();

	vtkSmartPointer<vtkCellPicker> picker =
		vtkSmartPointer<vtkCellPicker>::New();
	picker->SetTolerance(0.0005);
	picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());
	if (picker->GetCellId() == -1) return;

	vtkIdType face_id = picker->GetCellId();
	emit(selectFace(face_id));
}
