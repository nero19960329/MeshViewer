#include "mesh_processing_interactor_style.h"

#include <vtkCellLocator.h>
#include <vtkCellPicker.h>
#include <vtkIdTypeArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <vtkProperty.h>
#include <vtkPropPicker.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSelectVisiblePoints.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkVectorOperators.h>

vtkStandardNewMacro(MeshProcessingInteractorStyle);

MeshProcessingInteractorStyle::MeshProcessingInteractorStyle() {
	this->start_pos[0] = this->start_pos[1] = 0;
	this->end_pos[0] = this->end_pos[1] = 0;
	this->moving = false;
	this->pixel_array = vtkUnsignedCharArray::New();
	this->is_right_down = false;
	this->mesh_processing_data_model_ = MeshProcessingDataModel::getInstance();
}

void MeshProcessingInteractorStyle::OnMouseMove() {
	if (this->mesh_processing_data_model_->pick_mode_ == MeshProcessingDataModel::MULTI_VERTEX && this->is_right_down)
		movePickingMultiVertex();
	else
		vtkInteractorStyleTrackballCamera::OnMouseMove();
}

void MeshProcessingInteractorStyle::OnRightButtonDown() {
	this->is_right_down = true;
	if (this->mesh_processing_data_model_->combined_mesh_ == nullptr) return;

	if (this->mesh_processing_data_model_->pick_mode_ == MeshProcessingDataModel::VERTEX)
		pickVertex();
	else if (this->mesh_processing_data_model_->pick_mode_ == MeshProcessingDataModel::FACE)
		pickFace();
	else if (this->mesh_processing_data_model_->pick_mode_ == MeshProcessingDataModel::MULTI_VERTEX)
		startPickingMultiVertex();
}

void MeshProcessingInteractorStyle::OnRightButtonUp() {
	this->is_right_down = false;
	if (this->mesh_processing_data_model_->pick_mode_ == MeshProcessingDataModel::MULTI_VERTEX)
		endPickingMultiVertex();
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

	vtkSmartPointer<vtkPropPicker> picker =
		vtkSmartPointer<vtkPropPicker>::New();
	picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());

	double * pos = picker->GetPickPosition();
	if (pos[0] == 0 && pos[1] == 0 && pos[2] == 0) return;

	vtkSmartPointer<vtkCellLocator> cellLocator =
		vtkSmartPointer<vtkCellLocator>::New();
	cellLocator->SetDataSet(this->mesh_processing_data_model_->combined_mesh_);
	cellLocator->BuildLocator();

	vtkIdType face_id = cellLocator->FindCell(pos);
	emit(selectFace(face_id));
}

void MeshProcessingInteractorStyle::startPickingMultiVertex() {
	if (!this->Interactor || this->mesh_processing_data_model_->combined_mesh_ == nullptr) return;

	this->moving = true;

	vtkRenderWindow * renWin = this->Interactor->GetRenderWindow();

	this->start_pos[0] = this->Interactor->GetEventPosition()[0];
	this->start_pos[1] = this->Interactor->GetEventPosition()[1];
	this->end_pos[0] = this->start_pos[0];
	this->end_pos[1] = this->start_pos[1];

	this->pixel_array->Initialize();
	this->pixel_array->SetNumberOfComponents(3);
	int * size = renWin->GetSize();
	this->pixel_array->SetNumberOfTuples(size[0] * size[1]);

	renWin->GetPixelData(0, 0, size[0] - 1, size[1] - 1, 1, this->pixel_array);
	this->InvokeEvent(vtkCommand::StartInteractionEvent);
}

void MeshProcessingInteractorStyle::movePickingMultiVertex() {
	if (!this->Interactor || !this->moving) return;

	this->end_pos[0] = this->Interactor->GetEventPosition()[0];
	this->end_pos[1] = this->Interactor->GetEventPosition()[1];
	int * size = this->Interactor->GetRenderWindow()->GetSize();
	if (this->end_pos[0] > size[0] - 1) this->end_pos[0] = size[0] - 1;
	if (this->end_pos[0] < 0) this->end_pos[0] = 0;
	if (this->end_pos[1] > size[1] - 1) this->end_pos[1] = size[1] - 1;
	if (this->end_pos[1] < 0) this->end_pos[1] = 0;

	this->drawPolygon();
}

void MeshProcessingInteractorStyle::endPickingMultiVertex() {
	if (!this->Interactor || !this->moving) return;

	int * size = this->Interactor->GetRenderWindow()->GetSize();
	unsigned char * pixels = this->pixel_array->GetPointer(0);
	this->Interactor->GetRenderWindow()->SetPixelData(0, 0, size[0] - 1, size[1] - 1, pixels, 0);
	this->Interactor->GetRenderWindow()->Frame();

	int x1 = start_pos[0], x2 = end_pos[0];
	int y1 = start_pos[1], y2 = end_pos[1];
	if (x1 > x2) std::swap(x1, x2);
	if (y1 > y2) std::swap(y1, y2);

	vtkSmartPointer<vtkPolyData> pts =
		vtkSmartPointer<vtkPolyData>::New();
	pts->SetPoints(this->mesh_processing_data_model_->combined_mesh_->GetPoints());
	pts->GetPointData()->AddArray(this->mesh_processing_data_model_->combined_mesh_->GetPointData()->GetArray("number"));

	vtkNew<vtkSelectVisiblePoints> selectFilter;
	selectFilter->SetInputData(pts);
	selectFilter->SetTolerance(1e-4);
	selectFilter->SelectionWindowOn();
	selectFilter->SetSelection(x1, x2, y1, y2);
	selectFilter->SelectInvisibleOff();
	selectFilter->SetRenderer(this->GetDefaultRenderer());
	selectFilter->Update();

	vtkIdTypeArray * numberScalarArray = vtkIdTypeArray::SafeDownCast(selectFilter->GetOutput()->GetPointData()->GetArray("number"));
	std::vector<vtkIdType> ids;
	for (int i = 0; i < selectFilter->GetOutput()->GetNumberOfPoints(); ++i)
		ids.push_back(numberScalarArray->GetValue(i));
	emit(selectMultiVertex(ids));

	this->moving = false;
	this->InvokeEvent(vtkCommand::SelectionChangedEvent);
	this->InvokeEvent(vtkCommand::EndInteractionEvent);
}

void MeshProcessingInteractorStyle::drawPolygon() {
	vtkNew<vtkUnsignedCharArray> tmpPixelArray;
	tmpPixelArray->DeepCopy(this->pixel_array);
	unsigned char * pixels = tmpPixelArray->GetPointer(0);
	int * size = this->Interactor->GetRenderWindow()->GetSize();

	vtkVector2i start(this->start_pos[0], this->start_pos[1]);
	vtkVector2i end(this->end_pos[0], this->end_pos[1]);
	this->drawPixels(start, end, pixels, size);

	this->Interactor->GetRenderWindow()->SetPixelData(0, 0, size[0] - 1, size[1] - 1, pixels, 0);
	this->Interactor->GetRenderWindow()->Frame();
}

void MeshProcessingInteractorStyle::drawPixels(const vtkVector2i & startPos, const vtkVector2i & endPos, unsigned char * pixels, int * size) {
	int x1 = startPos.GetX(), x2 = endPos.GetX();
	int y1 = startPos.GetY(), y2 = endPos.GetY();

	if (x1 > x2) std::swap(x1, x2);
	if (y1 > y2) std::swap(y1, y2);

	for (int x = x1; x <= x2; ++x) for (int k = 0; k < 3; ++k) {
		pixels[3 * (y1 * size[0] + x) + k] = 255 ^ pixels[3 * (y1 * size[0] + x) + k];
		pixels[3 * (y2 * size[0] + x) + k] = 255 ^ pixels[3 * (y2 * size[0] + x) + k];
	}

	for (int y = y1; y <= y2; ++y) for (int k = 0; k < 3; ++k) {
		pixels[3 * (y * size[0] + x1) + k] = 255 ^ pixels[3 * (y * size[0] + x1) + k];
		pixels[3 * (y * size[0] + x2) + k] = 255 ^ pixels[3 * (y * size[0] + x2) + k];
	}
}
