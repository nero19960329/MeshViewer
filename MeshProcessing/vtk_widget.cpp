#include "vtk_widget.h"

#include <vtkDataSetMapper.h>
#include <vtkExtractSelection.h>
#include <vtkIdTypeArray.h>
#include <vtkLineSource.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkSphereSource.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkUnstructuredGrid.h>

#include <sstream>

VTKWidget::VTKWidget(QWidget * parent) : QVTKWidget(parent) {
	this->mesh_processing_data_model_ = MeshProcessingDataModel::getInstance();

	this->renderer = vtkSmartPointer<vtkRenderer>::New();
	this->renderer->GradientBackgroundOn();
	this->renderer->SetBackground2(.1, .1, .4);
	this->renderer->SetBackground(.9, .9, .9);
	this->GetRenderWindow()->AddRenderer(this->renderer);

	this->style = vtkSmartPointer<MeshProcessingInteractorStyle>::New();
	this->style->SetDefaultRenderer(this->renderer);
	this->GetInteractor()->SetInteractorStyle(this->style);

	this->initTextActor();
	this->update();
}

VTKWidget::~VTKWidget() {}

vtkSmartPointer<vtkActor> VTKWidget::addActor(vtkSmartPointer<vtkPolyData> mesh) {
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(mesh);

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->GetProperty()->SetAmbient(0.2);
	actor->GetProperty()->SetDiffuse(0.6);
	actor->GetProperty()->SetSpecular(0.3);
	actor->SetMapper(mapper);
	actor->GetProperty()->SetInterpolationToFlat();
	actor->GetProperty()->SetColor(0.7, 0.7, 0.7);

	this->renderer->AddActor(actor);
	this->update();

	return actor;
}

vtkSmartPointer<vtkActor> VTKWidget::addWireFrameActor(vtkSmartPointer<vtkPolyData> mesh) {
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(mesh);
	mapper->SetScalarVisibility(false);

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->GetProperty()->SetAmbient(1.0);
	actor->GetProperty()->SetDiffuse(0.0);
	actor->GetProperty()->SetSpecular(0.0);
	actor->SetMapper(mapper);
	actor->GetProperty()->SetRepresentationToWireframe();
	actor->GetProperty()->SetColor(0.0, 0.0, 0.0);
	actor->GetProperty()->SetLineWidth(2);
	actor->PickableOff();

	this->renderer->AddActor(actor);
	this->update();

	return actor;
}

void VTKWidget::removeActor(vtkSmartPointer<vtkActor> actor) {
	this->renderer->RemoveActor(actor);
}

void VTKWidget::highlightMesh(vtkSmartPointer<vtkActor> actor) {
	if (this->mesh_processing_data_model_->display_mode_ == MeshProcessingDataModel::DEFAULT) {
		actor->GetProperty()->SetColor(0.0, 0.5, 1.0);
		actor->GetMapper()->SetScalarVisibility(false);
	} else {
		actor->GetProperty()->SetColor(0.0, 0.0, 0.0);
		actor->GetMapper()->SetScalarVisibility(true);
	}
}

void VTKWidget::unhighlightMesh(vtkSmartPointer<vtkActor> actor) {
	actor->GetProperty()->SetColor(0.4, 0.4, 0.4);
	if (this->mesh_processing_data_model_->display_mode_ != MeshProcessingDataModel::DEFAULT)
		actor->GetMapper()->SetScalarVisibility(false);
}

vtkSmartPointer<vtkActor> VTKWidget::highlightVertex(vtkSmartPointer<vtkPolyData> mesh, vtkIdType id) {
	vtkSmartPointer<vtkSphereSource> sphereSource =
		vtkSmartPointer<vtkSphereSource>::New();
	sphereSource->SetCenter(this->mesh_processing_data_model_->combined_mesh_->GetPoint(id));
	sphereSource->SetRadius(.2 * this->mesh_processing_data_model_->mean_edge_length);
	//sphereSource->SetRadius(2.0);
	sphereSource->SetThetaResolution(40);
	sphereSource->SetPhiResolution(40);
	sphereSource->Update();

	return this->addActor(sphereSource->GetOutput());
}

vtkSmartPointer<vtkActor> VTKWidget::highlightFace(vtkSmartPointer<vtkPolyData> mesh, vtkIdType id) {
	vtkSmartPointer<vtkIdTypeArray> ids =
		vtkSmartPointer<vtkIdTypeArray>::New();
	ids->SetNumberOfComponents(1);
	ids->InsertNextValue(id);

	vtkSmartPointer<vtkSelectionNode> selectionNode =
		vtkSmartPointer<vtkSelectionNode>::New();
	selectionNode->SetFieldType(vtkSelectionNode::CELL);
	selectionNode->SetContentType(vtkSelectionNode::INDICES);
	selectionNode->SetSelectionList(ids);

	vtkSmartPointer<vtkSelection> selection =
		vtkSmartPointer<vtkSelection>::New();
	selection->AddNode(selectionNode);

	vtkSmartPointer<vtkExtractSelection> extractSelection =
		vtkSmartPointer<vtkExtractSelection>::New();
	extractSelection->SetInputData(0, mesh);
	extractSelection->SetInputData(1, selection);
	extractSelection->Update();

	vtkSmartPointer<vtkUnstructuredGrid> selected =
		vtkSmartPointer<vtkUnstructuredGrid>::New();
	selected->ShallowCopy(extractSelection->GetOutput());

	vtkSmartPointer<vtkDataSetMapper> mapper =
		vtkSmartPointer<vtkDataSetMapper>::New();
	mapper->SetInputData(selected);
	mapper->SetScalarVisibility(false);

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->GetProperty()->SetAmbient(0.5);
	actor->GetProperty()->SetDiffuse(0.6);
	actor->GetProperty()->SetSpecular(0.3);
	actor->SetMapper(mapper);
	actor->GetProperty()->SetInterpolationToFlat();
	actor->GetProperty()->SetColor(0.8, 0.5, 0.2);

	this->renderer->AddActor(actor);
	this->update();

	return actor;
}

vtkSmartPointer<vtkActor> VTKWidget::addLine(double * p1, double * p2) {
	vtkNew<vtkSphereSource> sphereSource1;
	sphereSource1->SetCenter(p1);
	sphereSource1->SetRadius(.2);
	sphereSource1->SetThetaResolution(40);
	sphereSource1->SetPhiResolution(40);
	sphereSource1->Update();
	auto actor = this->addActor(sphereSource1->GetOutput());
	actor->GetProperty()->SetColor(.8, .2, .2);

	vtkNew<vtkSphereSource> sphereSource2;
	sphereSource2->SetCenter(p2);
	sphereSource2->SetRadius(.2);
	sphereSource2->SetThetaResolution(40);
	sphereSource2->SetPhiResolution(40);
	sphereSource2->Update();
	this->addActor(sphereSource2->GetOutput());

	vtkNew<vtkLineSource> lineSource;
	lineSource->SetPoint1(p1);
	lineSource->SetPoint2(p2);
	lineSource->Update();
	this->addActor(lineSource->GetOutput());

	return nullptr;
}

void VTKWidget::updateTopText() {
	if (this->mesh_processing_data_model_->pick_mode_ == MeshProcessingDataModel::OBSERVE)
		this->topTextActor->SetInput("Pick mode: Observe");
	else if (this->mesh_processing_data_model_->pick_mode_ == MeshProcessingDataModel::VERTEX)
		this->topTextActor->SetInput("Pick mode: Vertex");
	else if (this->mesh_processing_data_model_->pick_mode_ == MeshProcessingDataModel::FACE)
		this->topTextActor->SetInput("Pick mode: Face");
	else if (this->mesh_processing_data_model_->pick_mode_ == MeshProcessingDataModel::MULTI_VERTEX)
		this->topTextActor->SetInput("Pick mode : MultiVertex");
	this->update();
}

void VTKWidget::updateBottomText(int number_of_points, int number_of_faces, int number_of_edges) {
	using std::ostringstream;
	ostringstream oss;
	oss << "Number of points: " << number_of_points << "\n";
	oss << "Number of faces: " << number_of_faces << "\n";
	oss << "Number of edges: " << number_of_edges;
	this->bottomTextActor->SetInput(oss.str().c_str());
	this->update();
}

void VTKWidget::resetCamera() {
	this->renderer->ResetCamera();
}

void VTKWidget::initTextActor() {
	this->bottomTextActor = vtkSmartPointer<vtkTextActor>::New();
	this->bottomTextActor->SetInput("Number of points: 0\nNumber of faces: 0\nNumber of edges: 0");
	this->bottomTextActor->GetTextProperty()->SetFontSize(16);
	this->bottomTextActor->GetTextProperty()->SetColor(.2, .2, .2);
	this->bottomTextActor->GetTextProperty()->SetFontFamilyToTimes();
	this->bottomTextActor->GetTextProperty()->BoldOn();
	this->renderer->AddActor2D(bottomTextActor);

	this->topTextActor = vtkSmartPointer<vtkTextActor>::New();
	this->topTextActor->SetInput("Pick mode: Observe");
	this->topTextActor->SetPosition(0, 680);
	this->topTextActor->GetTextProperty()->SetFontSize(16);
	this->topTextActor->GetTextProperty()->SetColor(.8, .2, .2);
	this->topTextActor->GetTextProperty()->SetFontFamilyToTimes();
	this->topTextActor->GetTextProperty()->BoldOn();
	this->renderer->AddActor2D(topTextActor);
}