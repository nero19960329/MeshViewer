#include "mesh_processing.h"

#include <QtCore/QStringList>
#include <QtWidgets/QFileDialog>

#include <vtkAppendPolyData.h>
#include <vtkExtractEdges.h>
#include <vtkMath.h>
#include <vtkLineSource.h>
#include <vtkOBJReader.h>
#include <vtkPoints.h>
#include <vtkProperty.h>
#include <vtkSphereSource.h>

#include "mesh_operation.h"
#include "vtkOFFReader.h"

MeshProcessing::MeshProcessing(QWidget *parent) : QMainWindow(parent) {
	ui.setupUi(this);

	this->mesh_processing_data_model_ = MeshProcessingDataModel::getInstance();
	this->vtk_widget_ = ui.vtk_widget;
	this->open_file_action_ = ui.open_file_action;
	this->wireframe_mode_action_ = ui.wireframe_mode_action;
	this->observe_mode_action_ = ui.observe_mode_action;
	this->vertex_mode_action_ = ui.vertex_mode_action;
	this->face_mode_action_ = ui.face_mode_action;
	this->display_normal_action_ = ui.display_normal_action;
	this->list_widget_model_ = ui.list_widget_model;

	this->wireframe_mode_action_->setChecked(true);

	connect(this->open_file_action_, SIGNAL(triggered()), this, SLOT(OnOpenFile()));
	connect(this->wireframe_mode_action_, SIGNAL(triggered()), this, SLOT(OnWireframeMode()));
	connect(this->observe_mode_action_, SIGNAL(triggered()), this, SLOT(OnObserveMode()));
	connect(this->vertex_mode_action_, SIGNAL(triggered()), this, SLOT(OnVertexMode()));
	connect(this->face_mode_action_, SIGNAL(triggered()), this, SLOT(OnFaceMode()));
	connect(this->display_normal_action_, SIGNAL(triggered()), this, SLOT(OnDisplayNormal()));
	connect(this->list_widget_model_, SIGNAL(itemChanged(QListWidgetItem *)), this, SLOT(OnListWidgetModelItemChanged(QListWidgetItem *)));

	connect(this->vtk_widget_->style, SIGNAL(selectVertex(vtkIdType)), this, SLOT(OnSelectVertex(vtkIdType)));
	connect(this->vtk_widget_->style, SIGNAL(selectFace(vtkIdType)), this, SLOT(OnSelectFace(vtkIdType)));
}

void MeshProcessing::OnOpenFile() {
	QStringList file_paths = QFileDialog::getOpenFileNames(this, QString::fromLocal8Bit("打开文件"), "../", tr("OBJ files(*.obj);; OFF files(*.off)"));
	if (file_paths.size() == 0) return;

	this->removeMeshActors();
	this->removeVertexActors();
	this->removeFaceActors();
	this->mesh_processing_data_model_->mesh_vec_.clear();
	this->mesh_processing_data_model_->neighbor_vertex_actor_vec_.clear();
	this->mesh_processing_data_model_->actor_vec_.clear();
	this->mesh_processing_data_model_->highlight_vec_.clear();
	this->mesh_processing_data_model_->mean_edge_length_vec_.clear();
	this->list_widget_model_->clear();

	int number_of_points = 0, number_of_faces = 0, number_of_edges = 0;
	for (int i = 0; i < file_paths.size(); ++i) {
		QString file_name = file_paths[i];
		QString extension_name = file_name.split('/').back().split('.').back();

		vtkSmartPointer<vtkPolyData> mesh =
			vtkSmartPointer<vtkPolyData>::New();
		if (extension_name == QString::fromLocal8Bit("obj")) {
			vtkSmartPointer<vtkOBJReader> objReader =
				vtkSmartPointer<vtkOBJReader>::New();
			objReader->SetFileName(file_name.toLocal8Bit());
			objReader->Update();
			mesh->DeepCopy(objReader->GetOutput());
		} else {
			vtkSmartPointer<vtkOFFReader> offReader =
				vtkSmartPointer<vtkOFFReader>::New();
			offReader->SetFileName(file_name.toLocal8Bit());
			offReader->Update();
			mesh->DeepCopy(offReader->GetOutput());
		}

		QString only_name = file_name.split('/').back().split('.').front();
		QListWidgetItem * item = new QListWidgetItem(only_name);
		item->setCheckState(Qt::Checked);
		this->list_widget_model_->addItem(item);

		auto actor = this->vtk_widget_->addActor(mesh);
		this->vtk_widget_->highlightMesh(actor);

		auto wireframe_actor = this->vtk_widget_->addWireFrameActor(mesh);

		this->mesh_processing_data_model_->mesh_vec_.push_back(mesh);
		this->mesh_processing_data_model_->wireframe_actor_vec_.push_back(wireframe_actor);
		this->mesh_processing_data_model_->actor_vec_.push_back(actor);
		this->mesh_processing_data_model_->highlight_vec_.push_back(1);

		number_of_points += mesh->GetNumberOfPoints();
		number_of_faces += mesh->GetNumberOfCells();

		vtkSmartPointer<vtkExtractEdges> extractEdges =
			vtkSmartPointer<vtkExtractEdges>::New();
		extractEdges->SetInputData(mesh);
		extractEdges->Update();

		vtkSmartPointer<vtkPolyData> mesh_edge =
			vtkSmartPointer<vtkPolyData>::New();
		mesh_edge->DeepCopy(extractEdges->GetOutput());

		this->mesh_processing_data_model_->mesh_edge_vec_.push_back(mesh_edge);

		vtkSmartPointer<vtkIdList> pts =
			vtkSmartPointer<vtkIdList>::New();
		double mean_edge_length = 0.0;
		while (extractEdges->GetOutput()->GetLines()->GetNextCell(pts)) {
			double p0[3], p1[3];
			mesh_edge->GetPoint(pts->GetId(0), p0);
			mesh_edge->GetPoint(pts->GetId(1), p1);
			mean_edge_length += std::sqrt(vtkMath::Distance2BetweenPoints(p0, p1));
		}
		mean_edge_length /= mesh_edge->GetNumberOfLines();

		this->mesh_processing_data_model_->mean_edge_length_vec_.push_back(mean_edge_length);
		number_of_edges += mesh_edge->GetNumberOfLines();
	}

	this->resetParameters();

	this->vtk_widget_->updateBottomText(number_of_points, number_of_faces, number_of_edges);
	this->vtk_widget_->resetCamera();
}

void MeshProcessing::OnListWidgetModelItemChanged(QListWidgetItem * item) {
	int row_id = this->list_widget_model_->row(item);
	
	if (item->checkState() == Qt::Unchecked) {
		this->vtk_widget_->unhighlightMesh(this->mesh_processing_data_model_->actor_vec_[row_id]);
		this->mesh_processing_data_model_->highlight_vec_[row_id] = 0;
	} else {
		this->vtk_widget_->highlightMesh(this->mesh_processing_data_model_->actor_vec_[row_id]);
		this->mesh_processing_data_model_->highlight_vec_[row_id] = 1;
	}

	this->resetParameters();

	this->vtk_widget_->update();
}

void MeshProcessing::OnWireframeMode() {
	bool isChecked = this->wireframe_mode_action_->isChecked();
	for (auto actor : this->mesh_processing_data_model_->wireframe_actor_vec_)
		actor->SetVisibility(isChecked);
	this->vtk_widget_->update();
}

void MeshProcessing::OnObserveMode() {
	this->removeVertexActors();
	this->removeFaceActors();
	this->mesh_processing_data_model_->pick_mode_ = MeshProcessingDataModel::OBSERVE;
	this->vtk_widget_->updateTopText();
}

void MeshProcessing::OnVertexMode() {
	this->removeFaceActors();
	this->mesh_processing_data_model_->pick_mode_ = MeshProcessingDataModel::VERTEX;
	this->vtk_widget_->updateTopText();
}

void MeshProcessing::OnFaceMode() {
	this->removeVertexActors();
	this->mesh_processing_data_model_->pick_mode_ = MeshProcessingDataModel::FACE;
	this->vtk_widget_->updateTopText();
}

void MeshProcessing::OnDisplayNormal() {
	if (this->mesh_processing_data_model_->pick_mode_ != MeshProcessingDataModel::FACE ||
		this->mesh_processing_data_model_->selected_face_id_ == -1)
		return;

	int cellId = this->mesh_processing_data_model_->selected_face_id_;
	vtkCell * cell = this->mesh_processing_data_model_->combined_mesh_->GetCell(cellId);

	double p0[3], p1[3], p2[3];
	this->mesh_processing_data_model_->combined_mesh_->GetPoint(cell->GetPointId(0), p0);
	this->mesh_processing_data_model_->combined_mesh_->GetPoint(cell->GetPointId(1), p1);
	this->mesh_processing_data_model_->combined_mesh_->GetPoint(cell->GetPointId(2), p2);

	double c[3];
	for (int i = 0; i < 3; ++i) c[i] = (p0[i] + p1[i] + p2[i]) / 3;

	double p0p1[3], p0p2[3];
	vtkMath::Subtract(p1, p0, p0p1);
	vtkMath::Subtract(p2, p0, p0p2);

	double normal[3];
	vtkMath::Cross(p0p1, p0p2, normal);
	vtkMath::Normalize(normal);

	double constant = (
		std::sqrt(vtkMath::Distance2BetweenPoints(p0, p1)) +
		std::sqrt(vtkMath::Distance2BetweenPoints(p1, p2)) +
		std::sqrt(vtkMath::Distance2BetweenPoints(p2, p0))
	);
	double e[3];
	for (int i = 0; i < 3; ++i) e[i] = c[i] + constant * normal[i];

	vtkSmartPointer<vtkLineSource> lineSource =
		vtkSmartPointer<vtkLineSource>::New();
	lineSource->SetPoint1(c);
	lineSource->SetPoint2(e);
	lineSource->Update();

	this->mesh_processing_data_model_->selected_face_normal_actor_ = this->vtk_widget_->addActor(lineSource->GetOutput());
	this->mesh_processing_data_model_->selected_face_normal_actor_->GetProperty()->SetColor(.8, .8, .2);
	this->mesh_processing_data_model_->selected_face_normal_actor_->GetProperty()->SetLineWidth(5);

	this->vtk_widget_->update();
}

void MeshProcessing::OnSelectVertex(vtkIdType id) {
	this->removeVertexActors();

	this->mesh_processing_data_model_->selected_vertex_actor_ = this->vtk_widget_->highlightVertex(this->mesh_processing_data_model_->combined_mesh_, id);
	this->mesh_processing_data_model_->selected_vertex_actor_->GetProperty()->SetInterpolationToGouraud();
	this->mesh_processing_data_model_->selected_vertex_actor_->GetProperty()->SetColor(.8, .2, .2);

	auto neighborVertexIds = MeshOperation::getConnectedVertices(this->mesh_processing_data_model_->combined_mesh_, id);
	for (const auto & neighbor : neighborVertexIds) {
		this->mesh_processing_data_model_->neighbor_vertex_actor_vec_.push_back(this->vtk_widget_->highlightVertex(this->mesh_processing_data_model_->combined_mesh_, neighbor));
		this->mesh_processing_data_model_->neighbor_vertex_actor_vec_.back()->GetProperty()->SetInterpolationToGouraud();
		this->mesh_processing_data_model_->neighbor_vertex_actor_vec_.back()->GetProperty()->SetColor(.8, .8, .2);
	}

	auto neighborFaceIds = MeshOperation::getVertexConnectedFaces(this->mesh_processing_data_model_->combined_mesh_, id);
	for (const auto & neighbor : neighborFaceIds) {
		this->mesh_processing_data_model_->neighbor_face_actor_vec_.push_back(this->vtk_widget_->highlightFace(this->mesh_processing_data_model_->combined_mesh_, neighbor));
		this->mesh_processing_data_model_->neighbor_face_actor_vec_.back()->GetProperty()->SetColor(.8, .5, .2);
	}

	this->vtk_widget_->update();
}

void MeshProcessing::OnSelectFace(vtkIdType id) {
	this->removeFaceActors();

	this->mesh_processing_data_model_->selected_face_actor_ = this->vtk_widget_->highlightFace(this->mesh_processing_data_model_->combined_mesh_, id);
	this->mesh_processing_data_model_->selected_face_actor_->GetProperty()->SetColor(.8, .2, .2);

	auto neighborFaceIds = MeshOperation::getFaceConnectedFaces(this->mesh_processing_data_model_->combined_mesh_, id);
	for (const auto & neighbor : neighborFaceIds) {
		this->mesh_processing_data_model_->neighbor_face2_actor_vec_.push_back(this->vtk_widget_->highlightFace(this->mesh_processing_data_model_->combined_mesh_, neighbor));
		this->mesh_processing_data_model_->neighbor_face2_actor_vec_.back()->GetProperty()->SetColor(.8, .8, .2);
	}

	this->mesh_processing_data_model_->selected_face_id_ = id;

	this->vtk_widget_->update();
}

void MeshProcessing::resetParameters() {
	vtkSmartPointer<vtkAppendPolyData> appendFilter =
		vtkSmartPointer<vtkAppendPolyData>::New();
	for (int i = 0; i < this->mesh_processing_data_model_->highlight_vec_.size(); ++i) {
		if (this->mesh_processing_data_model_->highlight_vec_[i])
			appendFilter->AddInputData(this->mesh_processing_data_model_->mesh_vec_[i]);
	}
	appendFilter->Update();
	this->mesh_processing_data_model_->combined_mesh_->DeepCopy(appendFilter->GetOutput());

	int edge_count = 0;
	this->mesh_processing_data_model_->mean_edge_length = 0.0;
	for (int i = 0; i < this->mesh_processing_data_model_->highlight_vec_.size(); ++i) {
		if (this->mesh_processing_data_model_->highlight_vec_[i]) {
			this->mesh_processing_data_model_->mean_edge_length += 
				this->mesh_processing_data_model_->mesh_edge_vec_[i]->GetNumberOfLines() * this->mesh_processing_data_model_->mean_edge_length_vec_[i];
			edge_count += this->mesh_processing_data_model_->mesh_edge_vec_[i]->GetNumberOfLines();
		}
	}
	this->mesh_processing_data_model_->mean_edge_length /= edge_count;
}

void MeshProcessing::removeMeshActors() {
	for (const auto & actor : this->mesh_processing_data_model_->actor_vec_)
		this->vtk_widget_->removeActor(actor);
	for (const auto & actor : this->mesh_processing_data_model_->wireframe_actor_vec_)
		this->vtk_widget_->removeActor(actor);
}

void MeshProcessing::removeVertexActors() {
	this->vtk_widget_->removeActor(this->mesh_processing_data_model_->selected_vertex_actor_);
	for (const auto & actor : this->mesh_processing_data_model_->neighbor_vertex_actor_vec_)
		this->vtk_widget_->removeActor(actor);
	for (const auto & actor : this->mesh_processing_data_model_->neighbor_face_actor_vec_)
		this->vtk_widget_->removeActor(actor);
}

void MeshProcessing::removeFaceActors() {
	this->vtk_widget_->removeActor(this->mesh_processing_data_model_->selected_face_actor_);
	for (const auto & actor : this->mesh_processing_data_model_->neighbor_face2_actor_vec_)
		this->vtk_widget_->removeActor(actor);
	this->vtk_widget_->removeActor(this->mesh_processing_data_model_->selected_face_normal_actor_);
}
