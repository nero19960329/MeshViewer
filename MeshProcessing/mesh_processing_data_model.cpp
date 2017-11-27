#include "mesh_processing_data_model.h"

MeshProcessingDataModel * MeshProcessingDataModel::instance_ = nullptr;

MeshProcessingDataModel::MeshProcessingDataModel() {
	this->combined_mesh_ = vtkSmartPointer<vtkPolyData>::New();
	this->selected_face_normal_actor_ = nullptr;
	this->selected_face_id_ = -1;
	this->pick_mode_ = OBSERVE;
	this->display_mode_ = DEFAULT;
}

MeshProcessingDataModel::~MeshProcessingDataModel() {}

MeshProcessingDataModel * MeshProcessingDataModel::getInstance() {
	if (instance_ == nullptr)
		instance_ = new MeshProcessingDataModel;
	return instance_;
}

void MeshProcessingDataModel::deleteInstance() {
	if (instance_ != nullptr) {
		delete instance_;
		instance_ = nullptr;
	}
}