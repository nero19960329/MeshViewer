#include "mesh_processing_data_model.h"

MeshProcessingDataModel * MeshProcessingDataModel::instance_ = nullptr;

MeshProcessingDataModel::MeshProcessingDataModel() {
	this->combined_mesh_ = vtkSmartPointer<vtkPolyData>::New();
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