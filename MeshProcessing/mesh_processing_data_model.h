#pragma once

#include <vector>

#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

class MeshProcessingDataModel {
public:
	MeshProcessingDataModel();
	~MeshProcessingDataModel();

	static MeshProcessingDataModel * getInstance();
	static void deleteInstance();

private:
	static MeshProcessingDataModel * instance_;

public:
	enum PickMode { OBSERVE, VERTEX, FACE };

	std::vector<vtkSmartPointer<vtkPolyData>> mesh_vec_;
	std::vector<vtkSmartPointer<vtkPolyData>> mesh_edge_vec_;
	std::vector<vtkSmartPointer<vtkActor>> actor_vec_;
	std::vector<int> highlight_vec_;
	std::vector<double> mean_edge_length_vec_;

	vtkSmartPointer<vtkPolyData> combined_mesh_;

	double mean_edge_length;
	vtkSmartPointer<vtkActor> selected_vertex_actor_;
	std::vector<vtkSmartPointer<vtkActor>> neighbor_vertex_actor_vec_;
	std::vector<vtkSmartPointer<vtkActor>> neighbor_face_actor_vec_;
	vtkSmartPointer<vtkActor> selected_face_actor_;
	std::vector<vtkSmartPointer<vtkActor>> neighbor_face2_actor_vec_;

	PickMode pick_mode_;
};