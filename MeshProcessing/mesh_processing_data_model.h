#pragma once

#include <unordered_set>
#include <vector>

#include <vtkActor.h>
#include <vtkLookupTable.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include "mesh_segmenter.h"

class MeshProcessingDataModel {
public:
	MeshProcessingDataModel();
	~MeshProcessingDataModel();

	static MeshProcessingDataModel * getInstance();
	static void deleteInstance();

private:
	static MeshProcessingDataModel * instance_;

public:
	enum PickMode { OBSERVE, VERTEX, FACE, MULTI_VERTEX };
	enum DisplayMode { DEFAULT, DISCRETE, CONTINUOUS };

	std::vector<vtkSmartPointer<vtkPolyData>> mesh_vec_;
	std::vector<vtkSmartPointer<vtkPolyData>> mesh_edge_vec_;
	std::vector<vtkSmartPointer<vtkActor>> actor_vec_;
	std::vector<vtkSmartPointer<vtkActor>> wireframe_actor_vec_;
	std::vector<int> highlight_vec_;
	std::vector<double> mean_edge_length_vec_;

	vtkSmartPointer<vtkPolyData> combined_mesh_;

	double mean_edge_length;
	vtkSmartPointer<vtkActor> selected_vertex_actor_;
	std::vector<vtkSmartPointer<vtkActor>> neighbor_vertex_actor_vec_;
	std::vector<vtkSmartPointer<vtkActor>> neighbor_face_actor_vec_;
	vtkSmartPointer<vtkActor> selected_face_actor_;
	std::vector<vtkSmartPointer<vtkActor>> neighbor_face2_actor_vec_;

	vtkIdType selected_face_id_;
	vtkSmartPointer<vtkActor> selected_face_normal_actor_;

	std::unordered_set<vtkIdType> selected_multi_vertex_set_;
	std::vector<vtkSmartPointer<vtkActor>> selected_multi_vertex_actor_vec_;

	std::vector<vtkSmartPointer<vtkActor>> fill_region_face_actor_vec_;

	int source_id, target_id;

	MeshSegmenter * mesh_segmenter_;
	int segment_id;

	vtkSmartPointer<vtkLookupTable> hueLut;

	PickMode pick_mode_;
	DisplayMode display_mode_;
};