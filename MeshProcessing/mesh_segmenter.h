#pragma once

#include <vtkDoubleArray.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

class MeshSegmenter {
private:
	vtkSmartPointer<vtkPolyData> mesh;

	int seed_num;
	double phy_ratio;
	vtkSmartPointer<vtkMutableUndirectedGraph> dual_graph;
	vtkSmartPointer<vtkDoubleArray> centers, mesh_dis, edge_lens;
	std::vector<std::vector<int>> cluster_face_ids;
	std::vector<int> face_id_to_cluster;
	std::vector<std::vector<int>> cluster_steps;

public:
	MeshSegmenter();

	void setMesh(vtkSmartPointer<vtkPolyData> mesh_) { mesh = mesh_; }

	void segment();

	vtkSmartPointer<vtkDoubleArray> getSegmentScalar(int n);

private:
	vtkSmartPointer<vtkMutableUndirectedGraph> calcDualGraph(double phy_ratio);
	void randomSelectSeeds();
	std::vector<double> calcDijkstraTable(int face_id);
	void mergeClusters();
};