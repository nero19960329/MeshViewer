#pragma once

#include <vtkDoubleArray.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

class VTKWidget;

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
	VTKWidget * vtk_widget_;

public:
	MeshSegmenter(int seed_num_, double phy_ratio_);

	void setMesh(vtkSmartPointer<vtkPolyData> mesh_) { mesh = mesh_; }

	void segment();

	vtkSmartPointer<vtkDoubleArray> getSegmentScalar(int n);

private:
	vtkSmartPointer<vtkMutableUndirectedGraph> calcDualGraph(double phy_ratio);

	vtkSmartPointer<vtkPolyData> removeAbnormalRegion(vtkSmartPointer<vtkPolyData> mesh);

	void selectSeedsRandomly();
	void selectSeedsByOctree();
	
	std::vector<std::vector<double>> calcMultiDijkstraTable(std::vector<int> face_ids);
	std::vector<double> calcDijkstraTableWithMinHeap(int face_id);
	std::vector<double> calcDijkstraTableWithFibHeap(int face_id);
	std::vector<double> calcSPFATable(int face_id);

	void mergeClusters();
};