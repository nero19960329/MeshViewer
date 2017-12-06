#include "mesh_segmenter.h"
#include "min_heap.h"
#include "vec3.h"

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkEdgeListIterator.h>
#include <vtkIdList.h>
#include <vtkInEdgeIterator.h>
#include <vtkMath.h>
#include <vtkNew.h>
#include <vtkPolyDataNormals.h>
#include <vtkTriangle.h>

#include <array>
#include <set>
#include <tuple>
#include <unordered_map>

MeshSegmenter::MeshSegmenter() {
	using std::vector;

	this->phy_ratio = 0.03;
	this->seed_num = 64;
	this->cluster_face_ids = vector<vector<int>>(this->seed_num);
}

void MeshSegmenter::segment() {
	using std::vector;

	this->dual_graph = this->calcDualGraph(this->phy_ratio);
	this->centers = vtkDoubleArray::SafeDownCast(dual_graph->GetVertexData()->GetArray("Centers"));
	this->mesh_dis = vtkDoubleArray::SafeDownCast(dual_graph->GetEdgeData()->GetArray("Weights"));
	this->edge_lens = vtkDoubleArray::SafeDownCast(dual_graph->GetEdgeData()->GetArray("EdgeLens"));

	this->randomSelectSeeds();

	vector<vector<double>> dists;
	for (int i = 0; i < this->seed_num; ++i)
		dists.push_back(vector<double>());
	for (int i = 0; i < seed_num; ++i)
		dists[i] = this->calcDijkstraTable(this->cluster_face_ids[i].front());

	vector<vector<int>> min_dis_ids(this->seed_num);
	for (int i = 0; i < this->mesh->GetNumberOfCells(); ++i) {
		double min_dis = DBL_MAX;
		int min_dis_id;
		for (int j = 0; j < this->seed_num; ++j) {
			if (min_dis > dists[j][i]) {
				min_dis = dists[j][i];
				min_dis_id = j;
			}
		}

		if (min_dis != DBL_MAX)
			min_dis_ids[min_dis_id].push_back(i);
	}

	face_id_to_cluster = vector<int>(this->mesh->GetNumberOfCells(), -1);
	for (int i = 0; i < this->seed_num; ++i) {
		this->cluster_face_ids[i] = vector<int>();
		for (int j = 0; j < min_dis_ids[i].size(); ++j) {
			this->cluster_face_ids[i].push_back(min_dis_ids[i][j]);
			this->face_id_to_cluster[min_dis_ids[i][j]] = i;
		}
	}

	this->mergeClusters();
}

vtkSmartPointer<vtkDoubleArray> MeshSegmenter::getSegmentScalar(int n) {
	vtkSmartPointer<vtkDoubleArray> scalars =
		vtkSmartPointer<vtkDoubleArray>::New();
	scalars->SetNumberOfValues(this->mesh->GetNumberOfCells());
	for (int i = 0; i < this->mesh->GetNumberOfCells(); ++i) {
		int cluster_id = this->cluster_steps[this->seed_num - n][this->face_id_to_cluster[i]];
		scalars->SetValue(i, cluster_id * 1.0 / n);
	}

	return scalars;
}

vtkSmartPointer<vtkMutableUndirectedGraph> MeshSegmenter::calcDualGraph(double phy_ratio) {
	using std::vector;
	using std::sqrt;

	vtkSmartPointer<vtkMutableUndirectedGraph> dual_graph =
		vtkSmartPointer<vtkMutableUndirectedGraph>::New();

	int cell_num = mesh->GetNumberOfCells();

	vtkSmartPointer<vtkDoubleArray> centers =
		vtkSmartPointer<vtkDoubleArray>::New();
	centers->SetName("Centers");
	centers->SetNumberOfComponents(3);

	vtkSmartPointer<vtkDoubleArray> areas =
		vtkSmartPointer<vtkDoubleArray>::New();
	areas->SetName("Areas");
	areas->SetNumberOfComponents(1);

	for (int i = 0; i < cell_num; ++i)
		dual_graph->AddVertex();

	for (int i = 0; i < cell_num; ++i) {
		vtkSmartPointer<vtkIdList> face_index =
			vtkSmartPointer<vtkIdList>::New();
		mesh->GetCellPoints(i, face_index);

		Vec3 p0, p1, p2;
		mesh->GetPoint(face_index->GetId(0), p0.data());
		mesh->GetPoint(face_index->GetId(1), p1.data());
		mesh->GetPoint(face_index->GetId(2), p2.data());

		double area = vtkTriangle::TriangleArea(p0.data(), p1.data(), p2.data());
		areas->InsertNextValue(area);

		Vec3 center = (p0 + p1 + p2) / 3;
		centers->InsertNextTuple(center.data());
	}

	vtkNew<vtkPolyDataNormals> normalFilter;
	normalFilter->SetInputData(mesh);
	normalFilter->ComputePointNormalsOff();
	normalFilter->ComputeCellNormalsOn();
	normalFilter->Update();

	vtkDataArray * normals = normalFilter->GetOutput()->GetCellData()->GetNormals();

	vtkSmartPointer<vtkDoubleArray> mesh_dis =
		vtkSmartPointer<vtkDoubleArray>::New();
	mesh_dis->SetName("Weights");
	mesh_dis->SetNumberOfComponents(1);

	vtkSmartPointer<vtkDoubleArray> phy_dis =
		vtkSmartPointer<vtkDoubleArray>::New();
	phy_dis->SetNumberOfComponents(1);

	vtkSmartPointer<vtkDoubleArray> angle_dis =
		vtkSmartPointer<vtkDoubleArray>::New();
	angle_dis->SetNumberOfComponents(1);

	vtkSmartPointer<vtkDoubleArray> edge_dis =
		vtkSmartPointer<vtkDoubleArray>::New();
	edge_dis->SetName("EdgeLens");
	edge_dis->SetNumberOfComponents(1);

	double phy_dis_avg = 0.0, angle_dis_avg = 0.0;
	for (int i = 0; i < cell_num; ++i) {
		vtkSmartPointer<vtkIdList> face_index =
			vtkSmartPointer<vtkIdList>::New();
		mesh->GetCellPoints(i, face_index);

		vector<int> neighbors;

		Vec3 p0, p1, p2;
		mesh->GetPoint(face_index->GetId(0), p0.data());
		mesh->GetPoint(face_index->GetId(1), p1.data());
		mesh->GetPoint(face_index->GetId(2), p2.data());

		double lateral[3];
		lateral[0] = sqrt(vtkMath::Distance2BetweenPoints(p0.data(), p1.data()));
		lateral[1] = sqrt(vtkMath::Distance2BetweenPoints(p1.data(), p2.data()));
		lateral[2] = sqrt(vtkMath::Distance2BetweenPoints(p2.data(), p0.data()));

		for (int j = 0; j < 3; ++j) {
			vtkSmartPointer<vtkIdList> id_list =
				vtkSmartPointer<vtkIdList>::New();
			id_list->InsertNextId(face_index->GetId(j));
			id_list->InsertNextId(face_index->GetId((j + 1) % 3));

			vtkSmartPointer<vtkIdList> neighbor_cell_ids =
				vtkSmartPointer<vtkIdList>::New();
			mesh->GetCellNeighbors(i, id_list, neighbor_cell_ids);
			for (int k = 0; k < neighbor_cell_ids->GetNumberOfIds(); ++k) {
				int neighbor_cell_id = neighbor_cell_ids->GetId(k);
				if (i >= neighbor_cell_id) continue;

				neighbors.push_back(neighbor_cell_id);

				double a, b;
				a = 2.0 * areas->GetValue(i) / (3 * lateral[j]);
				b = 2.0 * areas->GetValue(neighbor_cell_id) / (3 * lateral[j]);

				Vec3 n0, n1;
				normals->GetTuple(i, n0.data());
				normals->GetTuple(neighbor_cell_id, n1.data());

				Vec3 c0, c1;
				centers->GetTuple(i, c0.data());
				centers->GetTuple(neighbor_cell_id, c1.data());

				Vec3 w = c1 - c0;
				double phy = a + b;
				double angle = 0.0;
				if (n0 * w >= 0)
					angle = 1 - n0 * n1;

				phy_dis->InsertNextValue(phy);
				angle_dis->InsertNextValue(angle);
				edge_dis->InsertNextValue(lateral[j]);

				phy_dis_avg += phy;
				angle_dis_avg += angle;
			}
		}

		for (int j = 0; j < neighbors.size(); ++j)
			dual_graph->AddEdge(i, neighbors[j]);
	}

	int edge_num = phy_dis->GetNumberOfTuples();
	phy_dis_avg /= edge_num;
	angle_dis_avg /= edge_num;
	for (int i = 0; i < edge_num; ++i) {
		mesh_dis->InsertNextValue(
			phy_ratio * phy_dis->GetValue(i) / phy_dis_avg +
			(1 - phy_ratio) * angle_dis->GetValue(i) / angle_dis_avg
		);
	}

	dual_graph->GetEdgeData()->AddArray(mesh_dis);
	dual_graph->GetVertexData()->AddArray(centers);
	dual_graph->GetEdgeData()->AddArray(edge_dis);

	return dual_graph;
}

void MeshSegmenter::randomSelectSeeds() {
	using std::vector;

	vtkMath::RandomSeed(time(nullptr));

	vector<int> seed_map(this->mesh->GetNumberOfCells(), 0);
	for (int i = 0; i < this->seed_num; ++i) {
		int seed_id = (int)vtkMath::Random(0, this->mesh->GetNumberOfCells());
		while (seed_map[seed_id])
			seed_id = (int)vtkMath::Random(0, this->mesh->GetNumberOfCells());
		seed_map[seed_id] = 1;
		cluster_face_ids[i].push_back(seed_id);
	}
}

std::vector<double> MeshSegmenter::calcDijkstraTable(int face_id) {
	using std::vector;
	using std::make_tuple;
	using std::get;
	using HeapElem = std::tuple<int, double>;
	class HeapElemComp {
	public:
		bool operator() (const HeapElem & A, const HeapElem & B) {
			return get<1>(A) < get<1>(B) || (get<1>(A) == get<1>(B) && get<0>(A) < get<0>(B));
		}
	};

	vector<double> dijkstra_table(this->mesh->GetNumberOfCells(), DBL_MAX);
	vtkNew<vtkInEdgeIterator> it;
	this->dual_graph->GetInEdges(face_id, it.GetPointer());
	while (it->HasNext()) {
		vtkInEdgeType edge = it->Next();
		dijkstra_table[edge.Source] = mesh_dis->GetValue(edge.Id);
	}
	dijkstra_table[face_id] = 0.0;

	vector<int> S(this->mesh->GetNumberOfCells(), 0);
	S[face_id] = 1;

	vector<HeapElem> dis_pairs(this->mesh->GetNumberOfCells());
	for (int i = 0; i < this->mesh->GetNumberOfCells(); ++i)
		dis_pairs[i] = make_tuple(i, dijkstra_table[i]);
	
	MinHeap<HeapElem, HeapElemComp> min_heap(dis_pairs, this->mesh->GetNumberOfCells());
	min_heap.ExtractMin();

	vtkSmartPointer<vtkInEdgeIterator> uIt =
		vtkSmartPointer<vtkInEdgeIterator>::New();
	while (min_heap.Size()) {
		int u = get<0>(min_heap.ExtractMin());
		S[u] = 1;

		this->dual_graph->GetInEdges(u, uIt);
		while (uIt->HasNext()) {
			vtkInEdgeType uEdge = uIt->Next();
			vtkIdType v = uEdge.Source;
			if (S[v]) continue;

			double tmp = dijkstra_table[u] + mesh_dis->GetValue(uEdge.Id);
			if (dijkstra_table[v] > tmp) {
				dijkstra_table[v] = tmp;
				min_heap.DecreaseKey(make_tuple(v, tmp));
			}
		}
	}

	return dijkstra_table;
}

int computeHashValue(int a, int b, int seedCnt) {
	if (a < b) return a * seedCnt + b;
	else return b * seedCnt + a;
}

void MeshSegmenter::mergeClusters() {
	using std::abs;
	using std::array;
	using std::set;
	using std::unordered_map;
	using std::vector;
	using std::make_tuple;
	using std::get;
	
	using HeapElem = std::tuple<int, double>;
	class HeapElemComp {
	public:
		bool operator() (const HeapElem & A, const HeapElem & B) const {
			return get<1>(A) < get<1>(B) || (get<1>(A) == get<1>(B) && get<0>(A) < get<0>(B));
		}
	};

	struct Edge {
		int a, b;
		Edge(int a_, int b_) : a(a_), b(b_) {}
		bool operator == (const Edge & e) const { return (a == e.a && b == e.b) || (a == e.b && b == e.a); }
	};
	auto EdgeHash = [&](const Edge & e) {
		if (e.a < e.b) return e.a * this->seed_num + e.b;
		else return e.b * this->seed_num + e.a;
	};

	unordered_map<Edge, array<double, 5>, decltype(EdgeHash)> util_values(this->seed_num, EdgeHash);
	vtkSmartPointer<vtkEdgeListIterator> edgeIt =
		vtkSmartPointer<vtkEdgeListIterator>::New();

	this->dual_graph->GetEdges(edgeIt);
	while (edgeIt->HasNext()) {
		vtkEdgeType edge = edgeIt->Next();
		int cluster_num_a, cluster_num_b;
		cluster_num_a = this->face_id_to_cluster[edge.Source];
		cluster_num_b = this->face_id_to_cluster[edge.Target];

		Edge edge_ab(cluster_num_a, cluster_num_b);

		if (cluster_num_a == cluster_num_b ||
			cluster_num_a == -1 ||
			cluster_num_b == -1)
			continue;

		if (util_values.find(edge_ab) == util_values.end()) {
			util_values.insert(make_pair(edge_ab, array<double, 5>()));
			for (int i = 0; i < 5; ++i)
				util_values[edge_ab][i] = 0.0;
		}

		util_values[edge_ab][0] += edge_lens->GetValue(edge.Id) * mesh_dis->GetValue(edge.Id);
		util_values[edge_ab][1] += edge_lens->GetValue(edge.Id);
	}

	vector<array<double, 2>> sum_values(this->seed_num);
	for (int i = 0; i < this->seed_num; ++i) {
		double sum_D = 0.0, sum_L = 0.0;
		for (int j = 0; j < this->seed_num; ++j) {
			Edge edge_ij(i, j);
			if (util_values.find(edge_ij) != util_values.end()) {
				sum_D += util_values[edge_ij][0];
				sum_L += util_values[edge_ij][1];
			}
		}
		sum_values[i][0] = sum_D;
		sum_values[i][1] = sum_L;
	}

	set<HeapElem, HeapElemComp> min_heap;
	for (int i = 0; i < this->seed_num; ++i) for (int j = i + 1; j < this->seed_num; ++j) {
		Edge edge_ij(i, j);
		if (util_values.find(edge_ij) != util_values.end()) {
			util_values[edge_ij][2] = sum_values[i][0] + sum_values[j][0] - 2 * util_values[edge_ij][0];
			util_values[edge_ij][3] = sum_values[i][1] + sum_values[j][1] - 2 * util_values[edge_ij][1];
			util_values[edge_ij][4] = (util_values[edge_ij][0] / util_values[edge_ij][1]) / (util_values[edge_ij][2] / util_values[edge_ij][3]);
			min_heap.insert(make_tuple(i * this->seed_num + j, util_values[edge_ij][4]));
		}
	}

	this->cluster_steps = vector<vector<int>>(this->seed_num);
	for (int i = 0; i < this->seed_num; ++i) {
		this->cluster_steps[i] = vector<int>(this->seed_num, -1);
		this->cluster_steps[0][i] = i;
	}

	int remain_cluster_num = this->seed_num;
	while (remain_cluster_num > 2) {
		int cluster_num_a, cluster_num_b;
		cluster_num_a = get<0>(*min_heap.begin()) / this->seed_num;
		cluster_num_b = get<0>(*min_heap.begin()) % this->seed_num;

		Edge edge_ab(cluster_num_a, cluster_num_b);

		sum_values[cluster_num_a][0] = util_values[edge_ab][2];
		sum_values[cluster_num_a][1] = util_values[edge_ab][3];
		for (int i = 0; i < this->seed_num; ++i) {
			if (i == cluster_num_a ||
				i == cluster_num_b)
				continue;

			Edge edge_ai(cluster_num_a, i);
			Edge edge_bi(cluster_num_b, i);

			if (util_values.find(edge_ai) != util_values.end() && 
				util_values.find(edge_bi) != util_values.end()) {

				util_values[edge_ai][0] += util_values[edge_bi][0];
				util_values[edge_ai][1] += util_values[edge_bi][1];
				util_values[edge_ai][2] = sum_values[cluster_num_a][0] + sum_values[i][0] - 2 * util_values[edge_ai][0];
				util_values[edge_ai][3] = sum_values[cluster_num_a][1] + sum_values[i][1] - 2 * util_values[edge_ai][1];

				min_heap.erase(min_heap.find(make_tuple(computeHashValue(cluster_num_a, i, this->seed_num), util_values[edge_ai][4])));
				min_heap.erase(min_heap.find(make_tuple(computeHashValue(cluster_num_b, i, this->seed_num), util_values[edge_bi][4])));
				util_values.erase(edge_bi);
			} else if (util_values.find(edge_ai) != util_values.end() &&
				util_values.find(edge_bi) == util_values.end()) {

				util_values[edge_ai][2] = sum_values[cluster_num_a][0] + sum_values[i][0] - 2 * util_values[edge_ai][0];
				util_values[edge_ai][3] = sum_values[cluster_num_a][1] + sum_values[i][1] - 2 * util_values[edge_ai][1];

				min_heap.erase(min_heap.find(make_tuple(computeHashValue(cluster_num_a, i, this->seed_num), util_values[edge_ai][4])));
			} else if (util_values.find(edge_ai) == util_values.end() &&
				util_values.find(edge_bi) != util_values.end()) {

				util_values.insert(make_pair(edge_ai, array<double, 5>()));
				util_values[edge_ai][0] = util_values[edge_bi][0];
				util_values[edge_ai][1] = util_values[edge_bi][1];
				util_values[edge_ai][2] = sum_values[cluster_num_a][0] + sum_values[i][0] - 2 * util_values[edge_ai][0];
				util_values[edge_ai][3] = sum_values[cluster_num_a][1] + sum_values[i][1] - 2 * util_values[edge_ai][1];

				min_heap.erase(min_heap.find(make_tuple(computeHashValue(cluster_num_b, i, this->seed_num), util_values[edge_bi][4])));
				util_values.erase(edge_bi);
			} else
				continue;

			if (abs(util_values[edge_ai][1] * util_values[edge_ai][2]) < 1e-3) {
				util_values[edge_ai][4] = DBL_MAX;
			} else {
				util_values[edge_ai][4] =
					(util_values[edge_ai][0] * util_values[edge_ai][3]) /
					(util_values[edge_ai][1] * util_values[edge_ai][2]);
			}
			util_values[edge_ai][4] = util_values[edge_ai][4];
			min_heap.insert(make_tuple(computeHashValue(cluster_num_a, i, this->seed_num), util_values[edge_ai][4]));
		}

		min_heap.erase(make_tuple(computeHashValue(cluster_num_a, cluster_num_b, this->seed_num), util_values[edge_ab][4]));
		util_values.erase(edge_ab);

		--remain_cluster_num;

		for (int i = 0; i < this->seed_num; ++i) {
			if (cluster_steps[this->seed_num - remain_cluster_num - 1][i] == cluster_num_b)
				cluster_steps[this->seed_num - remain_cluster_num][i] = cluster_num_a;
			else
				cluster_steps[this->seed_num - remain_cluster_num][i] = cluster_steps[seed_num - remain_cluster_num - 1][i];
		}
	}
}
