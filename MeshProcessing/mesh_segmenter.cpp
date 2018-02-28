#include "fibonacci_heap.h"
#include "mesh_segmenter.h"
#include "min_heap.h"
#include "timer.h"
#include "vec3.h"

#include <vtkBoundingBox.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkDataArray.h>
#include <vtkEdgeListIterator.h>
#include <vtkIdList.h>
#include <vtkInEdgeIterator.h>
#include <vtkMath.h>
#include <vtkNew.h>
#include <vtkOBBTree.h>
#include <vtkOctreePointLocator.h>
#include <vtkPointData.h>
#include <vtkPolyDataNormals.h>
#include <vtkTriangle.h>

#include <array>
#include <queue>
#include <set>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

#include "vtk_widget.h"

MeshSegmenter::MeshSegmenter(int seed_num_, double phy_ratio_) : seed_num(seed_num_), phy_ratio(phy_ratio_) {
	using std::vector;

	this->cluster_face_ids = vector<vector<int>>(this->seed_num);
}

void MeshSegmenter::segment() {
	using std::vector;

	this->dual_graph = this->calcDualGraph(this->phy_ratio);
	this->centers = vtkDoubleArray::SafeDownCast(dual_graph->GetVertexData()->GetArray("Centers"));
	this->mesh_dis = vtkDoubleArray::SafeDownCast(dual_graph->GetEdgeData()->GetArray("Weights"));
	this->edge_lens = vtkDoubleArray::SafeDownCast(dual_graph->GetEdgeData()->GetArray("EdgeLens"));

	Timer timer;
	timer.begin();
	//this->selectSeedsRandomly();
	this->selectSeedsByOctree();
	std::cout << "Select seeds time cost = " << timer.getDuration() << " ms" << std::endl;

	/*vector<vector<double>> dists;
	for (int i = 0; i < this->seed_num; ++i)
		dists.push_back(vector<double>());
	timer.begin();
	#pragma omp parallel for
	for (int i = 0; i < seed_num; ++i)
		dists[i] = this->calcDijkstraTableWithMinHeap(this->cluster_face_ids[i].front());
	double duration = timer.getDuration();
	std::cout << "all = " << duration << " ms" << std::endl;*/
	//std::cout << "average = " << (duration / this->seed_num) << " ms" << std::endl;

	vector<vector<double>> dists;
	timer.begin();
	int thread_num = 4;
	vector<int> seed_bounds;
	for (int i = 0; i < thread_num; ++i)
		seed_bounds.push_back(seed_num * i / thread_num);
	seed_bounds.push_back(seed_num);
	vector<vector<vector<double>>> dists_per_thread(thread_num);
	#pragma omp parallel for
	for (int i = 0; i < thread_num; ++i) {
		vector<int> face_ids;
		for (int j = seed_bounds[i]; j < seed_bounds[i + 1]; ++j)
			face_ids.push_back(this->cluster_face_ids[j].front());
		dists_per_thread[i] = this->calcMultiDijkstraTable(face_ids);
	}
	for (int i = 0; i < thread_num; ++i)
		dists.insert(dists.end(), dists_per_thread[i].begin(), dists_per_thread[i].end());
	double duration = timer.getDuration();
	std::cout << "all = " << duration << " ms" << std::endl;

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
	using std::unordered_map;
	using std::make_pair;

	vtkSmartPointer<vtkDoubleArray> scalars =
		vtkSmartPointer<vtkDoubleArray>::New();
	scalars->SetNumberOfValues(this->mesh->GetNumberOfCells());

	unordered_map<int, int> cluster_convert_map;
	for (int i = 0; i < this->mesh->GetNumberOfCells(); ++i) {
		int cluster_id = this->cluster_steps[this->seed_num - n][this->face_id_to_cluster[i]];
		if (cluster_convert_map.find(cluster_id) == cluster_convert_map.end())
			cluster_convert_map.insert(make_pair(cluster_id, cluster_convert_map.size()));
		scalars->SetValue(i, cluster_convert_map[cluster_id]);
	}

	for (int i = 0; i < this->mesh->GetNumberOfCells(); ++i)
		scalars->SetValue(i, scalars->GetValue(i)/* * 1.0 / cluster_convert_map.size()*/);

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

vtkSmartPointer<vtkPolyData> MeshSegmenter::removeAbnormalRegion(vtkSmartPointer<vtkPolyData> mesh) {
	using std::max;
	using std::vector;

	vtkNew<vtkPolyDataNormals> normalFilter;
	normalFilter->SetInputData(mesh);
	normalFilter->ComputePointNormalsOn();
	normalFilter->ComputeCellNormalsOn();
	normalFilter->Update();

	vtkDataArray * cellNormals = normalFilter->GetOutput()->GetCellData()->GetNormals();
	vtkDataArray * pointNormals = normalFilter->GetOutput()->GetPointData()->GetNormals();

	int n = mesh->GetNumberOfPoints();
	vtkBoundingBox boundingBox;
	for (int i = 0; i < n; ++i)
		boundingBox.AddPoint(mesh->GetPoint(i));
	double longest = boundingBox.GetMaxLength() * 1e-2;
	double eps = longest * 1e-3;
	std::cout << longest << std::endl;

	vtkNew<vtkOBBTree> tree;
	tree->SetDataSet(mesh);
	tree->BuildLocator();

	vector<vtkIdType> deletePointIds;
	for (int i = 0; i < n; ++i) {
		double pointNormal[3];
		pointNormals->GetTuple(i, pointNormal);

		double p0[3], p1[3];
		mesh->GetPoint(i, p0);
		for (int j = 0; j < 3; ++j) {
			p0[j] -= eps * pointNormal[j];
			p1[j] = p0[j] - longest * pointNormal[j];
		}

		vtkSmartPointer<vtkPoints> intersectPoints =
			vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkIdList> cellIds =
			vtkSmartPointer<vtkIdList>::New();
		tree->IntersectWithLine(p0, p1, intersectPoints, cellIds);

		if (intersectPoints->GetNumberOfPoints() == 0)
			continue;

		double minDist2 = DBL_MAX, minCellId = -1;
		for (int j = 0; j < cellIds->GetNumberOfIds(); ++j) {
			double p2[3];
			intersectPoints->GetPoint(j, p2);
			double dist2 = vtkMath::Distance2BetweenPoints(p0, p2);
			if (dist2 < minDist2) {
				minDist2 = dist2;
				minCellId = cellIds->GetId(j);
			}
		}

		double cellNormal[3];
		cellNormals->GetTuple(minCellId, cellNormal);

		double cosine = vtkMath::Dot(cellNormal, pointNormal) / (vtkMath::Norm(cellNormal) * vtkMath::Norm(pointNormal));
		if (cosine >= -0.95)
			continue;

		std::cout << cosine << ", ";
		deletePointIds.push_back(i);
		std::cout << i << std::endl;

		std::cout << "p0 = " << p0[0] << ", " << p0[1] << ", " << p0[2] << std::endl;
		std::cout << "p1 = " << p1[0] << ", " << p1[1] << ", " << p1[2] << std::endl;
		if (i >= 1000)
			break;
		this->vtk_widget_->addLine(p0, p1);
	}

	vtkSmartPointer<vtkPolyData> res =
		vtkSmartPointer<vtkPolyData>::New();
	res->DeepCopy(mesh);

	return res;
}

void MeshSegmenter::selectSeedsRandomly() {
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

void MeshSegmenter::selectSeedsByOctree() {
	using std::vector;
	using std::abs;
	using std::unordered_set;

	vtkSmartPointer<vtkOctreePointLocator> octree =
		vtkSmartPointer<vtkOctreePointLocator>::New();
	octree->SetDataSet(this->mesh);
	octree->SetMaximumPointsPerRegion(this->mesh->GetNumberOfPoints() / this->seed_num);
	octree->BuildLocator();

	int n = octree->GetNumberOfLeafNodes();
	/*std::cout << "octree->leaf_nodes_number = " << octree->GetNumberOfLeafNodes() << std::endl;
	int p = 0, cnt = 0;
	for (int i = 0; i < n; ++i) {
		int k = octree->GetPointsInRegion(i)->GetNumberOfValues();
		if (k != 0) {
			++cnt;
			double bounds[6];
			octree->GetRegionDataBounds(i, bounds);
			std::cout << i << " : " << k << std::endl;
			for (int j = 0; j < 6; ++j)
				std::cout << bounds[j] << " ";
			std::cout << std::endl;
		}
		p += k;
	}
	std::cout << "leaf_node_number_all = " << p << std::endl;
	std::cout << "non_zero_leaf_node_number = " << cnt << std::endl;
	std::cout << "octree->maximuum_points_per_region = " << octree->GetMaximumPointsPerRegion() << std::endl;*/

	/*vtkSmartPointer<vtkDoubleArray> scalars =
		vtkSmartPointer<vtkDoubleArray>::New();
	scalars->SetNumberOfValues(this->mesh->GetNumberOfPoints());
	for (int i = 0; i < mesh->GetNumberOfPoints(); ++i) {
		double p[3];
		mesh->GetPoint(i, p);
		scalars->SetValue(i, octree->GetRegionContainingPoint(p[0], p[1], p[2]) * 1.0 / n);
	}
	mesh->GetPointData()->SetScalars(scalars);*/

	double eps_length = this->vtk_widget_->mesh_processing_data_model_->mean_edge_length;

	vector<int> nonEmptyLeafs;
	for (int i = 0; i < n; ++i) {
		if (octree->GetPointsInRegion(i)->GetNumberOfValues() != 0)
			nonEmptyLeafs.push_back(i);
	}

	auto isBoundConnect = [&](double x1, double x2, double x1_, double x2_, double eps) -> bool {
		if (abs(x1 - x1_) < eps ||
			abs(x1 - x2_) < eps ||
			abs(x2 - x1_) < eps ||
			abs(x2 - x2_) < eps) return true;
		return false;
	};

	// Floyd algorithm
	int m = nonEmptyLeafs.size();
	vector<vector<int>> A;
	A.resize(m);
	for (int i = 0; i < m; ++i) {
		A[i].resize(m);
		for (int j = 0; j < m; ++j)
			A[i][j] = m;
		A[i][i] = 0;
	}

	int edge_cnt = 0;
	for (int i = 0; i < m; ++i) {
		double bounds1[6];
		octree->GetRegionDataBounds(nonEmptyLeafs[i], bounds1);
		for (int j = i + 1; j < m; ++j) {
			double bounds2[6];
			octree->GetRegionDataBounds(nonEmptyLeafs[j], bounds2);

			if (isBoundConnect(bounds1[0], bounds1[1], bounds2[0], bounds2[1], eps_length) &&
				isBoundConnect(bounds1[2], bounds1[3], bounds2[2], bounds2[3], eps_length) &&
				isBoundConnect(bounds1[4], bounds1[5], bounds2[4], bounds2[5], eps_length)) {
				A[i][j] = 1;
				A[j][i] = 1;
				++edge_cnt;
			}
		}
	}

	for (int k = 0; k < m; ++k) {
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < m; ++j) {
				if (A[i][j] > A[i][k] + A[k][j])
					A[i][j] = A[i][k] + A[k][j];
			}
		}
	}

	unordered_set<int> isoSet;
	for (int i = 0; i < m; ++i) {
		bool isIso = true;
		for (int j = 0; j < m; ++j) {
			if (A[i][j] < m) {
				isIso = false;
				break;
			}
		}

		if (isIso)
			isoSet.insert(i);
	}

	unordered_set<int> seedLeafIds;
	int max_dist = 0, seed1 = -1, seed2 = -1;
	for (int i = 0; i < m; ++i) {
		if (isoSet.find(i) != isoSet.end()) break;
		for (int j = 0; j < m; ++j) {
			if (isoSet.find(j) != isoSet.end()) break;
			if (A[i][j] >= m) continue;

			if (max_dist < A[i][j]) {
				max_dist = A[i][j];
				seed1 = i;
				seed2 = j;
			}
		}
	}
	seedLeafIds.insert(seed1);
	seedLeafIds.insert(seed2);
	std::cout << max_dist << std::endl;

	vtkMath::RandomSeed(time(nullptr));

	// farthest point sampling
	/*for (int k = 2; k < this->seed_num; ++k) {
		int max_dist_sum = 0;
		vector<int> seed_vec;
		vector<int> dist_sum_vec;
		dist_sum_vec.resize(m);

		for (int i = 0; i < m; ++i) {
			dist_sum_vec[i] = 0;

			if (isoSet.find(i) != isoSet.end()) break;
			if (seedLeafIds.find(i) != seedLeafIds.end()) continue;

			int dist_sum = 0;
			for (auto it = seedLeafIds.begin(); it != seedLeafIds.end(); ++it) {
				dist_sum += A[i][*it];
			}

			if (max_dist_sum < dist_sum)
				max_dist_sum = dist_sum;

			dist_sum_vec[i] = dist_sum;
		}

		for (int i = 0; i < m; ++i) {
			if (dist_sum_vec[i] == max_dist_sum)
				seed_vec.push_back(i);
		}

		int seed = seed_vec[(int)vtkMath::Random(0, seed_vec.size())];
		seedLeafIds.insert(seed);
	}*/

	// coarse sampling
	for (int k = 2; k < this->seed_num; ++k) {
		int max_min_dist = 0;
		vector<int> seed_vec;
		vector<double> min_dist_vec;
		min_dist_vec.resize(m);

		for (int i = 0; i < m; ++i) {
			min_dist_vec[i] = DBL_MAX;

			if (isoSet.find(i) != isoSet.end()) break;
			if (seedLeafIds.find(i) != seedLeafIds.end()) continue;

			int min_dist = INT_MAX;
			for (auto it = seedLeafIds.begin(); it != seedLeafIds.end(); ++it) {
				if (min_dist > A[i][*it])
					min_dist = A[i][*it];
			}

			if (max_min_dist < min_dist)
				max_min_dist = min_dist;

			min_dist_vec[i] = min_dist;
		}

		for (int i = 0; i < m; ++i) {
			if (min_dist_vec[i] == max_min_dist)
				seed_vec.push_back(i);
		}

		int seed = seed_vec[(int) vtkMath::Random(0, seed_vec.size())];
		seedLeafIds.insert(seed);
	}

	int i = 0;
	for (auto it = seedLeafIds.begin(); it != seedLeafIds.end(); ++it) {
		vtkIdTypeArray * pts = octree->GetPointsInRegion(nonEmptyLeafs[*it]);
		int k = pts->GetNumberOfValues();
		int seedPointId = pts->GetValue((int) vtkMath::Random(0, k));
		vtkSmartPointer<vtkIdList> cellIdList =
			vtkSmartPointer<vtkIdList>::New();
		mesh->GetPointCells(seedPointId, cellIdList);
		int seedCellId = cellIdList->GetId((int) vtkMath::Random(0, cellIdList->GetNumberOfIds()));
		cluster_face_ids[i].push_back(seedCellId);

		/*if (i < 3) {
			this->vtk_widget_->highlightVertex(mesh, seedPointId);
			this->vtk_widget_->highlightFace(mesh, seedCellId);
		}*/
		
		++i;
	}
}

std::vector<std::vector<double>> MeshSegmenter::calcMultiDijkstraTable(std::vector<int> face_ids) {
	using std::vector;
	using std::make_tuple;
	using std::get;
	using std::tie;
	using HeapElem = std::tuple<int, double>;
	class HeapElemComp {
	public:
		bool operator() (const HeapElem & A, const HeapElem & B) {
			return get<1>(A) < get<1>(B);
		}
	};

	int n = this->mesh->GetNumberOfCells();
	int k = face_ids.size();

	vector<vector<double>> dijkstra_table;
	dijkstra_table.resize(k);

	vector<vector<int>> previous_table;
	previous_table.resize(k);

	for (int i = 0; i < k; ++i) {
		dijkstra_table[i].resize(n);
		previous_table[i].resize(n);
		for (int j = 0; j < n; ++j) {
			dijkstra_table[i][j] = DBL_MAX;
			previous_table[i][j] = -1;
		}

		int face_id = face_ids[i];
		dijkstra_table[i][face_id] = 0.0;

		vtkNew<vtkInEdgeIterator> it;
		this->dual_graph->GetInEdges(face_id, it.GetPointer());
		while (it->HasNext()) {
			vtkInEdgeType edge = it->Next();
			dijkstra_table[i][edge.Source] = mesh_dis->GetValue(edge.Id);
			previous_table[i][edge.Source] = face_id;
		}

		vector<int> S(n, 0), C(n, 0);
		S[face_id] = 1;

		vector<HeapElem> dis_pairs(n);
		for (int j = 0; j < n; ++j)
			dis_pairs[j] = make_tuple(j, dijkstra_table[i][j]);

		MinHeap<HeapElem, HeapElemComp> min_heap(dis_pairs, n);
		min_heap.extractMin();

		vtkNew<vtkInEdgeIterator> uIt;
		while (min_heap.size()) {
			int u;
			double d;
			tie(u, d) = min_heap.extractMin();

			if (d == DBL_MAX) break;
			S[u] = 1;

			if (C[previous_table[i][u]]) continue;
			for (int j = 0; j < i; ++j) {
				if (dijkstra_table[i][u] > dijkstra_table[j][u]) {
					C[u] = 1;
					break;
				}
			}
			if (C[u]) continue;

			this->dual_graph->GetInEdges(u, uIt.GetPointer());
			while (uIt->HasNext()) {
				vtkInEdgeType uEdge = uIt->Next();
				vtkIdType v = uEdge.Source;
				if (S[v]) continue;

				double tmp = dijkstra_table[i][u] + mesh_dis->GetValue(uEdge.Id);
				if (dijkstra_table[i][v] > tmp) {
					dijkstra_table[i][v] = tmp;
					previous_table[i][v] = u;
					min_heap.decreaseKey(make_tuple(v, tmp));
				}
			}
		}
	}

	return dijkstra_table;
}

std::vector<double> MeshSegmenter::calcDijkstraTableWithMinHeap(int face_id) {
	using std::vector;
	using std::make_tuple;
	using std::get;
	using HeapElem = std::tuple<int, double>;
	class HeapElemComp {
	public:
		bool operator() (const HeapElem & A, const HeapElem & B) {
			return get<1>(A) < get<1>(B);
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
	min_heap.extractMin();

	vtkSmartPointer<vtkInEdgeIterator> uIt =
		vtkSmartPointer<vtkInEdgeIterator>::New();
	while (min_heap.size()) {
		int u = get<0>(min_heap.extractMin());
		S[u] = 1;

		this->dual_graph->GetInEdges(u, uIt);
		while (uIt->HasNext()) {
			vtkInEdgeType uEdge = uIt->Next();
			vtkIdType v = uEdge.Source;
			if (S[v]) continue;

			double tmp = dijkstra_table[u] + mesh_dis->GetValue(uEdge.Id);
			if (dijkstra_table[v] > tmp) {
				dijkstra_table[v] = tmp;
				min_heap.decreaseKey(make_tuple(v, tmp));
			}
		}
	}

	return dijkstra_table;
}

std::vector<double> MeshSegmenter::calcDijkstraTableWithFibHeap(int face_id) {
	using std::vector;
	using std::make_tuple;
	using std::get;
	using HeapElem = std::tuple<int, double>;
	class HeapElemComp {
	public:
		bool operator() (const HeapElem & A, const HeapElem & B) {
			return get<1>(A) < get<1>(B);
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
	
	FibHeap<HeapElem, HeapElemComp> fib_heap(dis_pairs, this->mesh->GetNumberOfCells());
	fib_heap.extractMin();

	vtkSmartPointer<vtkInEdgeIterator> uIt =
		vtkSmartPointer<vtkInEdgeIterator>::New();
	while (fib_heap.size()) {
		int u = get<0>(fib_heap.extractMin()->key);
		S[u] = 1;

		this->dual_graph->GetInEdges(u, uIt);
		while (uIt->HasNext()) {
			vtkInEdgeType uEdge = uIt->Next();
			vtkIdType v = uEdge.Source;
			if (S[v]) continue;

			double tmp = dijkstra_table[u] + mesh_dis->GetValue(uEdge.Id);
			if (dijkstra_table[v] > tmp) {
				dijkstra_table[v] = tmp;
				fib_heap.decreaseKey(make_tuple(v, tmp));
			}
		}
	}

	return dijkstra_table;
}

std::vector<double> MeshSegmenter::calcSPFATable(int face_id) {
	using std::vector;
	using std::queue;

	vector<double> spfa_table(this->mesh->GetNumberOfCells(), DBL_MAX);
	vector<int> vis(this->mesh->GetNumberOfCells(), 0);
	spfa_table[face_id] = 0;
	vis[face_id] = 1;

	queue<int> q;
	q.push(face_id);

	vtkSmartPointer<vtkInEdgeIterator> vIt =
		vtkSmartPointer<vtkInEdgeIterator>::New();
	while (!q.empty()) {
		int v = q.front();
		q.pop();
		vis[v] = 0;

		this->dual_graph->GetInEdges(v, vIt);
		while (vIt->HasNext()) {
			vtkInEdgeType vEdge = vIt->Next();
			vtkIdType i = vEdge.Source;
			if (spfa_table[i] > spfa_table[v] + mesh_dis->GetValue(vEdge.Id)) {
				spfa_table[i] = spfa_table[v] + mesh_dis->GetValue(vEdge.Id);
				if (vis[i] == 0) {
					vis[i] = 1;
					q.push(i);
				}
			}
		}
	}

	return spfa_table;
}

void MeshSegmenter::mergeClusters() {
	using std::abs;
	using std::array;
	using std::set;
	using std::unordered_map;
	using std::vector;
	using std::make_tuple;
	using std::get;

	struct Edge {
		int a, b;
		Edge(int a_, int b_) : a(a_), b(b_) {}
		bool operator == (const Edge & e) const { return (a == e.a && b == e.b) || (a == e.b && b == e.a); }
	};
	auto EdgeHash = [&](const Edge & e) {
		if (e.a < e.b) return e.a * this->seed_num + e.b;
		else return e.b * this->seed_num + e.a;
	};

	using HeapElem = std::tuple<Edge, double>;
	class HeapElemComp {
	public:
		bool operator() (const HeapElem & A, const HeapElem & B) const {
			return get<1>(A) < get<1>(B);
		}
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
			min_heap.insert(make_tuple(edge_ij, util_values[edge_ij][4]));
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
		cluster_num_a = std::min(get<0>(*min_heap.begin()).a, get<0>(*min_heap.begin()).b);
		cluster_num_b = std::max(get<0>(*min_heap.begin()).a, get<0>(*min_heap.begin()).b);

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

				min_heap.erase(make_tuple(edge_ai, util_values[edge_ai][4]));
				min_heap.erase(make_tuple(edge_bi, util_values[edge_bi][4]));
				util_values.erase(edge_bi);
			} else if (util_values.find(edge_ai) != util_values.end() &&
				util_values.find(edge_bi) == util_values.end()) {

				util_values[edge_ai][2] = sum_values[cluster_num_a][0] + sum_values[i][0] - 2 * util_values[edge_ai][0];
				util_values[edge_ai][3] = sum_values[cluster_num_a][1] + sum_values[i][1] - 2 * util_values[edge_ai][1];

				min_heap.erase(make_tuple(edge_ai, util_values[edge_ai][4]));
			} else if (util_values.find(edge_ai) == util_values.end() &&
				util_values.find(edge_bi) != util_values.end()) {

				util_values.insert(make_pair(edge_ai, array<double, 5>()));
				util_values[edge_ai][0] = util_values[edge_bi][0];
				util_values[edge_ai][1] = util_values[edge_bi][1];
				util_values[edge_ai][2] = sum_values[cluster_num_a][0] + sum_values[i][0] - 2 * util_values[edge_ai][0];
				util_values[edge_ai][3] = sum_values[cluster_num_a][1] + sum_values[i][1] - 2 * util_values[edge_ai][1];

				min_heap.erase(make_tuple(edge_bi, util_values[edge_bi][4]));
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
			min_heap.insert(make_tuple(edge_ai, util_values[edge_ai][4]));
		}

		min_heap.erase(make_tuple(edge_ab, util_values[edge_ab][4]));
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
