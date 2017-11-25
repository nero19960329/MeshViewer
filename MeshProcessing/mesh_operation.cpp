#include "mesh_operation.h"

#include <vtkIdList.h>

#include <unordered_set>

std::vector<vtkIdType> MeshOperation::getConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, vtkIdType id) {
	using std::vector;
	using std::unordered_set;

	vector<vtkIdType> connectedVertices;
	unordered_set<vtkIdType> used_set;
	used_set.insert(id);

	vtkSmartPointer<vtkIdList> cellIdList =
		vtkSmartPointer<vtkIdList>::New();
	mesh->GetPointCells(id, cellIdList);

	for (vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); ++i) {
		vtkSmartPointer<vtkIdList> pointIdList =
			vtkSmartPointer<vtkIdList>::New();
		mesh->GetCellPoints(cellIdList->GetId(i), pointIdList);

		for (vtkIdType j = 0; j < pointIdList->GetNumberOfIds(); ++j) {
			int target_id = pointIdList->GetId(j);
			if (used_set.find(target_id) == used_set.end()) {
				used_set.insert(target_id);
				connectedVertices.push_back(target_id);
			}
		}
	}

	return connectedVertices;
}

std::vector<vtkIdType> MeshOperation::getVertexConnectedFaces(vtkSmartPointer<vtkPolyData> mesh, vtkIdType id) {
	using std::vector;

	vector<vtkIdType> connectedFaces;

	vtkSmartPointer<vtkIdList> cellIdList =
		vtkSmartPointer<vtkIdList>::New();
	mesh->GetPointCells(id, cellIdList);

	for (vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); ++i)
		connectedFaces.push_back(cellIdList->GetId(i));

	return connectedFaces;
}

std::vector<vtkIdType> MeshOperation::getFaceConnectedFaces(vtkSmartPointer<vtkPolyData> mesh, vtkIdType id) {
	using std::vector;

	vector<vtkIdType> connectedFaces;

	vtkSmartPointer<vtkIdList> cellPointIds =
		vtkSmartPointer<vtkIdList>::New();
	mesh->GetCellPoints(id, cellPointIds);

	for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); ++i) {
		vtkSmartPointer<vtkIdList> idList =
			vtkSmartPointer<vtkIdList>::New();
		idList->InsertNextId(cellPointIds->GetId(i));

		if (i + 1 == cellPointIds->GetNumberOfIds())
			idList->InsertNextId(cellPointIds->GetId(0));
		else
			idList->InsertNextId(cellPointIds->GetId(i + 1));

		vtkSmartPointer<vtkIdList> neighborCellIds =
			vtkSmartPointer<vtkIdList>::New();
		mesh->GetCellNeighbors(id, idList, neighborCellIds);

		for (vtkIdType j = 0; j < neighborCellIds->GetNumberOfIds(); ++j)
			connectedFaces.push_back(neighborCellIds->GetId(j));
	}

	return connectedFaces;
}
