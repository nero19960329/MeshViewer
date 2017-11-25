#pragma once

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <vector>

class MeshOperation {
public:
	static std::vector<vtkIdType> getConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, vtkIdType id);
	static std::vector<vtkIdType> getVertexConnectedFaces(vtkSmartPointer<vtkPolyData> mesh, vtkIdType id);
	static std::vector<vtkIdType> getFaceConnectedFaces(vtkSmartPointer<vtkPolyData> mesh, vtkIdType id);
};