#include "color_table_reader.h"

#include <fstream>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>

bool ColorTableReader::read() {
	using std::ifstream;

	ifstream fin(this->color_table_name.toLocal8Bit());
	if (!fin.good()) return false;

	color_value_vec.clear();
	while (!fin.eof()) {
		double tmp;
		fin >> tmp;
		if (tmp > this->max_scalar)
			this->max_scalar = tmp;
		if (tmp < this->min_scalar)
			this->min_scalar = tmp;
		color_value_vec.push_back(tmp);
	}
	color_value_vec.pop_back();

	if (color_value_vec.size() != this->mesh->GetNumberOfCells())
		return false;

	return true;
}

vtkSmartPointer<vtkPolyData> ColorTableReader::turnToContinuous() {
	if (mesh == nullptr) return nullptr;

	vtkSmartPointer<vtkDoubleArray> scalars =
		vtkSmartPointer<vtkDoubleArray>::New();
	scalars->SetNumberOfValues(this->mesh->GetNumberOfCells());
	for (int i = 0; i < this->mesh->GetNumberOfCells(); ++i)
		scalars->SetValue(i, color_value_vec[i]);
	mesh->GetCellData()->SetScalars(scalars);
	return mesh;
}

vtkSmartPointer<vtkPolyData> ColorTableReader::turnToDiscrete() {
	if (mesh == nullptr) return nullptr;

	vtkSmartPointer<vtkDoubleArray> scalars =
		vtkSmartPointer<vtkDoubleArray>::New();
	scalars->SetNumberOfValues(this->mesh->GetNumberOfCells());
	for (int i = 0; i < this->mesh->GetNumberOfCells(); ++i)
		scalars->SetValue(i, color_value_vec[i] - min_scalar);
	mesh->GetCellData()->SetScalars(scalars);
	return mesh;
}
