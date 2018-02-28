#include "color_table_handler.h"

#include <fstream>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>

bool ColorTableHandler::read() {
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
	fin.close();

	if (color_value_vec.size() != this->mesh->GetNumberOfCells())
		return false;

	return true;
}

bool ColorTableHandler::readColorValueVec(const std::vector<double>& color_value_vec_) {
	color_value_vec.clear();
	for (const auto & val : color_value_vec_) {
		if (val > this->max_scalar)
			this->max_scalar = val;
		if (val < this->min_scalar)
			this->min_scalar = val;
		color_value_vec.push_back(val);
	}

	if (color_value_vec.size() != this->mesh->GetNumberOfCells())
		return false;

	return true;
}

bool ColorTableHandler::write() {
	using std::ofstream;
	using std::endl;

	ofstream fout(this->color_table_name.toLocal8Bit());
	for (const auto & val : color_value_vec)
		fout << val << endl;
	fout.flush();
	fout.close();

	return true;
}

vtkSmartPointer<vtkPolyData> ColorTableHandler::turnToContinuous() {
	if (mesh == nullptr) return nullptr;

	vtkSmartPointer<vtkDoubleArray> scalars =
		vtkSmartPointer<vtkDoubleArray>::New();
	scalars->SetNumberOfValues(this->mesh->GetNumberOfCells());
	for (int i = 0; i < this->mesh->GetNumberOfCells(); ++i)
		scalars->SetValue(i, color_value_vec[i]);
	mesh->GetCellData()->SetScalars(scalars);
	return mesh;
}

vtkSmartPointer<vtkPolyData> ColorTableHandler::turnToDiscrete() {
	if (mesh == nullptr) return nullptr;

	vtkSmartPointer<vtkDoubleArray> scalars =
		vtkSmartPointer<vtkDoubleArray>::New();
	scalars->SetNumberOfValues(this->mesh->GetNumberOfCells());
	for (int i = 0; i < this->mesh->GetNumberOfCells(); ++i)
		scalars->SetValue(i, color_value_vec[i] - min_scalar);
	mesh->GetCellData()->SetScalars(scalars);
	return mesh;
}
