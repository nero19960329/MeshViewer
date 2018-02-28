#pragma once

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <QtCore/QString>

#include <vector>

class ColorTableHandler {
private:
	vtkSmartPointer<vtkPolyData> mesh;
	QString color_table_name;
	std::vector<double> color_value_vec;
	double max_scalar, min_scalar;

public:
	ColorTableHandler() : mesh(nullptr), max_scalar(DBL_MIN), min_scalar(DBL_MAX) {}
	
	void setMesh(vtkSmartPointer<vtkPolyData> mesh_) { mesh = mesh_; }
	void setColorTableName(const QString & color_table_name_) { color_table_name = color_table_name_; }

	bool read();
	bool readColorValueVec(const std::vector<double> & color_value_vec_);
	bool write();

	vtkSmartPointer<vtkPolyData> turnToContinuous();
	vtkSmartPointer<vtkPolyData> turnToDiscrete();

	double maxScalar() const { return max_scalar; }
	double minScalar() const { return min_scalar; }
};