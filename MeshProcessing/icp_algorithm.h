#pragma once

#include <vtkCellLocator.h>
#include <vtkMatrix4x4.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>

class ICPAlgorithm {
private:
	vtkSmartPointer<vtkPolyData> source, target;
	vtkSmartPointer<vtkTransform> transform;
	bool is_move_center;
	int max_iter;

	vtkSmartPointer<vtkCellLocator> locator;
	double source_center[3], target_center[3];

public:
	ICPAlgorithm();

	void setSource(vtkSmartPointer<vtkPolyData> source_);
	void setTarget(vtkSmartPointer<vtkPolyData> target_);
	void moveCenterOn();
	void moveCenterOff();
	void setMaxIter(int max_iter_);

	void registration();

	vtkMatrix4x4 * getTransformMatrix();

private:
	void move_center();
	void step_registration();
};