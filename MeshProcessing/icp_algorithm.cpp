#include "icp_algorithm.h"

#include <vtkCenterOfMass.h>
#include <vtkMath.h>
#include <vtkMatrix3x3.h>
#include <vtkMatrix4x4.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkTransformPolyDataFilter.h>

ICPAlgorithm::ICPAlgorithm() {
	this->source = nullptr;
	this->target = nullptr;
	this->is_move_center = false;
	this->max_iter = 0;
	this->min_error = 0.0;
	this->transform = vtkSmartPointer<vtkTransform>::New();
	this->locator = vtkSmartPointer<vtkCellLocator>::New();
}

void ICPAlgorithm::setSource(vtkSmartPointer<vtkPolyData> source_) {
	this->source = vtkSmartPointer<vtkPolyData>::New();
	this->source->DeepCopy(source_);
}

void ICPAlgorithm::setTarget(vtkSmartPointer<vtkPolyData> target_) {
	this->target = vtkSmartPointer<vtkPolyData>::New();
	this->target->DeepCopy(target_);
}

void ICPAlgorithm::moveCenterOn() {
	this->is_move_center = true;
}

void ICPAlgorithm::moveCenterOff() {
	this->is_move_center = false;
}

void ICPAlgorithm::setMaxIter(int max_iter_) {
	this->max_iter = max_iter_;
}

void ICPAlgorithm::setMinError(double min_error_) {
	this->min_error = min_error_;
}

void ICPAlgorithm::registration() {
	this->transform->Identity();
	this->transform->PostMultiply();

	vtkNew<vtkCenterOfMass> centerOfMass;
	centerOfMass->SetInputData(this->source);
	centerOfMass->Update();
	centerOfMass->GetCenter(this->source_center);
	centerOfMass->SetInputData(this->target);
	centerOfMass->Update();
	centerOfMass->GetCenter(this->target_center);

	if (this->is_move_center) {
		this->transform->Translate(-source_center[0], -source_center[1], -source_center[2]);
		this->transform->Update();
		move_center();
	}

	this->locator->SetDataSet(this->target);
	this->locator->SetNumberOfCellsPerBucket(1);
	this->locator->BuildLocator();

	for (iter_num = 0; iter_num < this->max_iter; ++iter_num) {
		step_registration();
		if (this->error < this->min_error)
			break;
	}

	if (this->is_move_center) {
		this->transform->Translate(target_center[0], target_center[1], target_center[2]);
		this->transform->Update();
	}
}

vtkMatrix4x4 * ICPAlgorithm::getTransformMatrix() {
	return this->transform->GetMatrix();
}

int ICPAlgorithm::getIterNum() {
	return this->iter_num;
}

double ICPAlgorithm::getError() {
	return this->error;
}

void ICPAlgorithm::move_center() {
	vtkSmartPointer<vtkTransform> center_transform =
		vtkSmartPointer<vtkTransform>::New();
	center_transform->Translate(-source_center[0], -source_center[1], -source_center[2]);
	center_transform->Update();
	vtkNew<vtkTransformPolyDataFilter> transformFilter;
	transformFilter->SetInputData(this->source);
	transformFilter->SetTransform(center_transform);
	transformFilter->Update();
	this->source->DeepCopy(transformFilter->GetOutput());

	center_transform->Identity();
	center_transform->Translate(-target_center[0], -target_center[1], -target_center[2]);
	center_transform->Update();
	transformFilter->SetInputData(this->target);
	transformFilter->SetTransform(center_transform);
	transformFilter->Update();
	this->target->DeepCopy(transformFilter->GetOutput());
}

void ICPAlgorithm::step_registration() {
	vtkSmartPointer<vtkPoints> sourcePoints = this->source->GetPoints();
	vtkSmartPointer<vtkPoints> targetPoints =
		vtkSmartPointer<vtkPoints>::New();

	int n = this->source->GetNumberOfPoints();
	for (int i = 0; i < n; ++i) {
		double closestPoint[3];
		vtkIdType cellId;
		int subId;
		double dist2;
		this->locator->FindClosestPoint(this->source->GetPoint(i), closestPoint, cellId, subId, dist2);
		targetPoints->InsertNextPoint(closestPoint);
	}

	double M[3][3] = { 0 };
	for (int k = 0; k < n; ++k) {
		double p[3], q[3];
		sourcePoints->GetPoint(k, p);
		targetPoints->GetPoint(k, q);

		for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j)
			M[i][j] += p[i] * q[j];
	}

	double U[3][3], W[3], VT[3][3];
	vtkMath::SingularValueDecomposition3x3(M, U, W, VT);

	double R[3][3] = { 0 };
	for (int k = 0; k < 3; ++k) {
		for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j)
			R[i][j] += U[i][k] * VT[k][j];
	}

	vtkMatrix4x4 * R_mat = vtkMatrix4x4::New();
	R_mat->Zero();
	for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j)
		R_mat->SetElement(i, j, R[j][i]);
	R_mat->SetElement(3, 3, 1.0);
	
	this->transform->Concatenate(R_mat);
	this->transform->Update();

	vtkSmartPointer<vtkTransform> step_transform =
		vtkSmartPointer<vtkTransform>::New();
	step_transform->SetMatrix(R_mat);
	step_transform->Update();
	vtkNew<vtkTransformPolyDataFilter> transformFilter;
	transformFilter->SetInputData(this->source);
	transformFilter->SetTransform(step_transform);
	transformFilter->Update();
	this->source->DeepCopy(transformFilter->GetOutput());

	sourcePoints = this->source->GetPoints();
	this->error = 0.0;
	for (int i = 0; i < n; ++i) {
		double p[3], q[3];
		sourcePoints->GetPoint(i, p);
		targetPoints->GetPoint(i, q);
		this->error += vtkMath::Distance2BetweenPoints(p, q);
	}
	this->error = std::sqrt(this->error / n);

	R_mat->Delete();
}
