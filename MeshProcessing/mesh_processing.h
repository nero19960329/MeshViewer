#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_mesh_processing.h"

#include "vtk_widget.h"

class MeshProcessing : public QMainWindow {
	Q_OBJECT

public:
	MeshProcessing(QWidget *parent = nullptr);

private:
	Ui::MeshProcessingClass ui;

	VTKWidget * vtk_widget_;
};
