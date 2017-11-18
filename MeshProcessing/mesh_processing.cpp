#include "mesh_processing.h"

MeshProcessing::MeshProcessing(QWidget *parent) : QMainWindow(parent) {
	ui.setupUi(this);

	this->vtk_widget_ = ui.vtk_widget;
}
