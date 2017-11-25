#pragma once

#include <QtWidgets/QAction>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QPushButton>

#include "mesh_processing_data_model.h"
#include "ui_mesh_processing.h"
#include "vtk_widget.h"

class MeshProcessing : public QMainWindow {
	Q_OBJECT

public:
	MeshProcessing(QWidget *parent = nullptr);

private:
	Ui::MeshProcessingClass ui;

	MeshProcessingDataModel * mesh_processing_data_model_;

	VTKWidget * vtk_widget_;

	QAction * open_file_action_;
	QAction * observe_mode_action_;
	QAction * vertex_mode_action_;
	QAction * face_mode_action_;

	QListWidget * list_widget_model_;

public slots:
	void OnOpenFile();
	void OnListWidgetModelItemChanged(QListWidgetItem *);
	void OnObserveMode();
	void OnVertexMode();
	void OnFaceMode();

	void OnSelectVertex(vtkIdType id);
	void OnSelectFace(vtkIdType id);

private:
	void resetParameters();
	void removeMeshActors();
	void removeVertexActors();
	void removeFaceActors();
};
