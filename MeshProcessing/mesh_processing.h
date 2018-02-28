#pragma once

#include <QtWidgets/QAction>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QScrollArea.h>
#include <QtWidgets/QSlider>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QWidget>

#include "color_table_handler.h"
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
	ColorTableHandler * color_table_handler_;

	VTKWidget * vtk_widget_;

	QAction * open_file_action_;
	QAction * read_color_table_action_;
	QAction * write_color_table_action_;
	QAction * wireframe_mode_action_;
	QAction * observe_mode_action_;
	QAction * vertex_mode_action_;
	QAction * face_mode_action_;
	QAction * multi_vertex_mode_action_;
	QAction * display_normal_action_;
	QAction * default_mode_action_;
	QAction * discrete_mode_action_;
	QAction * continuous_mode_action_;
	QAction * icp_registration_action_;
	QAction * segment_action_;
	QAction * fill_region_three_vertices_action_;
	QAction * fill_region_two_vertices_action_;

	QListWidget * list_widget_model_;

	QScrollArea * register_scroll_area_;
	QScrollArea * segment_scroll_area_;

	QTabWidget * tab_widget_;
	QSpinBox * max_iter_spin_box_;
	QDoubleSpinBox * min_error_double_spin_box_;
	QRadioButton * move_center_radio_button_;
	QPushButton * run_icp_button_;
	QLabel * iter_num_label_;
	QLabel * error_label_;
	QLabel * matrix_label_;
	QPushButton * exit_icp_button_;
	QPushButton * cancel_icp_button_;

	QSpinBox * seed_num_spin_box_;
	QDoubleSpinBox * phy_ratio_double_spin_box_;
	QSlider * cluster_num_slider_;
	QPushButton * run_segment_button_;
	QPushButton * exit_segment_button_;
	QPushButton * cancel_segment_button_;

public slots:
	void OnOpenFile();
	void OnReadColorTable();
	void OnWriteColorTable();
	void OnWireframeMode();
	void OnObserveMode();
	void OnVertexMode();
	void OnFaceMode();
	void OnMultiVertexMode();
	void OnDisplayNormal();
	void OnDefaultMode();
	void OnDiscreteMode();
	void OnContinuousMode();
	void OnICPRegistration();
	void OnSegment();
	void OnFillRegionThreeVertices();
	void OnFillRegionTwoVertices();
	void OnListWidgetModelItemChanged(QListWidgetItem *);
	void OnRunICP();
	void OnExitICP();
	void OnCancelICP();
	void OnClusterNumChanged(int n);
	void OnRunSegment();
	void OnExitSegment();
	void OnCancelSegment();

	void OnSelectVertex(vtkIdType id);
	void OnSelectFace(vtkIdType id);
	void OnSelectMultiVertex(const std::vector<vtkIdType> & ids);

private:
	void resetParameters();
	void removeMeshActors();
	void removeVertexActors();
	void removeFaceActors();
	void removeMultiVertexActors();
	void removeFillRegionFaceActors();
	void enableAllActions();
	void disableAllActions();
};
