#include "mesh_processing.h"

#include <QtCore/QStringList>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QMessageBox>

#include <vtkAppendPolyData.h>
#include <vtkExtractEdges.h>
#include <vtkIdTypeArray.h>
#include <vtkMath.h>
#include <vtkLineSource.h>
#include <vtkLookupTable.h>
#include <vtkMatrix4x4.h>
#include <vtkNew.h>
#include <vtkOBJReader.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkScalarBarActor.h>
#include <vtkSphereSource.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

#include <iomanip>
#include <sstream>

#include "icp_algorithm.h"
#include "mesh_operation.h"
#include "vtkOFFReader.h"

MeshProcessing::MeshProcessing(QWidget *parent) : QMainWindow(parent) {
	ui.setupUi(this);

	this->mesh_processing_data_model_ = MeshProcessingDataModel::getInstance();
	this->color_table_reader_ = new ColorTableReader;
	this->vtk_widget_ = ui.vtk_widget;
	this->open_file_action_ = ui.open_file_action;
	this->read_color_table_action_ = ui.read_color_table_action;
	this->wireframe_mode_action_ = ui.wireframe_mode_action;
	this->observe_mode_action_ = ui.observe_mode_action;
	this->vertex_mode_action_ = ui.vertex_mode_action;
	this->face_mode_action_ = ui.face_mode_action;
	this->multi_vertex_mode_action_ = ui.multi_vertex_mode_action;
	this->display_normal_action_ = ui.display_normal_action;
	this->default_mode_action_ = ui.default_mode_action;
	this->discrete_mode_action_ = ui.discrete_mode_action;
	this->continuous_mode_action_ = ui.continuous_mode_action;
	this->icp_registration_action_ = ui.icp_registraion_action;
	this->fill_region_three_vertices_action_ = ui.fill_region_three_vertices_action;
	this->fill_region_two_vertices_action_ = ui.fill_region_two_vertices_action;
	this->list_widget_model_ = ui.list_widget_model;
	this->tab_widget_ = ui.tab_widget;
	this->max_iter_spin_box_ = ui.max_iter_spin_box;
	this->min_error_double_spin_box_ = ui.min_error_double_spin_box;
	this->move_center_radio_button_ = ui.move_center_radio_button;
	this->run_icp_button_ = ui.run_icp_button;
	this->iter_num_label_ = ui.iter_num_label;
	this->error_label_ = ui.error_label;
	this->matrix_label_ = ui.matrix_label;
	this->exit_icp_button_ = ui.exit_icp_button;
	this->cancel_icp_button_ = ui.cancel_icp_button;

	this->wireframe_mode_action_->setChecked(true);
	this->default_mode_action_->setEnabled(false);
	this->discrete_mode_action_->setEnabled(false);
	this->continuous_mode_action_->setEnabled(false);
	this->tab_widget_->setTabEnabled(1, false);

	connect(this->open_file_action_, SIGNAL(triggered()), this, SLOT(OnOpenFile()));
	connect(this->read_color_table_action_, SIGNAL(triggered()), this, SLOT(OnReadColorTable()));
	connect(this->wireframe_mode_action_, SIGNAL(triggered()), this, SLOT(OnWireframeMode()));
	connect(this->observe_mode_action_, SIGNAL(triggered()), this, SLOT(OnObserveMode()));
	connect(this->vertex_mode_action_, SIGNAL(triggered()), this, SLOT(OnVertexMode()));
	connect(this->face_mode_action_, SIGNAL(triggered()), this, SLOT(OnFaceMode()));
	connect(this->multi_vertex_mode_action_, SIGNAL(triggered()), this, SLOT(OnMultiVertexMode()));
	connect(this->display_normal_action_, SIGNAL(triggered()), this, SLOT(OnDisplayNormal()));
	connect(this->default_mode_action_, SIGNAL(triggered()), this, SLOT(OnDefaultMode()));
	connect(this->discrete_mode_action_, SIGNAL(triggered()), this, SLOT(OnDiscreteMode()));
	connect(this->continuous_mode_action_, SIGNAL(triggered()), this, SLOT(OnContinuousMode()));
	connect(this->icp_registration_action_, SIGNAL(triggered()), this, SLOT(OnICPRegistration()));
	connect(this->fill_region_three_vertices_action_, SIGNAL(triggered()), this, SLOT(OnFillRegionThreeVertices()));
	connect(this->fill_region_two_vertices_action_, SIGNAL(triggered()), this, SLOT(OnFillRegionTwoVertices()));
	connect(this->list_widget_model_, SIGNAL(itemChanged(QListWidgetItem *)), this, SLOT(OnListWidgetModelItemChanged(QListWidgetItem *)));
	connect(this->run_icp_button_, SIGNAL(clicked()), this, SLOT(OnRunICP()));
	connect(this->exit_icp_button_, SIGNAL(clicked()), this, SLOT(OnExitICP()));
	connect(this->cancel_icp_button_, SIGNAL(clicked()), this, SLOT(OnCancelICP()));

	connect(this->vtk_widget_->style, SIGNAL(selectVertex(vtkIdType)), this, SLOT(OnSelectVertex(vtkIdType)));
	connect(this->vtk_widget_->style, SIGNAL(selectFace(vtkIdType)), this, SLOT(OnSelectFace(vtkIdType)));
	connect(this->vtk_widget_->style, SIGNAL(selectMultiVertex(const std::vector<vtkIdType> &)), this, SLOT(OnSelectMultiVertex(const std::vector<vtkIdType> &)));
}

void MeshProcessing::OnOpenFile() {
	QStringList file_paths = QFileDialog::getOpenFileNames(this, QString::fromLocal8Bit("打开文件"), "./objs/", tr("OBJ files(*.obj);; OFF files(*.off)"));
	if (file_paths.size() == 0) return;

	this->default_mode_action_->setEnabled(false);
	this->discrete_mode_action_->setEnabled(false);
	this->continuous_mode_action_->setEnabled(false);

	this->removeMeshActors();
	this->removeVertexActors();
	this->removeFaceActors();
	this->removeMultiVertexActors();
	this->removeFillRegionFaceActors();
	this->mesh_processing_data_model_->mesh_vec_.clear();
	this->mesh_processing_data_model_->neighbor_vertex_actor_vec_.clear();
	this->mesh_processing_data_model_->actor_vec_.clear();
	this->mesh_processing_data_model_->highlight_vec_.clear();
	this->mesh_processing_data_model_->mean_edge_length_vec_.clear();
	this->list_widget_model_->clear();

	this->mesh_processing_data_model_->display_mode_ = MeshProcessingDataModel::DEFAULT;

	int number_of_points = 0, number_of_faces = 0, number_of_edges = 0;
	for (int i = 0; i < file_paths.size(); ++i) {
		QString file_name = file_paths[i];
		QString extension_name = file_name.split('/').back().split('.').back();

		vtkSmartPointer<vtkPolyData> mesh =
			vtkSmartPointer<vtkPolyData>::New();
		if (extension_name == QString::fromLocal8Bit("obj")) {
			vtkSmartPointer<vtkOBJReader> objReader =
				vtkSmartPointer<vtkOBJReader>::New();
			objReader->SetFileName(file_name.toLocal8Bit());
			objReader->Update();
			mesh->DeepCopy(objReader->GetOutput());
		} else {
			vtkSmartPointer<vtkOFFReader> offReader =
				vtkSmartPointer<vtkOFFReader>::New();
			offReader->SetFileName(file_name.toLocal8Bit());
			offReader->Update();
			mesh->DeepCopy(offReader->GetOutput());
		}

		QString only_name = file_name.split('/').back().split('.').front();
		QListWidgetItem * item = new QListWidgetItem(only_name);
		item->setCheckState(Qt::Checked);
		this->list_widget_model_->addItem(item);

		auto actor = this->vtk_widget_->addActor(mesh);
		this->vtk_widget_->highlightMesh(actor);

		auto wireframe_actor = this->vtk_widget_->addWireFrameActor(mesh);

		this->mesh_processing_data_model_->wireframe_actor_vec_.push_back(wireframe_actor);
		this->mesh_processing_data_model_->mesh_vec_.push_back(mesh);
		this->mesh_processing_data_model_->actor_vec_.push_back(actor);
		this->mesh_processing_data_model_->highlight_vec_.push_back(1);

		if (!this->wireframe_mode_action_->isChecked())
			this->mesh_processing_data_model_->wireframe_actor_vec_.back()->VisibilityOff();

		number_of_points += mesh->GetNumberOfPoints();
		number_of_faces += mesh->GetNumberOfCells();

		vtkSmartPointer<vtkExtractEdges> extractEdges =
			vtkSmartPointer<vtkExtractEdges>::New();
		extractEdges->SetInputData(mesh);
		extractEdges->Update();

		vtkSmartPointer<vtkPolyData> mesh_edge =
			vtkSmartPointer<vtkPolyData>::New();
		mesh_edge->DeepCopy(extractEdges->GetOutput());

		this->mesh_processing_data_model_->mesh_edge_vec_.push_back(mesh_edge);

		vtkSmartPointer<vtkIdList> pts =
			vtkSmartPointer<vtkIdList>::New();
		double mean_edge_length = 0.0;
		while (extractEdges->GetOutput()->GetLines()->GetNextCell(pts)) {
			double p0[3], p1[3];
			mesh_edge->GetPoint(pts->GetId(0), p0);
			mesh_edge->GetPoint(pts->GetId(1), p1);
			mean_edge_length += std::sqrt(vtkMath::Distance2BetweenPoints(p0, p1));
		}
		mean_edge_length /= mesh_edge->GetNumberOfLines();

		this->mesh_processing_data_model_->mean_edge_length_vec_.push_back(mean_edge_length);
		number_of_edges += mesh_edge->GetNumberOfLines();
	}

	this->resetParameters();

	this->vtk_widget_->updateBottomText(number_of_points, number_of_faces, number_of_edges);
	this->vtk_widget_->resetCamera();
}

void MeshProcessing::OnReadColorTable() {
	if (this->mesh_processing_data_model_->mesh_vec_.size() != 1) {
		QMessageBox::critical(this, QString::fromLocal8Bit("错误"), QString::fromLocal8Bit("请检查目前读入的网格数是否为1！"));
		return;
	}

	QString file_path = QFileDialog::getOpenFileName(this, QString::fromLocal8Bit("打开文件"), "./", tr("txt Files(*.txt)"));
	if (file_path.size() == 0) return;

	this->color_table_reader_->setColorTableName(file_path);
	this->color_table_reader_->setMesh(this->mesh_processing_data_model_->mesh_vec_[0]);
	if (this->color_table_reader_->read() == false) {
		QMessageBox::critical(this, QString::fromLocal8Bit("错误"), QString::fromLocal8Bit("请检查颜色表文件是否正确！"));
		return;
	}

	this->default_mode_action_->setEnabled(true);
	this->discrete_mode_action_->setEnabled(true);
	this->continuous_mode_action_->setEnabled(true);

	QMessageBox::information(this, QString::fromLocal8Bit("提示"), QString::fromLocal8Bit("颜色表读取完毕！请选择着色模式！"));
	return;
}

void MeshProcessing::OnWireframeMode() {
	bool isChecked = this->wireframe_mode_action_->isChecked();
	for (auto actor : this->mesh_processing_data_model_->wireframe_actor_vec_)
		actor->SetVisibility(isChecked);
	this->vtk_widget_->update();
}

void MeshProcessing::OnObserveMode() {
	this->removeVertexActors();
	this->removeFaceActors();
	this->removeMultiVertexActors();
	this->removeFillRegionFaceActors();
	this->mesh_processing_data_model_->pick_mode_ = MeshProcessingDataModel::OBSERVE;
	this->vtk_widget_->updateTopText();
}

void MeshProcessing::OnVertexMode() {
	this->removeFaceActors();
	this->removeMultiVertexActors();
	this->removeFillRegionFaceActors();
	this->mesh_processing_data_model_->pick_mode_ = MeshProcessingDataModel::VERTEX;
	this->vtk_widget_->updateTopText();
}

void MeshProcessing::OnFaceMode() {
	this->removeVertexActors();
	this->removeMultiVertexActors();
	this->removeFillRegionFaceActors();
	this->mesh_processing_data_model_->pick_mode_ = MeshProcessingDataModel::FACE;
	this->vtk_widget_->updateTopText();
}

void MeshProcessing::OnMultiVertexMode() {
	this->removeVertexActors();
	this->removeFaceActors();
	this->removeFillRegionFaceActors();
	this->mesh_processing_data_model_->pick_mode_ = MeshProcessingDataModel::MULTI_VERTEX;
	this->vtk_widget_->updateTopText();
}

void MeshProcessing::OnDisplayNormal() {
	if (this->mesh_processing_data_model_->pick_mode_ != MeshProcessingDataModel::FACE ||
		this->mesh_processing_data_model_->selected_face_id_ == -1) {
		QMessageBox::warning(this, QString::fromLocal8Bit("警告"), QString::fromLocal8Bit("请在面拾取模式中选取一个面片！"));
		return;
	}

	int cellId = this->mesh_processing_data_model_->selected_face_id_;
	vtkCell * cell = this->mesh_processing_data_model_->combined_mesh_->GetCell(cellId);

	double p0[3], p1[3], p2[3];
	this->mesh_processing_data_model_->combined_mesh_->GetPoint(cell->GetPointId(0), p0);
	this->mesh_processing_data_model_->combined_mesh_->GetPoint(cell->GetPointId(1), p1);
	this->mesh_processing_data_model_->combined_mesh_->GetPoint(cell->GetPointId(2), p2);

	double c[3];
	for (int i = 0; i < 3; ++i) c[i] = (p0[i] + p1[i] + p2[i]) / 3;

	double p0p1[3], p0p2[3];
	vtkMath::Subtract(p1, p0, p0p1);
	vtkMath::Subtract(p2, p0, p0p2);

	double normal[3];
	vtkMath::Cross(p0p1, p0p2, normal);
	vtkMath::Normalize(normal);

	double constant = (
		std::sqrt(vtkMath::Distance2BetweenPoints(p0, p1)) +
		std::sqrt(vtkMath::Distance2BetweenPoints(p1, p2)) +
		std::sqrt(vtkMath::Distance2BetweenPoints(p2, p0))
	);
	double e[3];
	for (int i = 0; i < 3; ++i) e[i] = c[i] + constant * normal[i];

	vtkSmartPointer<vtkLineSource> lineSource =
		vtkSmartPointer<vtkLineSource>::New();
	lineSource->SetPoint1(c);
	lineSource->SetPoint2(e);
	lineSource->Update();

	this->mesh_processing_data_model_->selected_face_normal_actor_ = this->vtk_widget_->addActor(lineSource->GetOutput());
	this->mesh_processing_data_model_->selected_face_normal_actor_->GetProperty()->SetColor(.8, .8, .2);
	this->mesh_processing_data_model_->selected_face_normal_actor_->GetProperty()->SetLineWidth(5);

	this->vtk_widget_->update();
}

void MeshProcessing::OnDefaultMode() {
	this->mesh_processing_data_model_->display_mode_ = MeshProcessingDataModel::DEFAULT;
	if (this->mesh_processing_data_model_->highlight_vec_[0] == 1)
		this->vtk_widget_->highlightMesh(this->mesh_processing_data_model_->actor_vec_[0]);
	else
		this->vtk_widget_->unhighlightMesh(this->mesh_processing_data_model_->actor_vec_[0]);
	this->vtk_widget_->update();
}

void MeshProcessing::OnDiscreteMode() {
	if (this->mesh_processing_data_model_->actor_vec_.size() != 1) {
		QMessageBox::critical(this, QString::fromLocal8Bit("错误"), QString::fromLocal8Bit("请检查目前读入的网格数是否为1！"));
		return;
	}

	this->mesh_processing_data_model_->display_mode_ = MeshProcessingDataModel::DISCRETE;
	this->mesh_processing_data_model_->mesh_vec_[0] = this->color_table_reader_->turnToDiscrete();

	this->mesh_processing_data_model_->hueLut = vtkSmartPointer<vtkLookupTable>::New();
	this->mesh_processing_data_model_->hueLut->SetNumberOfTableValues(this->color_table_reader_->maxScalar() - this->color_table_reader_->minScalar() + 1);
	this->mesh_processing_data_model_->hueLut->Build();

	for (int i = 0; i < this->color_table_reader_->maxScalar() - this->color_table_reader_->minScalar() + 1; ++i) {
		double hue = i * 1.0 / (this->color_table_reader_->maxScalar() - this->color_table_reader_->minScalar() + 1);
		double r, g, b;
		vtkMath::HSVToRGB(hue, 1.0, 1.0, &r, &g, &b);
		this->mesh_processing_data_model_->hueLut->SetTableValue(i, r, g, b, 1);
	}

	this->mesh_processing_data_model_->actor_vec_[0]->GetMapper()->SetScalarRange(0, this->color_table_reader_->maxScalar() - this->color_table_reader_->minScalar() + 1);
	this->mesh_processing_data_model_->actor_vec_[0]->GetMapper()->SetScalarModeToUseCellData();
	this->mesh_processing_data_model_->actor_vec_[0]->GetMapper()->SetColorModeToMapScalars();
	this->mesh_processing_data_model_->actor_vec_[0]->GetMapper()->SetLookupTable(this->mesh_processing_data_model_->hueLut);

	if (this->mesh_processing_data_model_->highlight_vec_[0] == 1)
		this->vtk_widget_->highlightMesh(this->mesh_processing_data_model_->actor_vec_[0]);
	else
		this->vtk_widget_->unhighlightMesh(this->mesh_processing_data_model_->actor_vec_[0]);
	this->vtk_widget_->update();
}

void MeshProcessing::OnContinuousMode() {
	if (this->mesh_processing_data_model_->actor_vec_.size() != 1) {
		QMessageBox::critical(this, QString::fromLocal8Bit("错误"), QString::fromLocal8Bit("请检查目前读入的网格数是否为1！"));
		return;
	}

	this->mesh_processing_data_model_->display_mode_ = MeshProcessingDataModel::CONTINUOUS;
	this->mesh_processing_data_model_->mesh_vec_[0] = this->color_table_reader_->turnToContinuous();

	this->mesh_processing_data_model_->hueLut = vtkSmartPointer<vtkLookupTable>::New();
	this->mesh_processing_data_model_->hueLut->SetTableRange(this->color_table_reader_->minScalar(), this->color_table_reader_->maxScalar() + 1);
	this->mesh_processing_data_model_->hueLut->SetHueRange(0, 1);
	this->mesh_processing_data_model_->hueLut->SetSaturationRange(1, 1);
	this->mesh_processing_data_model_->hueLut->SetValueRange(1, 1);
	this->mesh_processing_data_model_->hueLut->Build();

	this->mesh_processing_data_model_->actor_vec_[0]->GetMapper()->SetScalarRange(this->color_table_reader_->minScalar(), this->color_table_reader_->maxScalar() + 1);
	this->mesh_processing_data_model_->actor_vec_[0]->GetMapper()->SetScalarModeToUseCellData();
	this->mesh_processing_data_model_->actor_vec_[0]->GetMapper()->SetColorModeToMapScalars();
	this->mesh_processing_data_model_->actor_vec_[0]->GetMapper()->SetLookupTable(this->mesh_processing_data_model_->hueLut);

	if (this->mesh_processing_data_model_->highlight_vec_[0] == 1)
		this->vtk_widget_->highlightMesh(this->mesh_processing_data_model_->actor_vec_[0]);
	else
		this->vtk_widget_->unhighlightMesh(this->mesh_processing_data_model_->actor_vec_[0]);
	this->vtk_widget_->update();
}

void MeshProcessing::OnICPRegistration() {
	std::vector<int> active_ids;
	for (int i = 0; i < this->mesh_processing_data_model_->highlight_vec_.size(); ++i) {
		if (this->mesh_processing_data_model_->highlight_vec_[i])
			active_ids.push_back(i);
	}

	if (active_ids.size() != 2) {
		QMessageBox::warning(this, QString::fromLocal8Bit("警告"), QString::fromLocal8Bit("请检查目前选中的网格数是否为2！"));
		return;
	}

	this->wireframe_mode_action_->setChecked(false);
	OnWireframeMode();
	this->disableAllActions();
	this->OnObserveMode();

	this->mesh_processing_data_model_->source_id = active_ids.front();
	this->mesh_processing_data_model_->target_id = active_ids.back();

	this->mesh_processing_data_model_->actor_vec_[this->mesh_processing_data_model_->target_id]->GetProperty()->SetColor(.8, .2, .2);
	this->vtk_widget_->update();

	this->tab_widget_->setTabEnabled(0, false);
	this->tab_widget_->setTabEnabled(1, true);
	this->tab_widget_->setCurrentIndex(1);
}

void MeshProcessing::OnFillRegionThreeVertices() {
	this->removeFillRegionFaceActors();

	for (int i = 0; i < this->mesh_processing_data_model_->combined_mesh_->GetNumberOfCells(); ++i) {
		vtkSmartPointer<vtkIdList> ptIds =
			vtkSmartPointer<vtkIdList>::New();
		this->mesh_processing_data_model_->combined_mesh_->GetCellPoints(i, ptIds);

		int validNum = 0;
		for (int k = 0; k < 3; ++k)
			validNum += this->mesh_processing_data_model_->selected_multi_vertex_set_.find(ptIds->GetId(k)) != this->mesh_processing_data_model_->selected_multi_vertex_set_.end() ? 1 : 0;

		if (validNum == 3) {
			this->mesh_processing_data_model_->fill_region_face_actor_vec_.push_back(this->vtk_widget_->highlightFace(this->mesh_processing_data_model_->combined_mesh_, i));
			this->mesh_processing_data_model_->fill_region_face_actor_vec_.back()->GetProperty()->SetColor(.8, .8, .2);
		}
	}

	this->vtk_widget_->update();
}

void MeshProcessing::OnFillRegionTwoVertices() {
	this->removeFillRegionFaceActors();

	for (int i = 0; i < this->mesh_processing_data_model_->combined_mesh_->GetNumberOfCells(); ++i) {
		vtkSmartPointer<vtkIdList> ptIds =
			vtkSmartPointer<vtkIdList>::New();
		this->mesh_processing_data_model_->combined_mesh_->GetCellPoints(i, ptIds);

		int validNum = 0;
		for (int k = 0; k < 3; ++k)
			validNum += this->mesh_processing_data_model_->selected_multi_vertex_set_.find(ptIds->GetId(k)) != this->mesh_processing_data_model_->selected_multi_vertex_set_.end() ? 1 : 0;

		if (validNum >= 2) {
			this->mesh_processing_data_model_->fill_region_face_actor_vec_.push_back(this->vtk_widget_->highlightFace(this->mesh_processing_data_model_->combined_mesh_, i));
			this->mesh_processing_data_model_->fill_region_face_actor_vec_.back()->GetProperty()->SetColor(.8, .8, .2);
		}
	}

	this->vtk_widget_->update();
}

void MeshProcessing::OnListWidgetModelItemChanged(QListWidgetItem * item) {
	int row_id = this->list_widget_model_->row(item);

	if (item->checkState() == Qt::Unchecked) {
		this->vtk_widget_->unhighlightMesh(this->mesh_processing_data_model_->actor_vec_[row_id]);
		this->mesh_processing_data_model_->highlight_vec_[row_id] = 0;
	} else {
		this->vtk_widget_->highlightMesh(this->mesh_processing_data_model_->actor_vec_[row_id]);
		this->mesh_processing_data_model_->highlight_vec_[row_id] = 1;
	}

	this->resetParameters();

	this->vtk_widget_->update();
}

void MeshProcessing::OnRunICP() {
	ICPAlgorithm icp;
	icp.setSource(this->mesh_processing_data_model_->mesh_vec_[this->mesh_processing_data_model_->source_id]);
	icp.setTarget(this->mesh_processing_data_model_->mesh_vec_[this->mesh_processing_data_model_->target_id]);
	if (this->move_center_radio_button_->isChecked())
		icp.moveCenterOn();
	else
		icp.moveCenterOff();
	icp.setMaxIter(this->max_iter_spin_box_->value());
	icp.setMinError(this->min_error_double_spin_box_->value());
	icp.registration();

	vtkMatrix4x4 * transform_matrix = icp.getTransformMatrix();
	
	std::ostringstream oss;
	oss << setprecision(2);
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) oss << transform_matrix->GetElement(i, j) << "\t";
		oss << std::endl;
	}

	this->iter_num_label_->setText(QString::number(icp.getIterNum()));
	this->error_label_->setText(QString::number(icp.getError()));
	this->matrix_label_->setText(QString::fromStdString(oss.str()));

	vtkSmartPointer<vtkTransform> transform =
		vtkSmartPointer<vtkTransform>::New();
	transform->SetMatrix(transform_matrix);
	transform->Update();

	vtkNew<vtkTransformPolyDataFilter> transformFilter;
	transformFilter->SetInputData(this->mesh_processing_data_model_->mesh_vec_[0]);
	transformFilter->SetTransform(transform);
	transformFilter->Update();

	//this->mesh_processing_data_model_->mesh_vec_[0]->DeepCopy(transformFilter->GetOutput());

	vtkPolyDataMapper * mapper = vtkPolyDataMapper::SafeDownCast(this->mesh_processing_data_model_->actor_vec_[0]->GetMapper());
	mapper->SetInputData(transformFilter->GetOutput());

	this->vtk_widget_->resetCamera();
	this->vtk_widget_->update();
}

void MeshProcessing::OnExitICP() {
	this->tab_widget_->setTabEnabled(0, true);
	this->tab_widget_->setTabEnabled(1, false);
	this->tab_widget_->setCurrentIndex(0);

	this->iter_num_label_->setText("");
	this->error_label_->setText("");
	this->matrix_label_->setText("");

	this->enableAllActions();

	this->mesh_processing_data_model_->mesh_vec_[0]->DeepCopy(this->mesh_processing_data_model_->actor_vec_[0]->GetMapper()->GetInput());
	this->mesh_processing_data_model_->actor_vec_[this->mesh_processing_data_model_->target_id]->GetProperty()->SetColor(.0, .5, 1.);
	this->vtk_widget_->update();

	this->resetParameters();
}

void MeshProcessing::OnCancelICP() {
	this->tab_widget_->setTabEnabled(0, true);
	this->tab_widget_->setTabEnabled(1, false);
	this->tab_widget_->setCurrentIndex(0);

	this->iter_num_label_->setText("");
	this->error_label_->setText("");
	this->matrix_label_->setText("");

	this->enableAllActions();

	vtkPolyDataMapper * mapper = vtkPolyDataMapper::SafeDownCast(this->mesh_processing_data_model_->actor_vec_[0]->GetMapper());
	mapper->SetInputData(this->mesh_processing_data_model_->mesh_vec_[0]);

	this->mesh_processing_data_model_->actor_vec_[this->mesh_processing_data_model_->target_id]->GetProperty()->SetColor(.0, .5, 1.);
	this->vtk_widget_->update();

	this->resetParameters();
}

void MeshProcessing::OnSelectVertex(vtkIdType id) {
	if (id == -1) return;

	this->removeVertexActors();

	this->mesh_processing_data_model_->selected_vertex_actor_ = this->vtk_widget_->highlightVertex(this->mesh_processing_data_model_->combined_mesh_, id);
	this->mesh_processing_data_model_->selected_vertex_actor_->GetProperty()->SetInterpolationToGouraud();
	this->mesh_processing_data_model_->selected_vertex_actor_->GetProperty()->SetColor(.8, .2, .2);

	auto neighborVertexIds = MeshOperation::getConnectedVertices(this->mesh_processing_data_model_->combined_mesh_, id);
	for (const auto & neighbor : neighborVertexIds) {
		this->mesh_processing_data_model_->neighbor_vertex_actor_vec_.push_back(this->vtk_widget_->highlightVertex(this->mesh_processing_data_model_->combined_mesh_, neighbor));
		this->mesh_processing_data_model_->neighbor_vertex_actor_vec_.back()->GetProperty()->SetInterpolationToGouraud();
		this->mesh_processing_data_model_->neighbor_vertex_actor_vec_.back()->GetProperty()->SetColor(.8, .8, .2);
	}

	auto neighborFaceIds = MeshOperation::getVertexConnectedFaces(this->mesh_processing_data_model_->combined_mesh_, id);
	for (const auto & neighbor : neighborFaceIds) {
		this->mesh_processing_data_model_->neighbor_face_actor_vec_.push_back(this->vtk_widget_->highlightFace(this->mesh_processing_data_model_->combined_mesh_, neighbor));
		this->mesh_processing_data_model_->neighbor_face_actor_vec_.back()->GetProperty()->SetColor(.8, .5, .2);
	}

	this->vtk_widget_->update();
}

void MeshProcessing::OnSelectFace(vtkIdType id) {
	if (id == -1) return;

	this->removeFaceActors();

	this->mesh_processing_data_model_->selected_face_actor_ = this->vtk_widget_->highlightFace(this->mesh_processing_data_model_->combined_mesh_, id);
	this->mesh_processing_data_model_->selected_face_actor_->GetProperty()->SetColor(.8, .2, .2);

	auto neighborFaceIds = MeshOperation::getFaceConnectedFaces(this->mesh_processing_data_model_->combined_mesh_, id);
	for (const auto & neighbor : neighborFaceIds) {
		this->mesh_processing_data_model_->neighbor_face2_actor_vec_.push_back(this->vtk_widget_->highlightFace(this->mesh_processing_data_model_->combined_mesh_, neighbor));
		this->mesh_processing_data_model_->neighbor_face2_actor_vec_.back()->GetProperty()->SetColor(.8, .8, .2);
	}

	this->mesh_processing_data_model_->selected_face_id_ = id;

	this->vtk_widget_->update();
}

void MeshProcessing::OnSelectMultiVertex(const std::vector<vtkIdType>& ids) {
	this->removeMultiVertexActors();

	this->mesh_processing_data_model_->selected_multi_vertex_set_.clear();
	for (const auto & id : ids) {
		this->mesh_processing_data_model_->selected_multi_vertex_actor_vec_.push_back(this->vtk_widget_->highlightVertex(this->mesh_processing_data_model_->combined_mesh_, id));
		this->mesh_processing_data_model_->selected_multi_vertex_actor_vec_.back()->GetProperty()->SetInterpolationToGouraud();
		this->mesh_processing_data_model_->selected_multi_vertex_actor_vec_.back()->GetProperty()->SetColor(.8, .2, .2);

		this->mesh_processing_data_model_->selected_multi_vertex_set_.insert(id);
	}

	this->vtk_widget_->update();
}

void MeshProcessing::resetParameters() {
	vtkSmartPointer<vtkAppendPolyData> appendFilter =
		vtkSmartPointer<vtkAppendPolyData>::New();
	int cnt = 0;
	int num_edges = 0;
	for (int i = 0; i < this->mesh_processing_data_model_->highlight_vec_.size(); ++i) {
		if (this->mesh_processing_data_model_->highlight_vec_[i]) {
			appendFilter->AddInputData(this->mesh_processing_data_model_->mesh_vec_[i]);
			num_edges += this->mesh_processing_data_model_->mesh_edge_vec_[i]->GetNumberOfLines();
			++cnt;
		}
	}
	if (cnt > 0) {
		appendFilter->Update();
		this->mesh_processing_data_model_->combined_mesh_ =
			vtkSmartPointer<vtkPolyData>::New();
		this->mesh_processing_data_model_->combined_mesh_->DeepCopy(appendFilter->GetOutput());

		vtkSmartPointer<vtkIdTypeArray> numberScalarArray =
			vtkSmartPointer<vtkIdTypeArray>::New();
		numberScalarArray->SetNumberOfComponents(1);
		numberScalarArray->SetName("number");
		numberScalarArray->SetNumberOfValues(this->mesh_processing_data_model_->combined_mesh_->GetNumberOfPoints());
		for (int i = 0; i < this->mesh_processing_data_model_->combined_mesh_->GetNumberOfPoints(); ++i)
			numberScalarArray->SetValue(i, i);

		this->mesh_processing_data_model_->combined_mesh_->GetPointData()->AddArray(numberScalarArray);
	} else this->mesh_processing_data_model_->combined_mesh_ = nullptr;

	int edge_count = 0;
	this->mesh_processing_data_model_->mean_edge_length = 0.0;
	for (int i = 0; i < this->mesh_processing_data_model_->highlight_vec_.size(); ++i) {
		if (this->mesh_processing_data_model_->highlight_vec_[i]) {
			this->mesh_processing_data_model_->mean_edge_length += 
				this->mesh_processing_data_model_->mesh_edge_vec_[i]->GetNumberOfLines() * this->mesh_processing_data_model_->mean_edge_length_vec_[i];
			edge_count += this->mesh_processing_data_model_->mesh_edge_vec_[i]->GetNumberOfLines();
		}
	}
	this->mesh_processing_data_model_->mean_edge_length /= edge_count;

	this->vtk_widget_->updateBottomText(
		this->mesh_processing_data_model_->combined_mesh_->GetNumberOfPoints(),
		this->mesh_processing_data_model_->combined_mesh_->GetNumberOfCells(),
		num_edges
	);
}

void MeshProcessing::removeMeshActors() {
	for (const auto & actor : this->mesh_processing_data_model_->actor_vec_)
		this->vtk_widget_->removeActor(actor);
	for (const auto & actor : this->mesh_processing_data_model_->wireframe_actor_vec_)
		this->vtk_widget_->removeActor(actor);
}

void MeshProcessing::removeVertexActors() {
	this->vtk_widget_->removeActor(this->mesh_processing_data_model_->selected_vertex_actor_);
	for (const auto & actor : this->mesh_processing_data_model_->neighbor_vertex_actor_vec_)
		this->vtk_widget_->removeActor(actor);
	for (const auto & actor : this->mesh_processing_data_model_->neighbor_face_actor_vec_)
		this->vtk_widget_->removeActor(actor);
}

void MeshProcessing::removeFaceActors() {
	this->vtk_widget_->removeActor(this->mesh_processing_data_model_->selected_face_actor_);
	for (const auto & actor : this->mesh_processing_data_model_->neighbor_face2_actor_vec_)
		this->vtk_widget_->removeActor(actor);
	this->vtk_widget_->removeActor(this->mesh_processing_data_model_->selected_face_normal_actor_);
}

void MeshProcessing::removeMultiVertexActors() {
	for (const auto & actor : this->mesh_processing_data_model_->selected_multi_vertex_actor_vec_)
		this->vtk_widget_->removeActor(actor);
}

void MeshProcessing::removeFillRegionFaceActors() {
	for (const auto & actor : this->mesh_processing_data_model_->fill_region_face_actor_vec_)
		this->vtk_widget_->removeActor(actor);
}

void MeshProcessing::enableAllActions() {
	this->open_file_action_->setEnabled(true);
	this->read_color_table_action_->setEnabled(true);
	this->wireframe_mode_action_->setEnabled(true);
	this->observe_mode_action_->setEnabled(true);
	this->vertex_mode_action_->setEnabled(true);
	this->face_mode_action_->setEnabled(true);
	this->multi_vertex_mode_action_->setEnabled(true);
	this->display_normal_action_->setEnabled(true);
	this->icp_registration_action_->setEnabled(true);
	this->fill_region_three_vertices_action_->setEnabled(true);
	this->fill_region_two_vertices_action_->setEnabled(true);
}

void MeshProcessing::disableAllActions() {
	this->open_file_action_->setEnabled(false);
	this->read_color_table_action_->setEnabled(false);
	this->wireframe_mode_action_->setEnabled(false);
	this->observe_mode_action_->setEnabled(false);
	this->vertex_mode_action_->setEnabled(false);
	this->face_mode_action_->setEnabled(false);
	this->multi_vertex_mode_action_->setEnabled(false);
	this->display_normal_action_->setEnabled(false);
	this->icp_registration_action_->setEnabled(false);
	this->fill_region_three_vertices_action_->setEnabled(false);
	this->fill_region_two_vertices_action_->setEnabled(false);
}
