#include "mesh_processing.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[]) {
	QApplication a(argc, argv);
	MeshProcessing w;
	w.show();
	return a.exec();
}
