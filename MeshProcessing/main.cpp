#include "mesh_processing.h"
#include <QtWidgets/QApplication>

#include "fibonacci_heap.h"

int main(int argc, char *argv[]) {
	QApplication a(argc, argv);
	MeshProcessing w;
	w.show();
	return a.exec();

	/*FibHeap<double> fibHeap;
	for (int i = 10; i >= 1; --i)
		fibHeap.push(i * 1.0 / 3);

	for (int i = 0; i < 10; ++i) {
		std::cout << fibHeap.top() << std::endl;
		fibHeap.pop();
	}*/
}
