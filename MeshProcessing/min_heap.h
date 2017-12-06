#pragma once

#include <vector>

template <class HeapElem, class HeapElemComp>
class MinHeap {
private:
	HeapElem *a;
	int *pos;
	int size;
	HeapElemComp cmpFun;

public:
	MinHeap() : a(nullptr), size(0) {}

	MinHeap(std::vector<HeapElem> input, int len) : size(len) {
		a = new HeapElem[size];
		for (int i = 0; i < size; ++i) a[i] = input[i];

		pos = new int[size];
		for (int i = 0; i < size; ++i) pos[i] = i;

		for (int i = parent(size); i >= 0; --i) min_heapify(i);
	}

	~MinHeap() {
		delete[] a;
		delete[] pos;
	}

	int Size() const { return size; }

	HeapElem GetMinimum() const { return a[0]; }

	HeapElem ExtractMin() {
		using std::swap;
		using std::get;

		HeapElem res = a[0];
		swap(pos[get<0>(a[0])], pos[get<0>(a[size - 1])]);
		a[0] = a[size - 1];
		--size;
		min_heapify(0);
		return res;
	}

	void DecreaseKey(const HeapElem& key) {
		using std::get;
		using std::swap;

		int i = pos[get<0>(key)];
		a[i] = key;
		int p = parent(i);
		while (i && cmpFun(a[i], a[p])) {
			swap(pos[get<0>(a[i])], pos[get<0>(a[p])]);
			swap(a[i], a[p]);
			i = p;
			p = parent(i);
		}
	}

private:
	void min_heapify(int i) {
		using std::get;
		using std::swap;

		int l, r;
		l = left(i);
		r = l + 1;

		int smallest;
		if (l < size && cmpFun(a[l], a[i])) smallest = l;
		else smallest = i;

		if (r < size && cmpFun(a[r], a[smallest])) smallest = r;

		if (smallest != i) {
			swap(pos[get<0>(a[i])], pos[get<0>(a[smallest])]);
			swap(a[i], a[smallest]);
			min_heapify(smallest);
		}
	}

	inline int left(int i) const {
		return 1 + (i << 1);
	}

	inline int parent(int i) const {
		return ((i - 1) >> 1);
	}
};