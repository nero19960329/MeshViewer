#pragma once

#include <vector>

template <class HeapElem, class HeapElemComp>
class MinHeap {
private:
	HeapElem *a;
	int *pos;
	int size_;
	HeapElemComp cmpFun;

public:
	MinHeap() : a(nullptr), size_(0) {}

	MinHeap(std::vector<HeapElem> input, int len) : size_(len) {
		a = new HeapElem[size_];
		for (int i = 0; i < size_; ++i) a[i] = input[i];

		pos = new int[size_];
		for (int i = 0; i < size_; ++i) pos[i] = i;

		for (int i = parent(size_); i >= 0; --i) minHeapify(i);
	}

	~MinHeap() {
		delete[] a;
		delete[] pos;
	}

	unsigned int size() const { return (unsigned int)size_; }

	HeapElem getMinimum() const { return a[0]; }

	HeapElem extractMin() {
		using std::swap;
		using std::get;

		HeapElem res = a[0];
		swap(pos[get<0>(a[0])], pos[get<0>(a[size_ - 1])]);
		a[0] = a[size_ - 1];
		--size_;
		minHeapify(0);
		return res;
	}

	void decreaseKey(const HeapElem& key) {
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
	void minHeapify(int i) {
		using std::get;
		using std::swap;

		int l, r;
		l = left(i);
		r = l + 1;

		int smallest;
		if (l < size_ && cmpFun(a[l], a[i])) smallest = l;
		else smallest = i;

		if (r < size_ && cmpFun(a[r], a[smallest])) smallest = r;

		if (smallest != i) {
			swap(pos[get<0>(a[i])], pos[get<0>(a[smallest])]);
			swap(a[i], a[smallest]);
			minHeapify(smallest);
		}
	}

	inline int left(int i) const {
		return 1 + (i << 1);
	}

	inline int parent(int i) const {
		return ((i - 1) >> 1);
	}
};