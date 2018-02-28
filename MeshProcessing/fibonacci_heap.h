#pragma once

#include <tuple>
#include <vector>

template <class HeapElem, class HeapElemComp>
class FibHeap {
public:
	class FibNode {
	public:
		HeapElem key;
		bool mark;
		FibNode * p, *left, *right, *child;
		int degree;

	public:
		FibNode(const HeapElem & k) : key(k), mark(false), p(nullptr), left(nullptr), right(nullptr), degree(-1) {}
	};

	int n;
	FibNode * min;
	std::vector<FibNode *> pos;
	HeapElemComp cmpFun;

public:
	FibHeap(std::vector<HeapElem> input, int len) {
		n = 0;
		min = nullptr;
		for (int i = 0; i < len; ++i)
			pos.push_back(this->push(input[i]));
	}

	~FibHeap() {
		deleteFibNodes(min);
	}

	void deleteFibNodes(FibNode * x) {
		if (!x) return;

		FibNode * cur = x;
		while (true) {
			if (cur->left && cur->left != x) {
				FibNode * tmp = cur;
				cur = cur->left;
				if (tmp->child)
					deleteFibNodes(tmp->child);
				delete tmp;
			} else {
				if (cur->child)
					deleteFibNodes(cur->child);
				delete cur;
				break;
			}
		}
	}

	void insert(FibNode * x) {
		x->degree = 0;
		x->p = nullptr;
		x->child = nullptr;
		x->mark = false;
		if (min == nullptr)
			min = x->left = x->right = x;
		else {
			min->left->right = x;
			x->left = min->left;
			min->left = x;
			x->right = min;
			if (cmpFun(x->key, min->key))
				min = x;
		}
		++n;
	}

	FibNode * minimum() {
		return min;
	}

	FibNode * extractMin() {
		FibNode * z, *x, *next;
		FibNode ** childList;

		z = min;
		if (z != nullptr) {
			x = z->child;
			if (x != nullptr) {
				childList = new FibNode *[z->degree];
				next = x;
				for (int i = 0; i < z->degree; ++i) {
					childList[i] = next;
					next = next->right;
				}
				for (int i = 0; i < z->degree; ++i) {
					x = childList[i];
					min->left->right = x;
					x->left = min->left;
					min->left = x;
					x->right = min;
					x->p = nullptr;
				}
				delete[] childList;
			}

			z->left->right = z->right;
			z->right->left = z->left;

			if (z == z->right)
				min = nullptr;
			else {
				min = z->right;
				consolidate();
			}
			--n;
		}

		return z;
	}

	void consolidate() {
		FibNode * w, * next, * x, * y, * temp;
		FibNode ** A, ** rootList;
		int d, rootSize;
		int max_degree = int(floor(log(n * 1.0) / log((1.0 + sqrt(5.0)) / 2)));

		A = new FibNode *[max_degree + 2];
		std::fill_n(A, max_degree + 2, nullptr);
		w = min;
		rootSize = 0;
		next = w;

		do {
			rootSize++;
			next = next->right;
		} while (next != w);

		rootList = new FibNode *[rootSize];

		for (int i = 0; i < rootSize; ++i) {
			rootList[i] = next;
			next = next->right;
		}

		for (int i = 0; i < rootSize; ++i) {
			w = rootList[i];
			x = w;
			d = x->degree;
			while (A[d] != nullptr) {
				y = A[d];
				if (!cmpFun(x->key, y->key)) {
					temp = x;
					x = y;
					y = temp;
				}
				fibHeapLink(y, x);
				A[d] = nullptr;
				++d;
			}
			A[d] = x;
		}
		delete[] rootList;
		min = nullptr;
		for (int i = 0; i < max_degree + 2; ++i) {
			if (A[i] != nullptr) {
				if (min == nullptr)
					min = A[i]->left = A[i]->right = A[i];
				else {
					min->left->right = A[i];
					A[i]->left = min->left;
					min->left = A[i];
					A[i]->right = min;
					if (cmpFun(A[i]->key, min->key))
						min = A[i];
				}
			}
		}

		delete[] A;
	}

	void fibHeapLink(FibNode * y, FibNode * x) {
		y->left->right = y->right;
		y->right->left = y->left;
		if (x->child != nullptr) {
			x->child->left->right = y;
			y->left = x->child->left;
			x->child->left = y;
			y->right = x->child;
		} else {
			x->child = y;
			y->right = y;
			y->left = y;
		}
		y->p = x;
		x->degree++;
		y->mark = false;
	}

	void decreaseKey(FibNode * x, const HeapElem & k) {
		FibNode * y;
		if (!cmpFun(k, x->key))
			return;
		x->key = k;
		y = x->p;
		if (y != nullptr && cmpFun(x->key, y->key)) {
			cut(x, y);
			cascadingCut(y);
		}
		if (cmpFun(x->key, min->key))
			min = x;
	}

	void decreaseKey(const HeapElem & k) {
		using std::get;
		this->decreaseKey(pos[get<0>(k)], k);
	}

	void cut(FibNode * x, FibNode * y) {
		if (x->right == x)
			y->child = nullptr;
		else {
			x->right->left = x->left;
			x->left->right = x->right;
			if (y->child == x)
				y->child = x->right;
		}
		y->degree--;

		min->right->left = x;
		x->right = min->right;
		min->right = x;
		x->left = min;
		x->p = nullptr;
		x->mark = false;
	}

	void cascadingCut(FibNode * y) {
		FibNode * z;
		z = y->p;
		if (z != nullptr) {
			if (y->mark == false)
				y->mark = true;
			else {
				cut(y, z);
				cascadingCut(z);
			}
		}
	}

	bool empty() const { return n == 0; }

	FibNode * topNode() { return minimum(); }

	HeapElem top() { return minimum()->key; }

	void pop() {
		if (empty()) return;
		FibNode * x = extractMin();
		if (x)
			delete x;
	}

	FibNode * push(const HeapElem & k) {
		FibNode * x = new FibNode(k);
		insert(x);
		return x;
	}

	unsigned int size() { return (unsigned int)n; }
};