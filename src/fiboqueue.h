/*
 * We need:
 * 	PUSH
 * 	POP
 * 	keyOf
 * 	contains
 * 	decreaseKey
 * 	isEmpty
 */

/**
 * Fibonacci Heap
 * Copyright (c) 2014, Emmanuel Benazera beniz@droidnik.fr, All rights reserved.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library.
 */
#ifndef FIBOHEAP_H
#define FIBOHEAP_H

#include <cstddef>
#include <math.h>
#include <limits>
#include <iostream>
#include <unordered_map>
#include <algorithm>

template<class T, class V>
class FibHeap {
private:
	// node
	class FibNode {
	public:
		FibNode(T k)
			:key(k),mark(false),p(nullptr),left(nullptr),right(nullptr),child(nullptr),degree(-1) {
		}
		T key;
		bool mark;
		FibNode *p;
		FibNode *left;
		FibNode *right;
		FibNode *child;
		bool present = true;
		int degree;
	}; // end FibNode

	typename std::unordered_map<T, FibNode*> 					 NodeMap;
	typename std::unordered_map<T, FibNode*>::iterator MapIterator;

	int n;
	int counter;
	FibNode *min;

	std::unordered_multimap<T,FibNode*> fstore;
	std::unique_ptr< FibNode[] > index;

	FibHeap(T N) :n(0),min(nullptr) {
		index  = std::unique_ptr< FibNode[] >(new FibNode[N]);
	}

	~FibHeap() {
		// delete all nodes.
		delete_fibnodes(min);
	}

	void delete_fibnodes(FibNode *x) {
		if (!x) return;
		FibNode *cur = x;
		while(true) {
			if (cur->left && cur->left != x) {
				FibNode *tmp = cur;
				cur = cur->left;
				if (tmp->child)
					delete_fibnodes(tmp->child);
				delete tmp;
			} else {
				if (cur->child) delete_fibnodes(cur->child);
				delete cur;
				break;
			}
		}
	}

	/*
	* insert(x)
	* 1. x.degree = 0
	* 2. x.p = NIL
	* 3. x.child = NIL
	* 4. x.mark = FALSE
	* 5. if H.min == NIL
	* 6. 	create a root list for H containing just x
	* 7. 	H.min = x
	* 8. else insert x into H's root list
	* 9. 	if x.key < H.min.key
	*10. 		H.min = x
	*11. H.n = H.n + 1
	*/
	void insert(FibNode *x) {
		x->degree = 0;
		x->p = nullptr;
		x->child = nullptr;
		x->mark = false;
		if ( min == nullptr)  min = x->left = x->right = x;
		else {
			min->left->right = x;
			x->left = min->left;
			min->left = x;
			x->right = min;
			if ( x->key < min->key ) min = x;
		}
		++n;
		index[counter++] = x;
	}

	/*
	* The minimum node of the heap.
	*/
	FibNode* minimum() const {
		return min;
	}

	/*
	* extract_min
	* 1. z = H.min
	* 2. if z != NIL
	* 3. 	for each child x of z
	* 4. 		add x to the root list of H
	* 5. 		x.p = NIL
	* 6. 	remove z from the root list of H
	* 7.		if z == z.right
	* 8. 		H.min = NIL
	* 9. 	else H.min = z.right
	*10. 		CONSOLIDATE(H)
	*11. 	H.n = H.n - 1
	*12. return z
	*/
	FibNode* extract_min() {
		FibNode *z, *x, *next;
		FibNode ** childList;

		// 1
		z = min;
		if ( z != nullptr ) {
			x = z->child;
			if ( x != nullptr ) {
				childList = new FibNode*[z->degree];
				next = x;
				for ( int i = 0; i < (int)z->degree; i++ ) {
					childList[i] = next;
					next = next->right;
				}
				for ( int i = 0; i < (int)z->degree; i++ ) {
					x = childList[i];
					min->left->right = x;
					x->left = min->left;
					min->left = x;
					x->right = min;
					x->p = nullptr;
				}
				delete [] childList;
			}
			z->left->right = z->right;
			z->right->left = z->left;
			if ( z == z->right ) min = nullptr;
			else {
				min = z->right;
				consolidate();
			}
			n--;
		}
		z -> present = false;
		return z;
	}

	/*
	* consolidate
	* 1. let A[0 . . D(H.n)] be a new array
	* 2. for i = 0 to D(H.n)
	* 3. 	A[i] = NIL
	* 4. for each node w in the root list of H
	* 5. 	x = w
	* 6. 	d = x.degree
	* 7. 	while A[d] != NIL
	* 8. 		y = A[d]
	* 9. 		if x.key > y.key
	*10.			exchange x with y
	*11. 		FIB-HEAP-LINK(H,y,x)
	*12. 		A[d] = NIL
	*13. 		d = d + 1
	*14. 	A[d] = x
	*15. H.min = NIL
	*16. for i = 0 to D(H.n)
	*17. 	if A[i] != NIL
	*18. 		if H.min == NIL
	*19. 			create a root list for H containing just A[i]
	*20. 			H.min = A[i]
	*21. 		else insert A[i] into H's root list
	*22. 			if A[i].key < H.min.key
	*23. 				H.min = A[i]
	*/
	void consolidate() {
		FibNode* w, * next, * x, * y, * temp;
		FibNode** A, ** rootList;
		// Max degree <= log base golden ratio of n
		int d, rootSize;
		int max_degree = static_cast<int>(floor(log(static_cast<double>(n))/log(static_cast<double>(1 + sqrt(static_cast<double>(5)))/2)));
		A = new FibNode*[max_degree+2]; // plus two both for indexing to max degree and so A[max_degree+1] == NIL
		std::fill_n(A, max_degree+2, nullptr);
		w = min;
		rootSize = 0;
		next = w;
		do {
			rootSize++;
			next = next->right;
		} while ( next != w );
		rootList = new FibNode*[rootSize];
		for ( int i = 0; i < rootSize; i++ ) {
			rootList[i] = next;
			next = next->right;
		}
		for ( int i = 0; i < rootSize; i++ ) {
			w = rootList[i];
			x = w;
			d = x->degree;
			while ( A[d] != nullptr )	{
				y = A[d];
				if ( x->key > y->key ) {
					temp = x;
					x = y;
					y = temp;
				}
				fib_heap_link(y,x);
				A[d] = nullptr;
				d++;
			}
			A[d] = x;
		}
		delete [] rootList;
		min = nullptr;
		for ( int i = 0; i < max_degree+2; i++ ) {
			if ( A[i] != nullptr ) {
				if ( min == nullptr ) min = A[i]->left = A[i]->right = A[i];
				else {
					min->left->right = A[i];
					A[i]->left = min->left;
					min->left = A[i];
					A[i]->right = min;
					if ( A[i]->key < min->key ) min = A[i];
				}
			}
		}
		delete [] A;
	}

	/*
	* fib_heap_link(y,x)
	* 1. remove y from the root list of heap
	* 2. make y a child of x, incrementing x.degree
	* 3. y.mark = FALSE
	*/
	void fib_heap_link( FibNode* y, FibNode* x ) {
		// 1
		y->left->right = y->right;
		y->right->left = y->left;
		// 2
		if ( x->child != nullptr ) {
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
		// 3
		y->mark = false;
	}



	/*
	* cut(x,y)
	* 1. remove x from the child list of y, decrementing y.degree
	* 2. add x to the root list of H
	* 3. x.p = NIL
	* 4. x.mark = FALSE
	*/
	void cut( FibNode* x, FibNode* y ) {
		if ( x->right == x ) {
			y->child = nullptr;
		}	else {
			x->right->left = x->left;
			x->left->right = x->right;
			if ( y->child == x ) {
				y->child = x->right;
			}
		}
		y->degree--;
		min->right->left = x;
		x->right = min->right;
		min->right = x;
		x->left = min;
		x->p = nullptr;
		x->mark = false;
	}

	/*
	* cascading_cut(y)
	* 1. z = y.p
	* 2. if z != NIL
	* 3. 	if y.mark == FALSE
	* 4. 		y.mark = TRUE
	* 5. 	else CUT(H,y,z)
	* 6. 		CASCADING-CUT(H,z)
	*/
	void cascading_cut( FibNode* y ) {
		FibNode* z;
		z = y->p;
		if ( z != nullptr ) {
			if ( y->mark == false ) {
				y->mark = true;
			} else {
				cut(y,z);
				cascading_cut(z);
			}
		}
	}

	/*
	* set to infinity so that it hits the top of the heap, then easily remove.
	*/
	void remove_fibnode( FibNode* x ) {
		decrease_key(x, - INFINITY);
		FibNode *fn = extract_min();
		delete fn;
	}

public:

	/*
	 * decrease_key(x,k)
	 * 1. if k > x.key
	 * 2. 	error "new key is greater than current key"
	 * 3. x.key = k
	 * 4. y = x.p
	 * 5. if y != NIL and x.key < y.key
	 * 6. 	CUT(H,x,y)
	 * 7. 	CASCADING-CUT(H,y)
	 * 8. if x.key < H.min.key
	 * 9. 	H.min = x
	 */
	void decrease_key( FibNode* x, T k ) {
		typename std::unordered_map<T, FibNode*>::iterator mit = find(x->key);
		fstore.erase(mit);
		fstore.insert(std::pair<T, FibNode*>(k,x));

		FibNode* y;

		if ( k > x->key ) return;
		x->key = k;
		y = x->p;
		if ( y != nullptr && x->key < y->key ) {
			cut(x,y);
			cascading_cut(y);
		}
		if ( x->key < min->key ) min = x;
	}

	/*
	* mapping operations to STL-compatible signatures.
	*/
	bool empty() const {
		return n == 0;
	}

	T minKey() const {
		return minimum()->key;
	}

	FibNode* push(T k) {
		FibNode *x = new FibNode(k);
		insert(x);
		fstore.insert(std::pair<T,FibNode>(k,x));
		return x;
	}

	unsigned int size() const {
		return (unsigned int) n;
	}

	FibNode* findIndex(V x) {
		return index[x];
	}

	typename std::unordered_map<T, FibNode*>::iterator find(T k) {
		typename std::unordered_map<T, FibNode*>::iterator mit = fstore.find(k);
		return mit;
	}

	void pop() {
		if (empty()) return;
		FibNode *x = extract_min();
		if (!x)  return; // should not happen.
		auto range = fstore.equal_range(x->key);
		auto mit = std::find_if(range.first, range.second,
                          [x](const std::pair<T, FibNode*> &ele){
                          	return ele.second == x;
                          }
		);
		if (mit != range.second) fstore.erase(mit);
		else std::cerr << "[Error]: key " << x->key << " cannot be found in FiboQueue fast store\n";
		delete x;
	}
};

#endif
