
#define LISTMAP(from, to, x, xs, fn) listmap< from, to >(xs, [&] (const from x) {return fn;})
#define LISTMAPREF(from, to, x, xs, fn) listmapref< from, to >(xs, [&] (const from &x) {return fn;})
#define LISTMAPPY(from, x, xs, fn) listmapPy< from >(xs, [&] (const from x) {return fn;})
#define LISTMAPPYREF(from, x, xs, fn) listmapPyRef< from >(xs, [&] (const from &x) {return fn;})
#define LISTMAPFROMPY(to, x, xs, fn) listmapFromPy< to >(xs, [&] (PyObject *x) {return fn;})
#ifndef LIST_MACRO_H
#define LIST_MACRO_H

#include "Python.h"
#include <vector>
#include <functional>
#include <iostream>
using namespace std;
template<class A, class B>
vector<B> listmap(const vector<A> &src, function<B (const A)> fn) {
	vector<B> bs;
	bs.reserve(src.size());
	for (auto it = src.begin(); it!=src.end(); it++) {
		bs.push_back(fn(*it));
	}
	return bs;
}

template<class A, class B>
vector<B> listmapref(const vector<A> &src, function<B (const A&)> fn) {
	vector<B> bs;
	bs.reserve(src.size());
	for (auto it=src.begin(); it!=src.end(); it++) {
		bs.push_back(fn(*it));
	}
	return bs;
}

template<class A>
PyObject *listmapPy(const vector<A> &src, function<PyObject *(const A)> fn) {
	PyObject *l = PyList_New(src.size());
	for (unsigned int i=0; i<src.size(); i++) {
		PyList_SetItem(l, i, fn(src[i]));
	}
	return l;
}

template<class A>
PyObject *listmapPyRef(const vector<A> &src, function<PyObject *(const A &)> fn) {
	PyObject *l = PyList_New(src.size());
	for (unsigned int i=0; i<src.size(); i++) {
		PyList_SetItem(l, i, fn(src[i]));
	}
	return l;
}

template<class A>
vector<A> listmapFromPy(PyObject *src, function<A (PyObject *)> fn) {
	int size = PyList_Size(src);
	if (size == -1) {
		cout << "Trying to map from bad python list" << endl;
	}
	vector<A> as;
	as.reserve(size);
	for (int i=0; i<size; i++) {
		PyObject *o = PyList_GetItem(src, i);
		as.push_back(fn(o));
	}
	return as;
}

#endif
