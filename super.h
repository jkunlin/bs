#ifndef MCQD
#define MCQD

#include <iostream>
#include <algorithm>
#include <assert.h>
#ifdef DBG
using namespace std;
#endif

class Maxclique {
	bool **e;
	int pk, level;
	const double Tlimit;
	class Vertices {
		class Vertex {
			int i, d;
		public:
			void set_i(const int ii)  { i = ii; }
			int get_i() const { return i; }
			void set_degree(int dd) { d = dd; }
			int get_degree() const { return d; }
		};
		Vertex *v;
		int sz;
		static bool desc_degree(const Vertex vi, const Vertex vj) { return (vi.get_degree() > vj.get_degree()); }
	public:
#ifdef DBG
		void dbg_v(const string msg="") const {
			std::cout << msg << " Vertices: [";
			for (int i=0; i < sz; i++) 
				std::cout << "(" << v[i].get_i() << "," << v[i].get_degree() << ") ";
			std::cout << "]" << std::endl;
		}
#endif
		Vertices(int size) : sz(0) { v = new Vertex[size]; }
		~Vertices () {}
		void dispose() { if (v) delete [] v; }
		void sort() { std::sort(v, v+sz, desc_degree); }
		void init_colors();
		void set_degrees(Maxclique&);
		int size() const { return sz; }
		void push(const int ii) { v[sz++].set_i(ii); };
		void pop() { sz--; };
		void pop(const int ii);
		Vertex& at(const int ii) const { return v[ii]; };
		Vertex& end() const { return v[sz - 1]; };
	};
	class ColorClass {
		int *i;
		int sz;
	public:
#ifdef DBG
		void dbg_i(const string msg="") const {
			std::cout << msg << " Class: [";
			for (int ii=0; ii < sz; ii++) 
				std::cout << i[ii] << " ";
			std::cout << "]" << std::endl;
		}
#endif
		ColorClass() : sz(0), i(0) {}
		ColorClass(const int sz) : sz(sz), i(0) { init(sz); }
		~ColorClass() { if (i) delete [] i;
		}
		void init(const int sz) { i = new int[sz]; rewind(); }
		void push(const int ii) { i[sz++] = ii; };
		void pop() { sz--; };
		void pop(const int index);
		void rewind() { sz = 0; };
		int size() const { return sz; }
		int& at(const int ii) const { return i[ii]; }
		int& end() const { return i[sz - 1]; }
		ColorClass& operator=(const ColorClass& dh) {
			for (int j = 0; j < dh.sz; j++) i[j] = dh.i[j];
			sz = dh.sz;
			return *this;
		}
	};
	Vertices V;
	ColorClass *C, QMAX, Q, CLIQUE_VERTEX;
	class StepCount {
		int i1, i2;
	public:
		StepCount() : i1(0), i2(0) {}
		void set_i1(const int ii)  { i1 = ii; }
		int get_i1() const { return i1; }
		void set_i2(const int ii)  { i2 = ii; }
		int get_i2() const { return i2; }
		void inc_i1()  { i1++; }
	};
	StepCount *S;
	bool connection(const int i, const int j) const { return e[i][j]; }
	bool conflict(const int, const ColorClass&);
	void cut(const Vertices&, Vertices&);
	void color_sort(Vertices&);
	void re_color_sort(Vertices&);
	int re_color(const int);
	int count_conflict(int, const ColorClass&, int&);
	bool detect_clique(Vertices&);
	void detect_clique_vertices(Vertices&);
	void expand(Vertices);
	void expand_dyn(Vertices);
	void _mcq(int*&, int&, bool);
	void degree_sort(Vertices &R) { R.set_degrees(*this); R.sort(); }
public:
#ifdef DBG
	void dbg_C() const {
		for (int i=0; i < V.size(); i++) {
			std::cout << "C["<< i << "] : ";
			C[i].dbg_i();
		}
	}
	void dbg_conn() const {
		for (int i=0; i < V.size(); i++) {
			for (int j=0; j < V.size(); j++) {
				std::cout <<e[i][j];
			}
			std::cout<< std::endl;
		}
	}
#endif
	Maxclique(bool **, const int, const double=0.025);
	int steps() const { return pk; }
	void mcq(int* &maxclique, int &sz) { _mcq(maxclique, sz, false); }
	void mcqdyn(int* &maxclique, int &sz) { _mcq(maxclique, sz, true); }
	void permute_vertices(); 
	~Maxclique() {
		if (C) delete [] C;
		if (S) delete [] S;
		V.dispose();
	};
};

Maxclique::Maxclique (bool **conn, const int sz, const double tt) : pk(0), level(1), Tlimit(tt), V(sz), Q(sz), CLIQUE_VERTEX(sz), QMAX(sz) {
	assert(conn!=0 && sz>0);
	for (int i=0; i < sz; i++) V.push(i);
	e = conn;
	C = new ColorClass[sz + 1];
	for (int i=0; i < sz + 1; i++) C[i].init(sz);
	S = new StepCount[sz + 1];
}

void Maxclique::permute_vertices() {
	int sz = V.size();
	bool **e_tmp = new bool *[sz];
	for (int i=0; i < sz; i++) {
		e_tmp[i] = new bool[sz];
	}
	for (int i = 0; i < sz; ++i) {
		for(int j = 0; j < sz; ++j)
			e_tmp[i][j] = e[i][j];
	}
	
	
	for (int i = 0; i < sz; ++i) {
		for (int j = 0; j < sz; ++j) {
			e[i][j] = e_tmp[V.at(i).get_i()][V.at(j).get_i()];
		}
	}
	
	for (int i = 0; i < sz; ++i) {
		V.at(i).set_i(i);
	}
	delete [] e_tmp;

}

void Maxclique::_mcq(int* &maxclique, int &sz, bool dyn) { 
	V.set_degrees(*this);
	V.sort();
	permute_vertices();
	V.init_colors();
	if (dyn) {
		for (int i=0; i < V.size() + 1; i++) {
			S[i].set_i1(0);
			S[i].set_i2(0);
		}
		expand_dyn(V);
	}
	else
		expand(V);
	maxclique = new int[QMAX.size()]; 
	for (int i=0; i<QMAX.size(); i++) { 
		maxclique[i] = QMAX.at(i);
	}
	sz = QMAX.size();
}

void Maxclique::Vertices::init_colors() { 
	const int max_degree = v[0].get_degree();
	for (int i = 0; i < max_degree; i++)
		v[i].set_degree(i + 1);
	for (int i = max_degree; i < sz; i++)
		v[i].set_degree(max_degree + 1);
}

void Maxclique::Vertices::set_degrees(Maxclique &m) { 
	for (int i=0; i < sz; i++) {
		int d = 0;
		for (int j=0; j < sz; j++)
			if (m.connection(v[i].get_i(), v[j].get_i())) d++;
		v[i].set_degree(d);
	}
}
void Maxclique::Vertices::pop(const int ii) {
	int iter = 0;
	while (v[iter++].get_i() != ii);
	for (; iter < sz; ++iter) {
		v[iter - 1] = v[iter];
		v[iter - 1].set_degree(v[iter].get_degree() - 1);
	}
	pop();
}

void Maxclique::ColorClass::pop(const int index){
	for (int iter = index + 1; iter < sz; ++iter) {
		i[iter - 1] = i[iter];
	}
	pop();
}

bool Maxclique::conflict(const int pi, const ColorClass &A) {
	for (int i = 0; i < A.size(); i++)
		if (connection(pi, A.at(i)))
			return true;
	return false;
}

void Maxclique::cut(const Vertices &A, Vertices &B) {
	for (int i = 0; i < A.size() - 1; i++) {
		if (connection(A.end().get_i(), A.at(i).get_i()))
			B.push(A.at(i).get_i());
	}
}

void Maxclique::color_sort(Vertices &R) {
	int j = 0;
	int maxno = 1;
	int min_k = QMAX.size() - Q.size() + 1;
	C[1].rewind();
	C[2].rewind();
	int k = 1;
	for (int i=0; i < R.size(); i++) {
		int pi = R.at(i).get_i();
		k = 1;
		while (conflict(pi, C[k]))
			k++;
		if (k > maxno) {
			maxno = k;
			C[maxno + 1].rewind();
		}
		C[k].push(pi);
		if (k < min_k) {
			R.at(j++).set_i(pi);
		}
	}

	/*
	j = R.size() - 1;
	if(min_k <= 0) min_k = 1;
	//normal branch
	for(k = maxno; k >= min_k; --k) {
		for(int i = C[k].size() - 1; i >= 0; --i) {
			R.at(j).set_i(C[k].at(i));
			R.at(j--).set_degree(k);
		} 
	}
	R.at(j).set_degree(0);
	for(; k >= 1; --k) {
		for(int i = C[k].size() - 1; i >= 0; --i) {
			R.at(j--).set_i(C[k].at(i));
		} 
	}
	
	//R.at(j).set_degree(0);
	for( k = min_k - 1; k > 0; --k) {
		for (int i = C[k].size() - 1; i >= 0; --i)
			R.at(j).set_i(C[k].at(i));
			R.at(j--).set_degree(k);
	}
	*/

	if (j > 0) R.at(j-1).set_degree(0);
	if (min_k <= 0) min_k = 1;
	for (k = min_k; k <= maxno; k++)
		for (int i = 0; i < C[k].size(); i++) {
			R.at(j).set_i(C[k].at(i));
			R.at(j++).set_degree(k);
		}
}

void Maxclique::re_color_sort(Vertices &R) {
	int j = 0;
	int maxno = 1;
	int min_k = QMAX.size() - Q.size() + 1;
	C[1].rewind();
	C[2].rewind();
	int k = 1;
	for (int i=0; i < R.size(); i++) {
		int pi = R.at(i).get_i();
		k = 1;
		while (conflict(pi, C[k]))
			k++;
		C[k].push(pi);
		if (k > maxno) {
			maxno = k;
			C[maxno + 1].rewind();
			
			if (k >= min_k) {
				k = re_color(k);
				if(C[maxno].size() == 0)
					maxno -= 1;
			}
			
		}
		if (k < min_k) {
			R.at(j++).set_i(pi);
		}
	}
	if (j > 0) R.at(j-1).set_degree(0);
	if (min_k <= 0) min_k = 1;
	for (k = min_k; k <= maxno; k++)
		for (int i = 0; i < C[k].size(); i++) {
			R.at(j).set_i(C[k].at(i));
			R.at(j++).set_degree(k);
		}
	/*
	int j = R.size() - 1;
	if(min_k <= 0) min_k = 1;
	//last color maybe better than this one, try it when no bug
	
	//rever branch
	for(k = maxno; k >= 1; --k) {
		for(int i = 0; i < C[k].size(); ++i) {
			R.at(j).set_i(C[k].at(i));
			R.at(j--).set_degree(k);
		} 
	}
	
	//normal branch
	for(k = maxno; k >= 1; --k) {
		for(int i = C[k].size() - 1; i >= 0; --i) {
			R.at(j).set_i(C[k].at(i));
			R.at(j--).set_degree(k);
		} 
	}
	*/	
}

int Maxclique::re_color(const int k) {
	int min_k = QMAX.size() - Q.size() + 1;
	int pi = C[k].end();
	int q_ci = -1;
	int qi = -1;
	for(int k1 = 1; k1 < min_k - 1; ++k1) {
		if(count_conflict(pi, C[k1], q_ci) == 1) {
			qi = C[k1].at(q_ci);
			for(int k2 = k1 + 1; k2 <= min_k - 1; ++k2) {
				if(!conflict(qi, C[k2])) {
					C[k].pop();
					C[k1].push(pi);
					C[k2].push(qi);
					C[k1].pop(q_ci);
					return k1;
				}
			}
		}
	}
	return k;
}

int Maxclique::count_conflict(int pi, const ColorClass& A, int &q_ci) {
	int num = 0;
	for(int i = 0; i < A.size(); ++i) {
		if(connection(pi, A.at(i))) {
				q_ci = i;
				++num;
		}
	}
	return num;
}


bool Maxclique::detect_clique(Vertices &R) {
	if (R.end().get_degree() == R.size()) {
		if (Q.size() + R.size() > QMAX.size()) {
			QMAX = Q;
			for (int i = 0; i < R.size(); ++i) {
				QMAX.push(R.at(i).get_i());
			}
			std::cout << "step = " << pk << " current max. clique size = " << QMAX.size() << std::endl; 
		}
		return true;
	}
	return false;
}

void Maxclique::detect_clique_vertices(Vertices &R) {
	int max_no = R.end().get_degree();
	int clq_v = -1;
	for (int k = 1; k <= max_no; ++k) {
		if (C[k].size() > 1)
			continue;
		clq_v = C[k].at(0);
		for (int l = 1; l < k; ++l) {
			for (int i = 0; i < C[l].size(); ++i) {
				if (!connection(clq_v, C[l].at(i)))
					goto next_color;
			}
		}
		Q.push(clq_v);
		CLIQUE_VERTEX.push(clq_v);
		R.pop(clq_v); //how to do ?
next_color:
		;
	}
	if (Q.size() > QMAX.size())
		QMAX = Q;
}



void Maxclique::expand(Vertices R) {
	while (R.size()) {
		if (Q.size() + R.end().get_degree() > QMAX.size()) {
			Q.push(R.end().get_i());
			Vertices Rp(R.size());
			cut(R, Rp);
			if (Rp.size()) {
				color_sort(Rp);
				pk++;
				expand(Rp);
			}
			else if (Q.size() > QMAX.size()) { 
				std::cout << "step = " << pk << " current max. clique size = " << Q.size() << std::endl; 
				QMAX = Q;
			}    
			Rp.dispose();
			Q.pop();
		}
		else {
			return;
		}
		R.pop();
	}
}

void Maxclique::expand_dyn(Vertices R) {
	S[level].set_i1(S[level].get_i1() + S[level - 1].get_i1() - S[level].get_i2());
	S[level].set_i2(S[level - 1].get_i1());
	while (R.size()) {
		if (Q.size() + R.end().get_degree() > QMAX.size()) {
			Q.push(R.end().get_i());
			Vertices Rp(R.size());
			cut(R, Rp);
			if (Rp.size()) {
				if ((float)S[level].get_i1()/++pk < Tlimit) {
					degree_sort(Rp);
				}
				re_color_sort(Rp);
			/*	
				if (detect_clique(Rp)) {
					goto next_vertex;
				}
				else
					detect_clique_vertices(Rp);
			*/		
				S[level].inc_i1();
				level++;
				expand_dyn(Rp);
				level--;
			}
			else if (Q.size() > QMAX.size()) { 
				std::cout << "step = " << pk << " current max. clique size = " << Q.size() << std::endl; 
				QMAX = Q;
			}    
next_vertex:
			Rp.dispose();
			while(CLIQUE_VERTEX.size() != 0 && Q.end() == CLIQUE_VERTEX.end()) {
				Q.pop();
				CLIQUE_VERTEX.pop();
			}
			Q.pop();
		}
		else {
			return;
		}
		R.pop();
	}
}

#endif
