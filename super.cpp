#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <string.h>
#include <map>
#include <assert.h>
#include "super.h"

using namespace std;

void read_dimacs(string name, bool** &conn, int &size) {
	ifstream f (name.c_str());
	string buffer;
	assert(f.is_open());
	set<int> v;
	multimap<int,int> e;
	while (!getline(f, buffer).eof()) {
		if (buffer[0] == 'e') {
			int vi, vj;
			sscanf(buffer.c_str(), "%*c %d %d", &vi, &vj);
			v.insert(vi);
			v.insert(vj);
			e.insert(make_pair(vi, vj));
		}
	}
	size = *v.rbegin();
	conn = new bool*[size];
	for (int i=0; i < size; i++) {
		conn[i] = new bool[size];
		memset(conn[i], 0, size * sizeof(bool));
	}
	for (multimap<int,int>::iterator it = e.begin(); it != e.end(); it++) {
		conn[it->first - 1][it->second - 1] = true;
		conn[it->second - 1][it->first - 1] = true;
	}
	cout << "|E| = " << e.size() << "  |V| = " << v.size() << " p = " << (double) e.size() / (v.size() * (v.size() - 1) / 2) << endl;
	f.close();
}


int main(int argc, char *argv[]) {
	assert(argc == 2);
	cout << "args = " << argv[1] << endl;
	bool **conn;
	int size;
	read_dimacs(argv[1], conn, size);
		clock_t start1 = time(NULL);
		clock_t start2 = clock();
	int *qmax;
	int qsize;
		start1 = time(NULL);
		start2 = clock();
	Maxclique md(conn, size, 0.025);  //(3rd parameter is optional - default is 0.025 - this heuristics parameter enables you to use dynamic resorting of vertices (time expensive)
	// on the part of the search tree that is close to the root - in this case, approximately 2.5% of the search tree -
	// you can probably find a more optimal value for your graphs
	md.mcqdyn(qmax, qsize);  // run max clique with improved coloring and dynamic sorting of vertices 
	sort(qmax, qmax + qsize);
	for(int i = 0; i < qsize - 1; ++i) {
		for(int j = i + 1; j < qsize; ++j) {
			if(!conn[qmax[i]][qmax[j]]) {
				cout<<'('<<qmax[i]<<','<<qmax[j]<<")Failure!!!"<<endl;
			}
		}
	}
	cout << "Maximum clique: ";
	for (int i = 0; i < qsize; i++) 
		cout << qmax[i] + 1<< " ";
	cout << endl;
	cout << "Size = " << qsize << endl;
	cout << "Number of steps = " << md.steps() << endl;
		cout << "Time = " << difftime(time(NULL), start1) << endl;
		cout << "Time (precise) = " << ((double) (clock() - start2)) / CLOCKS_PER_SEC << endl << endl;
	delete [] qmax;
	for (int i=0;i<size;i++)
		delete [] conn[i];
	delete [] conn;
	return 0;
}
