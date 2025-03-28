#include "sssq.h"
using namespace std;


int sssq::Process(char* inFile,char* offFile)
{
	ifstream fin1(inFile);
	int d,num;	//#data point
	double x, y;
	fin1 >> d >> num;
	while (fin1 >> x >> y) {
		dataset.push_back(make_pair(x, y));
	}
	char cmd[100] = "D:/bin/qvoronoi o < ";
	strcat(cmd, inFile);
	strcat(cmd, " TO ");
	//char offFile[20] = File;
	strcat(cmd, offFile);
	system(cmd);
	//cout <<"cmd: "<<cmd<<endl;
	ifstream fin(offFile);
	/*char outFile[20];
	strcpy(outFile,inFile);
	strcat(outFile,".out");
	ofstream fout(outFile);*/
	int n, m;	//n: #vertex; m: #edge
	int l, v2;
	fin >> n >> n;
	fin >> m;
	fin >> l;
	for (int i = 0; i < n; i++) {
		fin >> x >> y;
		vertex.push_back(make_pair(x, y));
	}
	vector<int> zero;
	//ofstream fout0("zerotest2");
	for (int i = 0; i < m; i++) {
		fin >> l;
		if (l) {
			bentry be;
			be.minY = be.minX = 10000;
			be.maxY = be.maxX = -10000;
			int j;
			for (j = 0; j < l; j++) {
				fin >> v2;
				if (v2 == 0) {
					continue;
				}
				if (be.maxX < vertex[v2].first) be.maxX = vertex[v2].first;
				if (be.minX > vertex[v2].first) be.minX = vertex[v2].first;
				if (be.maxY < vertex[v2].second) be.maxY = vertex[v2].second;
				if (be.minY > vertex[v2].second) be.minY = vertex[v2].second;
			}
			if (be.maxX > BORDER) be.maxX = BORDER;
			if (be.maxY > BORDER) be.maxY = BORDER;
			if (be.minX < 0) be.minX = 0;
			if (be.minY < 0) be.minY = 0;
			//cout << be.maxX << " " << be.maxY << " " << be.minX << " " << be.minY << endl;
			border.push_back(be);
			for (j = j + 1; j < l; j++)
				fin >> v2;
		}
		else {
			zero.push_back(i);	//l==0: approximate point exits
		}
	}
	vector<pair<double, double> >::iterator begin = dataset.begin();
	for (int i = zero.size() - 1; i >= 0; --i)
		dataset.erase(begin + zero[i]);	//remove the approximate point from the dataset
	//fout <<dataset.size()<<endl;	//#effective data point
	//fout <<n<<endl;			//#vertex
	fin.close();
	fin1.close();
	return 0;
}
/*int main(int argc, char* argv[])
{
	cout<<1<<endl;
	sssq dp;
	dp.Process(argv[1]);
	return 0;
}*/
