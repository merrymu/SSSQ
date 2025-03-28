#include <cstring>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <limits>
#include <math.h>
#include <set>
#include <map>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <chrono>
#include "sssq.h"
using namespace std;

sssq dp;
struct Point
{
	double x;
	double y;
};
struct Rect {
	vector<double> rmax;
	vector<double> rmin;
};
struct Cell
{
	vector<Point> pX;	//vertex sorted by X
	vector<Point> pY;
	vector<int> pi;	//point index
	Point ll;		//lower left corner
	Point ur;
};
struct TreeNode
{
	Cell* c;
	int dim;
	double cut;
        ore_ciphertext ctxt;
	TreeNode* left;
	TreeNode* right;
};
Rect r;
set<int> ans;
vector<vector<double>> Q;
vector<vector<double>> dis;
vector<int> sp;
ore_token rmaxf[2];
ore_token rminf[2];
ore_token rmaxs[2];
ore_token rmins[2];
uint32_t nbits = 32;
element_t k1, k21, k22;
ore_params params;
pairing_t pairing;
element_t g1, g2;
int cnt;
bool toss()
{
	int yingbi = rand() % 100; 
	if (yingbi >= 0 && yingbi <= 49)
	{
		return true;
	}
	else if (yingbi >= 50 && yingbi <= 99)
	{
		return false;
	}
}
bool sortX(Point a, Point b)
{
	return (a.x < b.x);
}
bool sortY(Point a, Point b)
{
	return (a.y < b.y);
}
struct entryComp
{
	bool operator()(Cell a, Cell b)
	{
		if (a.pX.size() < b.pX.size()) return true;
		else return false;
	}
};
void parDCell(vector<int>& cell, double sepLine, vector<int>& left, vector<int>& right, char t)
{
	//cout << sepLine << endl;
	if (t == 'x') {		//using Y
		for (int i = 0; i < cell.size(); ++i) {
			if (dp.dataset[cell[i]].first < sepLine) {
				left.push_back(cell[i]);
				if (dp.border[cell[i]].maxX > sepLine) {
					right.push_back(cell[i]);
				}
			}
			else if (dp.dataset[cell[i]].first > sepLine) {
				right.push_back(cell[i]);
				if (dp.border[cell[i]].minX < sepLine) {
					left.push_back(cell[i]);
				}
			}
			else {
				left.push_back(cell[i]);
				right.push_back(cell[i]);
			}
		}
	}
	else {		//using X
		for (int i = 0; i < cell.size(); ++i) {
			if (dp.dataset[cell[i]].second < sepLine) {
				left.push_back(cell[i]);
				if (dp.border[cell[i]].maxY > sepLine) {
					right.push_back(cell[i]);
				}
			}
			else if (dp.dataset[cell[i]].second > sepLine) {
				right.push_back(cell[i]);
				if (dp.border[cell[i]].minY < sepLine) {
					left.push_back(cell[i]);
				}
			}
			else {
				left.push_back(cell[i]);
				right.push_back(cell[i]);
			}
		}
	}
}
void build_tree(TreeNode* root, double num){
	if (num == 1) {
		root->c->pX.clear();
		root->c->pY.clear();
		/*cout << root->c->ll.x << " " << root->c->ll.y << endl;
		cout << root->c->ur.x << " " << root->c->ur.y << endl;
		for (int i = 0; i < root->c->pi.size(); i++) {
			//cout << root->c->pi[i] << endl;
			cout << dp.dataset[root->c->pi[i]].first << " " << dp.dataset[root->c->pi[i]].second << endl;
		}
		cout << endl;*/
		return;
	}
	int sign = toss();
	int glbInd = -1;
	vector<Point> left, right;
	//in each cell, split it in the "middle"
	int size;
	int l, r, m;
	double middle;
	//cout << root->c->pX.size() << endl;
	//cout << root->c->pY.size() << endl;
	if (sign == 1 && root->c->pX.size() <= 2)
		sign = 0;
	else if (sign == 0 && root->c->pY.size() <= 2)
		sign = 1;
	//cout << sign << endl;
	if (sign == 1) {
		size = root->c->pX.size();
		l = r = m = (size - 1) / 2;
		middle = root->c->pX[m].x;
		if ((size % 2 == 1) || (middle == root->c->pX[m + 1].x)) { //size是单数
			glbInd = m;
		}
		else { //size是双数
			double middle2 = root->c->pX[m + 1].x;
			for (int i = m + 1; i < size; ++i)
				if (root->c->pX[i].x != middle2) {
					r = i;
					break;
				}
			for (int i = m - 1; i > 0; --i)
				if (root->c->pX[i].x != middle) {
					l = i;
					break;
				}
			int n1 = m - l;
			int n2 = r - (m + 1);
			glbInd = (n1 > n2) ? m : m + 1;
		}
	}
	else {
		size = root->c->pY.size();	//using X
		l = r = m = (size - 1) / 2;
		middle = root->c->pY[m].y;
		if ((size % 2 == 1) || (middle == root->c->pY[m + 1].y)) {
			glbInd = m;
		}
		else {
			double middle2 = root->c->pY[m + 1].y;
			for (int i = m + 1; i < size; ++i)
				if (root->c->pY[i].y != middle2) {
					r = i;
					break;
				}
			for (int i = m - 1; i > 0; --i)
				if (root->c->pY[i].y != middle) {
					l = i;
					break;
				}
			int n1 = m - l;
			int n2 = r - (m + 1);
			glbInd = (n1 > n2) ? m : m + 1;
		}
	}
	if(sign == 1){
	root->dim=0;
	root->cut=root->c->pX[glbInd].x;
	}
	else{
	root->dim=1;
	root->cut=root->c->pY[glbInd].y;
	}
	init_ore_ciphertext(root->ctxt, params, pairing, g1, k21);
	ore_encryption(root->ctxt, (uint64_t)(root->cut*1000000), pairing, k1);
//cout<<(uint64_t)(root->cut*1000000)<<endl;
	Cell* a = new Cell;
	Cell* b = new Cell;
	if (sign == 1)	//using Y
	{
        a->ll = root->c->ll;
		a->ur.x = root->c->pX[glbInd].x;
		a->ur.y = root->c->ur.y;
		b->ll.x = root->c->pX[glbInd].x;
		b->ll.y = root->c->ll.y;
		b->ur = root->c->ur;
		left.insert(left.begin(), root->c->pX.begin(), root->c->pX.begin() + glbInd + 1);
		right.insert(right.begin(), root->c->pX.begin() + glbInd, root->c->pX.end());
		for (int i = glbInd + 1; i < root->c->pX.size(); ++i) {
			if (root->c->pX[i].x != root->c->pX[glbInd].x) break;
			left.push_back(root->c->pX[i]);
		}
		for (int i = glbInd - 1; i >= 0; --i) {
			if (root->c->pX[i].x != root->c->pX[glbInd].x) break;
			right.insert(right.begin(), root->c->pX[i]);
		}
		a->pX = left;
		sort(left.begin(), left.end(), sortY);
		a->pY = left;
		b->pX = right;
		sort(right.begin(), right.end(), sortY);
		b->pY = right;
		parDCell(root->c->pi, root->c->pX[glbInd].x, a->pi, b->pi, 'x');
	}
	else
	{
		a->ll = root->c->ll;
		a->ur.x = root->c->ur.x;
		a->ur.y = root->c->pY[glbInd].y;
		b->ll.x = root->c->ll.x;  
		b->ll.y = root->c->pY[glbInd].y;
		b->ur = root->c->ur;
		left.insert(left.begin(), root->c->pY.begin(), root->c->pY.begin() + glbInd + 1);
		right.insert(right.begin(), root->c->pY.begin() + glbInd, root->c->pY.end());
		for (int i = glbInd + 1; i < root->c->pY.size(); ++i) {
			if (root->c->pY[i].y != root->c->pY[glbInd].y) break;
			left.push_back(root->c->pY[i]);
		}
		for (int i = glbInd - 1; i >= 0; --i) {
			if (root->c->pY[i].y != root->c->pY[glbInd].y) break;
			right.insert(right.begin(), root->c->pY[i]);
		}
		a->pY = left;
		sort(left.begin(), left.end(), sortX);
		a->pX = left;
		b->pY = right;
		sort(right.begin(), right.end(), sortX);
		b->pX = right;
		parDCell(root->c->pi, root->c->pY[glbInd].y, a->pi, b->pi, 'y');
	}
	/*for (int i = 0; i < left.size(); i++) {
		cout << "left " << left[i].x << " " << left[i].y << endl;
	}
	for (int i = 0; i < right.size(); i++) {
		cout << "right " << right[i].x << " " << right[i].y << endl;
	}*/
	TreeNode* lt = new TreeNode;
	TreeNode* rt = new TreeNode;
	lt->left = NULL;
	lt->right = NULL;
	rt->left = NULL;
	rt->right = NULL;
	root->left = lt;
	root->right = rt;
	lt->c = a;
	rt->c = b;
	/*int tag = toss();
	if (tag == 1) {
		lt->c = a;
		rt->c = b;
		root->lll = a->ll;
		root->lur = a->ur;
		root->rll = b->ll;
		root->rur = b->ur;
	}
	else {
		lt->c = b;
		rt->c = a;
		root->lll = b->ll;
		root->lur = b->ur;
		root->rll = a->ll;
		root->rur = a->ur;
	}*/
	delete root->c;
	//cout << "left" << endl;
	build_tree(root->left, num / 2);
	//cout << "right" << endl;
	build_tree(root->right, num / 2);
}
void query1(TreeNode* root){
	if (root->left == NULL || root->right == NULL) {
		for (int i = 0; i < root->c->pi.size(); i++)
			ans.insert(root->c->pi[i]);
		return;
	}
	int res1, res2;
	//cout<<res1<<" "<<res2<<endl;
	ore_compare(&res1, root->ctxt, rminf[root->dim], pairing);
	ore_compare(&res2, root->ctxt, rmaxf[root->dim], pairing);
//cout<<root->ctxt<<" "<<rminf[root->dim]<<" "<<rminf[root->dim]<<endl;
	//ore_compare(&res3, root->ctxt, root->ctxt);
	//cout<<res1<<" "<<root->cut<<" "<<r.rmin[root->dim]<<endl;
	//cout<<res2<<" "<<root->cut<<" "<<r.rmax[root->dim]<<endl;
	if (res1==-1) {/*root->ctxt < r.rminc[root->dim]*/
		query1(root->right);
	}
	else if (res2==1) {/*root->ctxt > r.rmaxc[root->dim]*/
		query1(root->left);
	}
	else {
		query1(root->left);
		query1(root->right);
	}
}
void query2(TreeNode* root){
	if (root->left == NULL || root->right == NULL) {
		for (int i = 0; i < root->c->pi.size(); i++)
			ans.insert(root->c->pi[i]);
		return;
	}
	int res1, res2;
	//cout<<res1<<" "<<res2<<endl;
	ore_compare(&res1, root->ctxt, rmins[root->dim], pairing);
	ore_compare(&res2, root->ctxt, rmaxs[root->dim], pairing);
//cout<<root->ctxt<<" "<<rmins[root->dim]<<endl;
	//ore_compare(&res3, root->ctxt, root->ctxt);
	//cout<<res1<<" "<<root->cut<<" "<<r.rmin[root->dim]<<endl;
	//cout<<res2<<" "<<root->cut<<" "<<r.rmax[root->dim]<<endl;
	if (res1==-1) {/*root->ctxt < r.rminc[root->dim]*/
		query2(root->right);
	}
	else if (res2==1) {/*root->ctxt > r.rmaxc[root->dim]*/
		query2(root->left);
	}
	else {
		query2(root->left);
		query2(root->right);
	}
}
void get_distance() {
	double e_d;
	auto it = ans.begin();
	for (; it != ans.end(); it++) {
		vector<double> tmp;
		for (int j = 0; j < Q.size(); j++) {
			e_d = 0;
			e_d += pow(dp.dataset[(*it)].first - Q[j][0], 2);
			e_d += pow(dp.dataset[(*it)].second - Q[j][1], 2);
			tmp.push_back(sqrt(e_d));
		}
		dis.push_back(tmp);
	}
	return;
}
void skyline() {
	int tag;
	for (int i = 0; i < sp.size(); i++) {
		if (sp[i] == -1)
			continue;
		for (int j = 0; j < sp.size(); j++) {
			tag = 0;
			//cout << i << " " << j << endl;
			if (i == j || sp[j] == -1)
				continue;
cnt++;
			for (int k = 0; k < Q.size(); k++) {
				if (dis[sp[i]][k] <= dis[sp[j]][k]) {
					tag++;
				}
			}
			if (tag == Q.size()) {
				//cout << sp[i] << " " << tmp[sp[i]][0] << " " << tmp[sp[i]][1] << endl;
				//cout << sp[j] << " " << tmp[sp[j]][0] << " " << tmp[sp[j]][1] << endl;
				sp[j] = -1;
			}
		}
	}
	set<int> s(sp.begin(), sp.end());
	s.erase(-1);
	sp.assign(s.begin(), s.end());
}
void get_rect(Rect& r, int i) {
	double max, min;
	for (int k = 0; k < 2; k++) {
		max = 0;
		min = BORDER;
		for (int j = 0; j < dis[i].size(); j++) {
			if (Q[j][k] - dis[i][j] < min) {
				min = Q[j][k] - dis[i][j];
			}
			if (Q[j][k] + dis[i][j] > max) {
				max = Q[j][k] + dis[i][j];
			}
		}
		r.rmax.push_back(max);
		r.rmin.push_back(min);
	}
}
bool in_rect(int i, vector<double> p) {
	for (int k = 0; k < 2; k++) {
		if (p[k]<r.rmin[k] || p[k]>r.rmax[k])
			return false;
	}
	return true;
}
void update_rect(int i) {
	Rect rtmp;
	get_rect(rtmp, i);
	for (int k = 0; k < 2; k++) {
		for (int j = 0; j < dis[i].size(); j++) {
			if (rtmp.rmin[k] > r.rmin[k]) {
				r.rmin[k] = rtmp.rmin[k];
			}
			if (rtmp.rmax[k] < r.rmax[k]) {
				r.rmax[k] = rtmp.rmax[k];
			}
		}
		//cout << r.rmin[k] << " " << r.rmax[k] << endl;
	}
}
int main(int argc, char** argv)
{
	double upperbd = BORDER;
	double lowerbd = 0;
cnt=0;
	int mv = 256;//using 4 value to present one cell

	init_pairing(pairing, g1, g2);
	init_ore_params(params, nbits);
	element_init_Zr(k1, pairing);
	element_init_Zr(k21, pairing);
	element_init_Zr(k22, pairing);
	element_random(k1);
	element_random(k21);
	element_random(k22);

	ifstream fin(argv[1]);
	dp.Process(argv[1],argv[2]);
	Cell* c = new Cell;
	Point cor;
	srand(time(NULL));
	for (int i = 1; i < dp.vertex.size(); ++i)
	{
		cor.x = dp.vertex[i].first;
		cor.y = dp.vertex[i].second;
		c->pX.push_back(cor);
		c->pY.push_back(cor);
	}
	cor.x = 0; cor.y = 0;
	c->pX.push_back(cor);
	c->pY.push_back(cor);
	cor.x = BORDER; cor.y = BORDER;
	c->pX.push_back(cor);
	c->pY.push_back(cor);
	sort(c->pX.begin(), c->pX.end(), sortX);
	sort(c->pY.begin(), c->pY.end(), sortY);
	for (int i = 0; i < dp.dataset.size(); ++i)
	{
		c->pi.push_back(i);
	}
	c->ll.x = c->ll.y = lowerbd;
	c->ur.x = c->ur.y = upperbd;
	TreeNode* root = new TreeNode;
	root->c = c;
time_t start, end;
	start = clock();
	build_tree(root, mv);
end = clock();
	cout << "Build tree time: " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
	Q = { {0.47, 0.47},
{0.55, 0.46},
{0.46, 0.54},
{0.53, 0.54},
{0.54, 0.54},
{0.47, 0.53},
{0.53, 0.53},
{0.55, 0.47},
{0.46, 0.53},
{0.54, 0.45} };
	r.rmin.push_back(0.46); r.rmin.push_back(0.45);
	r.rmax.push_back(0.55); r.rmax.push_back(0.54);
//cout<<(uint64_t)(r.rmin[0]*100000000)<<endl;
	init_ore_token(rminf[0], params, pairing, g2, k22);
	init_ore_token(rminf[1], params, pairing, g2, k22);
	init_ore_token(rmaxf[0], params, pairing, g2, k22);
	init_ore_token(rmaxf[1], params, pairing, g2, k22);
	ore_token_gen(rminf[0], (uint64_t)(r.rmin[0]*1000000), pairing, k1);
        ore_token_gen(rminf[1], (uint64_t)(r.rmin[1]*1000000), pairing, k1);
        ore_token_gen(rmaxf[0], (uint64_t)(r.rmax[0]*1000000), pairing, k1);
        ore_token_gen(rmaxf[1], (uint64_t)(r.rmax[1]*1000000), pairing, k1);
	start = clock();
	query1(root);
cout<<"Number1 of retrieved leaves: "<<ans.size()<<endl;
	/*for (auto i = ans.begin(); i != ans.end(); ++i)
	{
		cout << dp.dataset[*i].first << " " << dp.dataset[*i].second << endl;
	}
	cout << endl;*/
	get_distance();
	r.rmin.clear();
	r.rmax.clear();
	get_rect(r, 0);
	auto it = ans.begin();
	it++;
	vector<double> tmp;
	for (int i = 1; i < dis.size(); i++) {
		tmp.push_back(dp.dataset[*it].first);
		tmp.push_back(dp.dataset[*it].second);
		if (in_rect(i, tmp)) {
			update_rect(i);
		}
		tmp.clear();
		it++;
	}
	dis.clear();
//cout<<(uint64_t)(r.rmin[0]*100000000)<<endl;
	init_ore_token(rmins[0], params, pairing, g2, k22);
	init_ore_token(rmins[1], params, pairing, g2, k22);
	init_ore_token(rmaxs[0], params, pairing, g2, k22);
	init_ore_token(rmaxs[1], params, pairing, g2, k22);
	ore_token_gen(rmins[0], (uint64_t)(r.rmin[0]*1000000), pairing, k1);
        ore_token_gen(rmins[1], (uint64_t)(r.rmin[1]*1000000), pairing, k1);
        ore_token_gen(rmaxs[0], (uint64_t)(r.rmax[0]*1000000), pairing, k1);
        ore_token_gen(rmaxs[1], (uint64_t)(r.rmax[1]*1000000), pairing, k1);
	query2(root);
cout<<"Number2 of retrieved leaves: "<<ans.size()<<endl;
	get_distance();
	for (int i = 0; i < dis.size(); i++) {
		sp.push_back(i);
	}
	skyline();
	end = clock();
	cout << "DBSCAN time: " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
	//cout << sp.size() << endl;
cout << "Number of dominance tests: " << cnt << endl;
	vector<int> an;
	an.assign(ans.begin(),ans.end());
	/*for (auto i = sp.begin(); i != sp.end(); ++i)
	{
		//cout << *i << endl;
		//cout << an[*i] << endl;
        cout << dp.dataset[an[*i]].first << " " << dp.dataset[an[*i]].second << endl;
	}*/
        clear_ore_token(rminf[0]);
        clear_ore_token(rminf[1]);
	clear_ore_token(rmaxf[0]);
        clear_ore_token(rmaxf[1]);
	clear_ore_token(rmins[0]);
        clear_ore_token(rmins[1]);
	clear_ore_token(rmaxs[0]);
        clear_ore_token(rmaxs[1]);
        element_clear(k1);
	element_clear(k21);
	element_clear(k22);
	return 0;
}
