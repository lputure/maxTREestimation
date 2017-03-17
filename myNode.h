#pragma once
#include <vector>
using namespace std;
class myNode{
private:
	double position[3];
	double r;
	int N;
	vector<double> gamma;
	vector<double> delta;
	vector<vector<double>> direction;
	vector<myNode*> next;
public:
	myNode(double center[],double radius,int id1, int id2);
};