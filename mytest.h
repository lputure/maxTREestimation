#pragma once
#include "myMultiTree.h"
#include "myconvexhull.h"
#include "mysolution.h"
#include <fstream>
using namespace std;

#define MARGIN 0.3

// an external function
// given the convex hull vertex counterclockwise, find the center of the Minimum Enclosing Circle
void myminiball(myPoint &center, double &radius, vector<myPoint> convexhull); 

class myTest
{
private:
	vector<myPoint> pointlist;  //pointlist of the orginal point
	vector<myPoint> convexhull; //pointlist of the convex hull (the first point is different with the last one)
	myTree tree;				//the tree structure of the feasible solutionset
	vector<myPoint> testpoint;	//pointlist used for test
	vector<mySolution> solutionset;//solution set corresponding to the test point set
public:
	myTest(double X[],double P[],int Xn,int Pn)
	{
		for(int i=0;i<Xn;i++)
		{
			pointlist.push_back(myPoint(X[2*i],X[2*i+1]));
		}
		for(int i=0;i<Pn;i++)
		{
			testpoint.push_back(myPoint(P[2*i],P[2*i+1]));
		}
	}
	void buildTreeStructure()
	{		
		convexhull = myconvexhull(pointlist); // construct the convex hull	
		myPoint center;
		double radius;
		myminiball(center,radius,convexhull); // the center and radius of miniball 
		tree = myTree(convexhull,center,radius); // constructing the tree struction for the solution
	}

	void findsolution(double Y[],double V[],double F[])//
	{
		for(unsigned int i=0;i<testpoint.size();i++)
		{
			mySolution mysl = tree.findsolution(testpoint[i]);
			solutionset.push_back(mysl); 
			Y[2*i] = mysl.pos.getx();
			Y[2*i+1] = mysl.pos.gety();
			V[i] = mysl.val;
			F[i] = mysl.flag==true? 1:0;

		}
	}
	mySolution findsolution(myPoint p)// find the optimal solution for a single test point
	{
		return tree.findsolution(p);
	}
};