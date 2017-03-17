// Implementation of Andrew's monotone chain 2D convex hull algorithm.
// Asymptotic complexity: O(n log n).
// Practical performance: 0.5-1.0 seconds for n=1000000 on a 1GHz machine.
#include "myPoint.h"
#include "myVec.h"
#include <algorithm>
#include <vector>
using namespace std;


// 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
// Returns a positive value, if OAB makes a counter-clockwise turn,
// negative for clockwise turn, and zero if the myPoints are collinear.
double cross(myPoint O, myPoint A, myPoint B)
{
	//return (long)(A.getx() - O.getx()) * (B.gety() - O.gety()) - (long)(A.gety() - O.gety()) * (B.getx() - O.getx());
	return (A.getx() - O.getx()) * (B.gety() - O.gety()) - (A.gety() - O.gety()) * (B.getx() - O.getx());
}

// Self Defined Function
// Noter: to delete the points which are the midpoints of an edge like (0,0) between (-1,0) and (1,0)

vector<myPoint> compact(vector<myPoint> convexhull)
{
	vector<myPoint> compactCH;
	int N = convexhull.size();
	for(int i=0;i<N;i++)
	{
		int id1 = i - 1;
		int id2 = i + 1;
		id1 = id1<0? id1+N:id1;
		id2 = id2<0? id2+N:id2;
		id1 = id1>=N? id1-N:id1;
		id2 = id2>=N? id2-N:id2;
		myPoint p = convexhull[i];
		myPoint p1 = convexhull[id1];
		myPoint p2 = convexhull[id2];
		myVec v1(p,p1);
		myVec v2(p,p2);
		if(abs(v1.cross(v2))>1e-4)
			compactCH.push_back(p);
	}
	return compactCH;
}

// Returns a list of myPoints on the convex hull in counter-clockwise order.
// Note: the last myPoint in the returned list is the same as the first one.
vector<myPoint> myconvexhull(vector<myPoint> P)
{
	int n = P.size(), k = 0;
	vector<myPoint> H(2*n);

	// Sort myPoints lexicographically
	sort(P.begin(), P.end());

	// Build lower hull
	for (int i = 0; i < n; ++i) {
		//while (k >= 2 && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
		while (k >= 2 && cross(H[k-2], H[k-1], P[i]) < 0) k--;
		H[k++] = P[i];
	}

	// Build upper hull
	for (int i = n-2, t = k+1; i >= 0; i--) {
		//while (k >= t && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
		while (k >= t && cross(H[k-2], H[k-1], P[i]) < 0) k--;
		if(i==0)
			break;
		H[k++] = P[i];
	}

	H.resize(k);
	
	return compact(H);
}