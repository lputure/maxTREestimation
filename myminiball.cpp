#include "miniball.hpp"
#include "myPoint.h"
#include <vector>
#include <list>
#include <iterator>
using namespace std;
void myminiball(myPoint &center, double &radius, vector<myPoint> convexhull)
{
	typedef double mytype;	
	list<vector<double>> lp;

	// transform convexhull to lp
	for(unsigned int i=0;i<convexhull.size();i++)
	{
		vector<double> tp;
		tp.push_back(convexhull[i].getx());
		tp.push_back(convexhull[i].gety());
		lp.push_back(tp);
	}

	// define the types of iterators through the points and their coordinates
	// ----------------------------------------------------------------------
	typedef list<vector<double> >::const_iterator PointIterator; 
	typedef vector<double>::const_iterator CoordIterator;
	// create an instance of Miniball
	// ------------------------------
	typedef Miniball::Miniball <Miniball::CoordAccessor<PointIterator, CoordIterator>> 	MB;
	MB mb(2, lp.begin(), lp.end());

	const double *c = mb.center();
	center.setx(*c);
	c++;
	center.sety(*c);
	//radius = mb.squared_radius();
	radius = sqrt(mb.squared_radius());
}
