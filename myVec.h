#pragma once
#include "myPoint.h"
#include "math.h"
class myVec{
private:
	double x;
	double y;
public:
	myVec(double x,double y)
	{
		this->x = x;
		this->y = y;
	}
	myVec(myPoint p1, myPoint p2)
	{
		x = p2.getx() - p1.getx();
		y = p2.gety() - p1.gety();
	}
	double getx(){return x;}
	double gety(){return y;}
	myVec operator +(myVec v)
	{
		return myVec(x+v.getx(),y+v.gety());
	}
	
	myVec operator -(myVec v)
	{
		return myVec(x-v.getx(),y-v.gety());
	}
	myVec operator *(double l)
	{
		return myVec(x*l,y*l);
	}

	double dot(myVec v)
	{
		return x*v.getx()+y*v.gety();
	}
	myVec norm()
	{
		if(x==0.0&&y==0.0)
			return myVec(0,0);
		else
		{
			double l = sqrt(x*x+y*y);
			return myVec(x/l,y/l);
		}
	}
	double length()
	{
		return sqrt(x*x+y*y);
	}
	double cross(myVec v)
	{
		return this->x*v.gety()-this->y*v.getx();
	}

};