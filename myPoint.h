#pragma once
#include "math.h"
class myPoint
{
private:
	double x;
	double y;
public:
	myPoint(){}
	myPoint(double x, double y)
	{
		this->x = x;
		this->y = y;
	}
	void setx(double x){this->x=x;}
	void sety(double y){this->y=y;}
	double getx(){return x;}
	double gety(){return y;}

	myPoint operator +(myPoint p)
	{
		return myPoint(x+p.getx(),y+p.gety());
	}
	myPoint operator /(double l)
	{
		return myPoint(x/l,y/l);
	}
	double distance(myPoint p)
	{
		return sqrt((x-p.getx())*(x-p.getx())+(y-p.gety())*(y-p.gety()));
	}
	bool operator <(const myPoint &p) const
	{
		return x < p.x || (x == p.x && y < p.y);
	}

};