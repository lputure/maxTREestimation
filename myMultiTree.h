#pragma once
#include "myPoint.h"
#include "myVec.h"
#include "mysolution.h"
#include <vector>
#include <iostream>
using namespace std;

#define EPS 1e-6

class myNode{
public:
	myPoint pos;			//postion of the node
	double radius;			//radius of the minimum circle with pos as the center
	int N;					//number of branches
	vector<double> gamma;	//distance to each supporting edges
	vector<double> delta;	//length of every branch
	vector<myVec> direction;//direction of each branch
	vector<myNode*> next;	//point to the child node
};

class myTree{
private:
	myNode *root;
	vector<myPoint> convexhull;
	vector<bool> flag; //indicate whether the root node lies on the diameter of MEC
public:
	myTree(){};
	myTree(vector<myPoint> convexhull, myPoint miniballcenter, double radius)
	{
		this->convexhull = convexhull;
		root = new myNode;
		root->pos = miniballcenter;
		root->radius = radius;		

		vector<int> supportid; //find the points lyiong on MEC
		for(unsigned int i=0;i<convexhull.size();i++)
		{
			double tdis = convexhull[i].distance(miniballcenter);
			if(tdis<radius+EPS && tdis>radius-EPS) //for numerical stability
			{
				supportid.push_back(i);
			}
		}
		root->N = supportid.size();

		for(unsigned int i=0;i<supportid.size();i++)//traverse all edges with both endpoints lying on MEC
		{
			int id1 = supportid[i];
			int id2 = supportid[(i+1)%supportid.size()];
			myPoint x1 = convexhull[id1];
			myPoint x2 = convexhull[id2];	//take out two points

			if(fabs(cross(x1,x2,miniballcenter))<EPS)//x1,x2 lies on the diameter
			{
				this->flag.push_back(true);
				myVec v(miniballcenter,x1);
				myVec vp(v.gety(),-v.getx());
				root->direction.push_back(vp.norm());
				root->gamma.push_back(0);
				if(id1+1==id2 || (id1+1)%convexhull.size()==id2) //there exist other points between these two points
				{
					root->next.push_back(NULL);
					root->delta.push_back(-1); //delta
				}else //not exist
				{
					myNode *leaf = buildleaf(root->pos,root->direction[i],id1,id2); //construct the leaf node
					root->next.push_back(leaf); //connect the leaf node with the root node
					root->delta.push_back(miniballcenter.distance(leaf->pos)); //delta
				}
			}else //x1,x2 do not lie on the diameter
			{
				this->flag.push_back(false);
				myVec v1(x1,miniballcenter);
				myVec v2(x2,miniballcenter);
				myVec v = v1 + v2;
				myVec dir = v.norm();
				root->direction.push_back(dir); //direction
				root->gamma.push_back(miniballcenter.distance((x1+x2)/2)); //gamma
				if(id1+1==id2 || (id1+1)%convexhull.size()==id2) //there is no other points between these two points
				{
					root->next.push_back(NULL);
					root->delta.push_back(-1); //an infinite branch
				}else //exist
				{
					myNode *leaf = buildleaf(root->pos,root->direction[i],id1,id2); //construct the leaf node
					root->next.push_back(leaf); //connect the leaf node with the root node
					root->delta.push_back(miniballcenter.distance(leaf->pos)); //delta
				}			
			}
		}
	}

	myNode *buildleaf(myPoint parentnodeposition, myVec v, int j1, int j2)
	{
		myNode *leafnode = new myNode;

		if(j2<j1) j2 = j2 + convexhull.size(); //modify the second point

		vector<double> tdelta; //compute the candidate position of the leaf node
		for(int i=j1+1;i<j2;i++)
		{
			myVec y1(convexhull[j1],parentnodeposition); //j1
			myVec y3(convexhull[i%convexhull.size()],parentnodeposition); //point i between j1,j2
			double tmp = (y3.length()*y3.length()-y1.length()*y1.length())/(2*v.dot(y1-y3));
			tdelta.push_back(tmp);
		}
		double vmin = tdelta[0]; //use the closest one as the postion of the leaf node
		if(tdelta.size()>1)
		{
			for(unsigned int i=1;i<tdelta.size();i++)
				if(tdelta[i]<vmin)
					vmin = tdelta[i];
		}

		leafnode->pos = myPoint(parentnodeposition.getx()+v.getx()*vmin,parentnodeposition.gety()+v.gety()*vmin);
		leafnode->radius = leafnode->pos.distance(convexhull[j1]);
		
		tdelta.push_back(vmin);//add the value of j2 to parameter
		int pid = j1; //pid point to j1
		int cnt = 0;
		for(unsigned int i=0;i<tdelta.size();i++)
		{			
 			if(tdelta[i]<vmin+EPS && tdelta[i]>vmin-EPS)//for numerical stability
			{				
				int id1 = pid;
				int id2 = (j1+1+i)%convexhull.size();
				myPoint x1 = convexhull[id1];
				myPoint x2 = convexhull[id2];
				myVec v1(x1,leafnode->pos);
				myVec v2(x2,leafnode->pos);
				myVec v = v1 + v2;
				myVec dir = v.norm();
				leafnode->direction.push_back(dir);
				cnt = cnt + 1;
				leafnode->gamma.push_back(leafnode->pos.distance((x1+x2)/2)); //same with the construction of root node
				if(id1+1==id2 || (id1+1)%convexhull.size()==id2)
				{
					leafnode->next.push_back(NULL);
					leafnode->delta.push_back(-1);
				}else
				{
					myNode *nextleaf = buildleaf(leafnode->pos,leafnode->direction[cnt-1],id1,id2); //recursion
					leafnode->next.push_back(nextleaf); 
					leafnode->delta.push_back(leafnode->pos.distance(nextleaf->pos));
				}
				pid = id2; //update pid
			}
		}
		leafnode->N = leafnode->next.size();
		return leafnode;
	}

	double cross(myPoint p0, myPoint p1, myPoint p2)
	{
		return (p1.getx()-p0.getx())*(p2.gety()-p0.gety())-(p2.getx()-p0.getx())*(p1.gety()-p0.gety());
	}

	bool isInConvexHull(myPoint p) //decide whether p lies inside convex hull, O(N) complexity
	{
		int n = convexhull.size();
		if(n<3) return false;
		for(int i=0;i<n;i++)//convex hull vertex is arranged counterclockwise
		{
			int id1 = i;
			int id2 = (i+1)%n;
			if(cross(convexhull[id1],p,convexhull[id2])>0)
				return false;
		}
		return true;
	}

	bool isInBranch(myPoint p,myNode node, int i) //decide whether p belongs to the i-th branch of node
	// return true, if p lies inside the supporting circle and outside the convex hull
	{
		double radius = node.radius*node.radius/2/node.gamma[i]; // r'=r*r/(2*gamma)
		myPoint center;
		center.setx(node.pos.getx()-radius*node.direction[i].getx());
		center.sety(node.pos.gety()-radius*node.direction[i].gety());
		myVec po(p, node.pos);
		return center.distance(p)<radius && node.gamma[i]<po.dot(node.direction[i]);
	}

	bool isInBranch(myPoint p,myNode node, int i,bool flag) ///decide whether p belongs to the i-th branch of node, used for root node
	{
		if(flag==false)
		{
			double radius = node.radius*node.radius/2/node.gamma[i]; // r'=r*r/(2*gamma)
			myPoint center;
			center.setx(node.pos.getx()-radius*node.direction[i].getx());
			center.sety(node.pos.gety()-radius*node.direction[i].gety());			
			myVec po(p, node.pos);
			return center.distance(p)<radius && node.gamma[i]<po.dot(node.direction[i]);
		}else
		{
			myVec po(p, node.pos);
			return node.gamma[i]<po.dot(node.direction[i]);
		}
	}
	
	mySolution findsolutionInTheNode(myNode node, myPoint p)//find the optimal solution from the present node
	{
		mySolution solution;
		for(int i=0;i<node.N;i++)//traverse all branches
		{
			if(isInBranch(p,node,i))//in the i-th branch
			{
				myVec po(p, node.pos);
				double kcosbeta = po.dot(node.direction[i]);
				double ksquare = po.length()*po.length();
				double a = node.gamma[i] - kcosbeta;
				double b = node.radius*node.radius - ksquare;
				double c = node.radius*node.radius*kcosbeta - ksquare*node.gamma[i]; //quadratic equation
					//(gamma-k*(cos(theta))lambda^2 + (r^2-k^2)lambda + (kcos(theta)*r^2 - gamma*k^2) = 0 
				double d = -b/(2*a) + sqrt((b*b-4*a*c)/(4*a*a)); //solution should be larger than 0				
				if(node.delta[i]==-1 || node.delta[i]>d) // solution does not exceed the present branch
				{
					solution.flag = false;
					solution.pos.setx(node.pos.getx()+node.direction[i].getx()*d);
					solution.pos.sety(node.pos.gety()+node.direction[i].gety()*d);
					solution.val = solution.pos.distance(p)/sqrt(d*d + node.radius*node.radius + 2*d*node.gamma[i]);
					return solution;
				}else // solution exceeds the present branch
				{
					solution = findsolutionInTheNode(*node.next[i],p); //find solution from the leaf node
					return solution;
				}
			}
		}
		//no optimal solution found in present discussion, the node position is the optimal
		solution.flag = false;
		solution.pos = node.pos;
		solution.val = solution.pos.distance(p)/node.radius;		
		return solution;
	}

	mySolution findsolutionInTheRootNode(myNode node, myPoint p, vector<bool> flag)//find the optimal solution from the root node
	{
		mySolution solution;
		for(int i=0;i<node.N;i++)
		{
			if(isInBranch(p,node,i,flag[i]))//the only diffenrence with findsolutionInTheNode
			{
				myVec po(p, node.pos);
				double kcosbeta = po.dot(node.direction[i]);
				double ksquare = po.length()*po.length();
				double a = node.gamma[i] - kcosbeta;
				double b = node.radius*node.radius - ksquare;
				double c = node.radius*node.radius*kcosbeta - ksquare*node.gamma[i];
				//(gamma-k*(cos(theta))lambda^2 + (r^2-k^2)lambda + (kcos(theta)*r^2 - gamma*k^2) = 0 
				double d = -b/(2*a) + sqrt((b*b-4*a*c)/(4*a*a)); 			
				if(node.delta[i]==-1 || node.delta[i]>d) 
				{
					solution.flag = false;
					solution.pos.setx(node.pos.getx()+node.direction[i].getx()*d);
					solution.pos.sety(node.pos.gety()+node.direction[i].gety()*d);
					solution.val = solution.pos.distance(p)/sqrt(d*d + node.radius*node.radius + 2*d*node.gamma[i]);
					return solution;
				}else 
				{
					solution = findsolutionInTheNode(*node.next[i],p); 
					return solution;
				}
			}
		}
		solution.flag = false;
		solution.pos = node.pos;
		solution.val = solution.pos.distance(p)/node.radius;		
		return solution;
	}

	mySolution findsolution(myPoint p)//find the optimal solution from the tree 
	{
		mySolution solution;
		if(isInConvexHull(p))//if p lies in the convexhull, the optimal solution is infinity and the optimal value is 1
		{
			solution.flag = true;
			solution.pos.setx(0);
			solution.pos.sety(0);
			solution.val = 1;			
		}else//if p lies outside the convexhull, the optimal solution is found from the tree
		{
			solution = findsolutionInTheRootNode(*root,p,flag);
		}
		return solution;
	}

	// void getbound(myPoint &p1, myPoint &p2) //根据root对应的分支确定边界
	// {
		// double xmin=root->pos.getx();
		// double xmax=root->pos.getx();
		// double ymin=root->pos.gety();
		// double ymax=root->pos.gety();			
		// for(int i=0;i<root->N;i++)
		// {
			// double x1,x2,y1,y2;
			// double radius = root->radius*root->radius/(2*root->gamma[i]);
			// x1 = root->pos.getx()-radius*root->direction[i].getx()-radius;
			// x2 = root->pos.getx()-radius*root->direction[i].getx()+radius;
			// y1 = root->pos.gety()-radius*root->direction[i].gety()-radius;
			// y2 = root->pos.gety()-radius*root->direction[i].gety()+radius;
			// xmin = (x1<xmin)? x1:xmin;
			// ymin = (y1<ymin)? y1:ymin;
			// xmax = (x2>xmax)? x2:xmax;
			// ymax = (y2>ymax)? y2:ymax;
		// }
		// p1.setx(xmin);
		// p1.sety(ymin);
		// p2.setx(xmax);
		// p2.sety(ymax);	
	// }
};