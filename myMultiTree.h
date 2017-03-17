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
	myPoint pos;			//圆心位置
	double radius;			//半径
	int N;					//分支数量（叶子节点数）
	vector<double> gamma;	//gamma参数，圆心到边的距离
	vector<double> delta;	//delta参数，分支上下一节点到当前节点的位置
	vector<myVec> direction;//指示分支方向的单位向量
	vector<myNode*> next;	//下一节点
};

class myTree{
private:
	myNode *root;
	vector<myPoint> convexhull;
	vector<bool> flag; //指示RootNode对应的Branch是否为直径
public:
	myTree(){};
	myTree(vector<myPoint> convexhull, myPoint miniballcenter, double radius)
	{
		this->convexhull = convexhull;
		root = new myNode;
		root->pos = miniballcenter;
		//cout<<root->pos.getx()<<' '<<root->pos.gety()<<endl; //临时的屏幕输出
		root->radius = radius;		

		vector<int> supportid; //确定位于最小外接圆上的点
		for(unsigned int i=0;i<convexhull.size();i++)
		{
			//if(convexhull[i].distance(miniballcenter)==radius) //先用严格相等筛选,有错
			double tdis = convexhull[i].distance(miniballcenter);
			if(tdis<radius+EPS && tdis>radius-EPS) //为数值稳定
			{
				supportid.push_back(i);
			}
		}
		root->N = supportid.size();

		for(unsigned int i=0;i<supportid.size();i++)//遍历凸包上位于最小外接圆的点组成的边
		{
			int id1 = supportid[i];
			int id2 = supportid[(i+1)%supportid.size()];
			myPoint x1 = convexhull[id1];
			myPoint x2 = convexhull[id2];	//取出两个点

			if(fabs(cross(x1,x2,miniballcenter))<EPS)//x1,x2位于直径上
			{
				this->flag.push_back(true);
				myVec v(miniballcenter,x1);
				myVec vp(v.gety(),-v.getx());
				root->direction.push_back(vp.norm());
				root->gamma.push_back(0);
				if(id1+1==id2 || (id1+1)%convexhull.size()==id2) //当前支撑点之间没有其他的点
				{
					root->next.push_back(NULL);
					root->delta.push_back(-1); //delta参数
				}else //当前支撑点之间存在其他的点
				{
					myNode *leaf = buildleaf(root->pos,root->direction[i],id1,id2); //建立叶子节点
					root->next.push_back(leaf); //建立从根节点到叶子节点的联系
					root->delta.push_back(miniballcenter.distance(leaf->pos)); //delta参数
				}
			}else //一般情形
			{
				this->flag.push_back(false);
				myVec v1(x1,miniballcenter);
				myVec v2(x2,miniballcenter);
				myVec v = v1 + v2;
				myVec dir = v.norm();
				root->direction.push_back(dir); //射线的方向
				root->gamma.push_back(miniballcenter.distance((x1+x2)/2)); //gamma参数
				if(id1+1==id2 || (id1+1)%convexhull.size()==id2) //当前支撑点之间没有其他的点
				{
					root->next.push_back(NULL);
					root->delta.push_back(-1); //delta参数
				}else //当前支撑点之间存在其他的点
				{
					myNode *leaf = buildleaf(root->pos,root->direction[i],id1,id2); //建立叶子节点
					root->next.push_back(leaf); //建立从根节点到叶子节点的联系
					root->delta.push_back(miniballcenter.distance(leaf->pos)); //delta参数
				}			
			}
		}
	}

	myNode *buildleaf(myPoint parentnodeposition, myVec v, int j1, int j2)
	{
		myNode *leafnode = new myNode;

		if(j2<j1) j2 = j2 + convexhull.size(); //校正第二个节点

		vector<double> tdelta; //确定j1,j2之间的点在射线方向上确定的拐点
		for(int i=j1+1;i<j2;i++)
		{
			myVec y1(convexhull[j1],parentnodeposition); //j1点
			myVec y3(convexhull[i%convexhull.size()],parentnodeposition); //j1,j2之间的i点
				// from |y1-tmp*v| = |y3-tmp*v| to tmp = (|y3|*|y3|-|y1|*|y1|)/(2v(y1-y3))
			double tmp = (y3.length()*y3.length()-y1.length()*y1.length())/(2*v.dot(y1-y3));
			tdelta.push_back(tmp);
		}
		double vmin = tdelta[0]; //以最小值为拐点
		if(tdelta.size()>1)
		{
			for(unsigned int i=1;i<tdelta.size();i++)
				if(tdelta[i]<vmin)
					vmin = tdelta[i];
		}

		leafnode->pos = myPoint(parentnodeposition.getx()+v.getx()*vmin,parentnodeposition.gety()+v.gety()*vmin);
		//cout<<leafnode->pos.getx()<<' '<<leafnode->pos.gety()<<endl;
		leafnode->radius = leafnode->pos.distance(convexhull[j1]);
		
		tdelta.push_back(vmin);//将j2点的值加入参数中
		int pid = j1; //pid指向第一个点
		int cnt = 0;
		for(unsigned int i=0;i<tdelta.size();i++)
		{			
 			if(tdelta[i]<vmin+EPS && tdelta[i]>vmin-EPS)//为数值稳定
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
				leafnode->gamma.push_back(leafnode->pos.distance((x1+x2)/2)); //与root节点相同的建立过程
				if(id1+1==id2 || (id1+1)%convexhull.size()==id2)
				{
					leafnode->next.push_back(NULL);
					leafnode->delta.push_back(-1);
				}else
				{
					//myNode *nextleaf = buildleaf(leafnode->pos,leafnode->direction[i],id1,id2); //递归调用叶子节点的建立
					myNode *nextleaf = buildleaf(leafnode->pos,leafnode->direction[cnt-1],id1,id2); //递归调用叶子节点的建立
					leafnode->next.push_back(nextleaf); //建议父节点到子节点的联系
					leafnode->delta.push_back(leafnode->pos.distance(nextleaf->pos)); //delta参数
				}
				pid = id2; //更新pid
			}
		}
		leafnode->N = leafnode->next.size();
		return leafnode;
	}

	double cross(myPoint p0, myPoint p1, myPoint p2)
	{
		return (p1.getx()-p0.getx())*(p2.gety()-p0.gety())-(p2.getx()-p0.getx())*(p1.gety()-p0.gety());
	}

	bool isInConvexHull(myPoint p) //判断点p是否在凸包内，O(N)复杂度
	{
		int n = convexhull.size();
		if(n<3) return false;
		for(int i=0;i<n;i++)//convex hull按右手序排列
		{
			int id1 = i;
			int id2 = (i+1)%n;
			if(cross(convexhull[id1],p,convexhull[id2])>0)
				return false;
		}
		return true;
	}

	bool isInBranch(myPoint p,myNode node, int i) //判断点p是否在节点node的第i个分支
		// 返回真的条件：点在约束圆之内 & 凸包之外
	{
		double radius = node.radius*node.radius/2/node.gamma[i]; // r'=r*r/(2*gamma)
		myPoint center;
		center.setx(node.pos.getx()-radius*node.direction[i].getx());
		center.sety(node.pos.gety()-radius*node.direction[i].gety());
		//多一个判定条件，点在之外（要之要之）
		myVec po(p, node.pos);
		return center.distance(p)<radius && node.gamma[i]<po.dot(node.direction[i]);
	}
	 
	mySolution findsolutionInTheNode(myNode node, myPoint p)//从树结构的节点中寻找最优解
	{
		mySolution solution;
		for(int i=0;i<node.N;i++)//遍历所有的分支
		{
			if(isInBranch(p,node,i))//在第i个分支内
			{
				myVec po(p, node.pos);
				double kcosbeta = po.dot(node.direction[i]);
				double ksquare = po.length()*po.length();
				double a = node.gamma[i] - kcosbeta;
				double b = node.radius*node.radius - ksquare;
				double c = node.radius*node.radius*kcosbeta - ksquare*node.gamma[i]; //二次方程
					//(gamma-k*(cos(theta))lambda^2 + (r^2-k^2)lambda + (kcos(theta)*r^2 - gamma*k^2) = 0 
				double d = -b/(2*a) + sqrt((b*b-4*a*c)/(4*a*a)); //取大于零的解
				//if(d<0) //有些情况下，d的值为负
				//{
				//	cout<<i<<' '<<p.getx()<<' '<<p.gety()<<' '<<a<<' '<<b<<' '<<c<<' '<<d<<endl;
				//}				
				if(node.delta[i]==-1 || node.delta[i]>d) //解在当前分支
				{
					solution.flag = false;
					solution.pos.setx(node.pos.getx()+node.direction[i].getx()*d);
					solution.pos.sety(node.pos.gety()+node.direction[i].gety()*d);
					solution.val = solution.pos.distance(p)/sqrt(d*d + node.radius*node.radius + 2*d*node.gamma[i]);
					return solution;
				}else //解超出当前分支
				{
					solution = findsolutionInTheNode(*node.next[i],p); //调用下一分支找最优解
					return solution;
				}
			}
		}
		//之前的遍历均没有结果，那么当前节点位置为最优解
		solution.flag = false;
		solution.pos = node.pos;
		solution.val = solution.pos.distance(p)/node.radius;		
		return solution;
	}

	bool isInBranch(myPoint p,myNode node, int i,bool flag) //判断点p是否在根节点的第i个分支
		// 返回真的条件：点在约束圆之内 & 凸包之外
	{
		if(flag==false)
		{
			double radius = node.radius*node.radius/2/node.gamma[i]; // r'=r*r/(2*gamma)
			myPoint center;
			center.setx(node.pos.getx()-radius*node.direction[i].getx());
			center.sety(node.pos.gety()-radius*node.direction[i].gety());
			//多一个判定条件，点在之外（要之要之）
			myVec po(p, node.pos);
			return center.distance(p)<radius && node.gamma[i]<po.dot(node.direction[i]);
		}else
		{
			myVec po(p, node.pos);
			return node.gamma[i]<po.dot(node.direction[i]);
		}
	}

	mySolution findsolutionInTheRootNode(myNode node, myPoint p, vector<bool> flag)//从树结构的节点中寻找最优解
	{
		mySolution solution;
		for(int i=0;i<node.N;i++)//遍历所有的分支
		{
			if(isInBranch(p,node,i,flag[i]))//在第i个分支内
			{
				myVec po(p, node.pos);
				double kcosbeta = po.dot(node.direction[i]);
				double ksquare = po.length()*po.length();
				double a = node.gamma[i] - kcosbeta;
				double b = node.radius*node.radius - ksquare;
				double c = node.radius*node.radius*kcosbeta - ksquare*node.gamma[i]; //二次方程
				//(gamma-k*(cos(theta))lambda^2 + (r^2-k^2)lambda + (kcos(theta)*r^2 - gamma*k^2) = 0 
				double d = -b/(2*a) + sqrt((b*b-4*a*c)/(4*a*a)); //取大于零的解
				//if(d<0) //有些情况下，d的值为负
				//{
				//	cout<<i<<' '<<p.getx()<<' '<<p.gety()<<' '<<a<<' '<<b<<' '<<c<<' '<<d<<endl;
				//}				
				if(node.delta[i]==-1 || node.delta[i]>d) //解在当前分支
				{
					solution.flag = false;
					solution.pos.setx(node.pos.getx()+node.direction[i].getx()*d);
					solution.pos.sety(node.pos.gety()+node.direction[i].gety()*d);
					solution.val = solution.pos.distance(p)/sqrt(d*d + node.radius*node.radius + 2*d*node.gamma[i]);
					return solution;
				}else //解超出当前分支
				{
					solution = findsolutionInTheNode(*node.next[i],p); //调用下一分支找最优解
					return solution;
				}
			}
		}
		//之前的遍历均没有结果，那么当前节点位置为最优解
		solution.flag = false;
		solution.pos = node.pos;
		solution.val = solution.pos.distance(p)/node.radius;		
		return solution;
	}

	mySolution findsolution(myPoint p)//遍历树结构，确定最优解和最优解的值
	{
		mySolution solution;
		if(isInConvexHull(p))//p在凸包内，最优解在无穷远，最优值为1
		{
			solution.flag = true;
			solution.pos.setx(0);
			solution.pos.sety(0);
			solution.val = 1;			
		}else//p在凸包外，从解的树结构中寻找最优解
		{
			solution = findsolutionInTheRootNode(*root,p,flag);
		}
		return solution;
	}

	void getbound(myPoint &p1, myPoint &p2) //根据root对应的分支确定边界
	{
		double xmin=root->pos.getx();
		double xmax=root->pos.getx();
		double ymin=root->pos.gety();
		double ymax=root->pos.gety();			
		for(int i=0;i<root->N;i++)
		{
			double x1,x2,y1,y2;
			double radius = root->radius*root->radius/(2*root->gamma[i]);
			x1 = root->pos.getx()-radius*root->direction[i].getx()-radius;
			x2 = root->pos.getx()-radius*root->direction[i].getx()+radius;
			y1 = root->pos.gety()-radius*root->direction[i].gety()-radius;
			y2 = root->pos.gety()-radius*root->direction[i].gety()+radius;
			xmin = (x1<xmin)? x1:xmin;
			ymin = (y1<ymin)? y1:ymin;
			xmax = (x2>xmax)? x2:xmax;
			ymax = (y2>ymax)? y2:ymax;
		}
		p1.setx(xmin);
		p1.sety(ymin);
		p2.setx(xmax);
		p2.sety(ymax);	
	}
};