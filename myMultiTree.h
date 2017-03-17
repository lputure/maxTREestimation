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
	myPoint pos;			//Բ��λ��
	double radius;			//�뾶
	int N;					//��֧������Ҷ�ӽڵ�����
	vector<double> gamma;	//gamma������Բ�ĵ��ߵľ���
	vector<double> delta;	//delta��������֧����һ�ڵ㵽��ǰ�ڵ��λ��
	vector<myVec> direction;//ָʾ��֧����ĵ�λ����
	vector<myNode*> next;	//��һ�ڵ�
};

class myTree{
private:
	myNode *root;
	vector<myPoint> convexhull;
	vector<bool> flag; //ָʾRootNode��Ӧ��Branch�Ƿ�Ϊֱ��
public:
	myTree(){};
	myTree(vector<myPoint> convexhull, myPoint miniballcenter, double radius)
	{
		this->convexhull = convexhull;
		root = new myNode;
		root->pos = miniballcenter;
		//cout<<root->pos.getx()<<' '<<root->pos.gety()<<endl; //��ʱ����Ļ���
		root->radius = radius;		

		vector<int> supportid; //ȷ��λ����С���Բ�ϵĵ�
		for(unsigned int i=0;i<convexhull.size();i++)
		{
			//if(convexhull[i].distance(miniballcenter)==radius) //�����ϸ����ɸѡ,�д�
			double tdis = convexhull[i].distance(miniballcenter);
			if(tdis<radius+EPS && tdis>radius-EPS) //Ϊ��ֵ�ȶ�
			{
				supportid.push_back(i);
			}
		}
		root->N = supportid.size();

		for(unsigned int i=0;i<supportid.size();i++)//����͹����λ����С���Բ�ĵ���ɵı�
		{
			int id1 = supportid[i];
			int id2 = supportid[(i+1)%supportid.size()];
			myPoint x1 = convexhull[id1];
			myPoint x2 = convexhull[id2];	//ȡ��������

			if(fabs(cross(x1,x2,miniballcenter))<EPS)//x1,x2λ��ֱ����
			{
				this->flag.push_back(true);
				myVec v(miniballcenter,x1);
				myVec vp(v.gety(),-v.getx());
				root->direction.push_back(vp.norm());
				root->gamma.push_back(0);
				if(id1+1==id2 || (id1+1)%convexhull.size()==id2) //��ǰ֧�ŵ�֮��û�������ĵ�
				{
					root->next.push_back(NULL);
					root->delta.push_back(-1); //delta����
				}else //��ǰ֧�ŵ�֮����������ĵ�
				{
					myNode *leaf = buildleaf(root->pos,root->direction[i],id1,id2); //����Ҷ�ӽڵ�
					root->next.push_back(leaf); //�����Ӹ��ڵ㵽Ҷ�ӽڵ����ϵ
					root->delta.push_back(miniballcenter.distance(leaf->pos)); //delta����
				}
			}else //һ������
			{
				this->flag.push_back(false);
				myVec v1(x1,miniballcenter);
				myVec v2(x2,miniballcenter);
				myVec v = v1 + v2;
				myVec dir = v.norm();
				root->direction.push_back(dir); //���ߵķ���
				root->gamma.push_back(miniballcenter.distance((x1+x2)/2)); //gamma����
				if(id1+1==id2 || (id1+1)%convexhull.size()==id2) //��ǰ֧�ŵ�֮��û�������ĵ�
				{
					root->next.push_back(NULL);
					root->delta.push_back(-1); //delta����
				}else //��ǰ֧�ŵ�֮����������ĵ�
				{
					myNode *leaf = buildleaf(root->pos,root->direction[i],id1,id2); //����Ҷ�ӽڵ�
					root->next.push_back(leaf); //�����Ӹ��ڵ㵽Ҷ�ӽڵ����ϵ
					root->delta.push_back(miniballcenter.distance(leaf->pos)); //delta����
				}			
			}
		}
	}

	myNode *buildleaf(myPoint parentnodeposition, myVec v, int j1, int j2)
	{
		myNode *leafnode = new myNode;

		if(j2<j1) j2 = j2 + convexhull.size(); //У���ڶ����ڵ�

		vector<double> tdelta; //ȷ��j1,j2֮��ĵ������߷�����ȷ���Ĺյ�
		for(int i=j1+1;i<j2;i++)
		{
			myVec y1(convexhull[j1],parentnodeposition); //j1��
			myVec y3(convexhull[i%convexhull.size()],parentnodeposition); //j1,j2֮���i��
				// from |y1-tmp*v| = |y3-tmp*v| to tmp = (|y3|*|y3|-|y1|*|y1|)/(2v(y1-y3))
			double tmp = (y3.length()*y3.length()-y1.length()*y1.length())/(2*v.dot(y1-y3));
			tdelta.push_back(tmp);
		}
		double vmin = tdelta[0]; //����СֵΪ�յ�
		if(tdelta.size()>1)
		{
			for(unsigned int i=1;i<tdelta.size();i++)
				if(tdelta[i]<vmin)
					vmin = tdelta[i];
		}

		leafnode->pos = myPoint(parentnodeposition.getx()+v.getx()*vmin,parentnodeposition.gety()+v.gety()*vmin);
		//cout<<leafnode->pos.getx()<<' '<<leafnode->pos.gety()<<endl;
		leafnode->radius = leafnode->pos.distance(convexhull[j1]);
		
		tdelta.push_back(vmin);//��j2���ֵ���������
		int pid = j1; //pidָ���һ����
		int cnt = 0;
		for(unsigned int i=0;i<tdelta.size();i++)
		{			
 			if(tdelta[i]<vmin+EPS && tdelta[i]>vmin-EPS)//Ϊ��ֵ�ȶ�
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
				leafnode->gamma.push_back(leafnode->pos.distance((x1+x2)/2)); //��root�ڵ���ͬ�Ľ�������
				if(id1+1==id2 || (id1+1)%convexhull.size()==id2)
				{
					leafnode->next.push_back(NULL);
					leafnode->delta.push_back(-1);
				}else
				{
					//myNode *nextleaf = buildleaf(leafnode->pos,leafnode->direction[i],id1,id2); //�ݹ����Ҷ�ӽڵ�Ľ���
					myNode *nextleaf = buildleaf(leafnode->pos,leafnode->direction[cnt-1],id1,id2); //�ݹ����Ҷ�ӽڵ�Ľ���
					leafnode->next.push_back(nextleaf); //���鸸�ڵ㵽�ӽڵ����ϵ
					leafnode->delta.push_back(leafnode->pos.distance(nextleaf->pos)); //delta����
				}
				pid = id2; //����pid
			}
		}
		leafnode->N = leafnode->next.size();
		return leafnode;
	}

	double cross(myPoint p0, myPoint p1, myPoint p2)
	{
		return (p1.getx()-p0.getx())*(p2.gety()-p0.gety())-(p2.getx()-p0.getx())*(p1.gety()-p0.gety());
	}

	bool isInConvexHull(myPoint p) //�жϵ�p�Ƿ���͹���ڣ�O(N)���Ӷ�
	{
		int n = convexhull.size();
		if(n<3) return false;
		for(int i=0;i<n;i++)//convex hull������������
		{
			int id1 = i;
			int id2 = (i+1)%n;
			if(cross(convexhull[id1],p,convexhull[id2])>0)
				return false;
		}
		return true;
	}

	bool isInBranch(myPoint p,myNode node, int i) //�жϵ�p�Ƿ��ڽڵ�node�ĵ�i����֧
		// �����������������Լ��Բ֮�� & ͹��֮��
	{
		double radius = node.radius*node.radius/2/node.gamma[i]; // r'=r*r/(2*gamma)
		myPoint center;
		center.setx(node.pos.getx()-radius*node.direction[i].getx());
		center.sety(node.pos.gety()-radius*node.direction[i].gety());
		//��һ���ж�����������֮�⣨Ҫ֮Ҫ֮��
		myVec po(p, node.pos);
		return center.distance(p)<radius && node.gamma[i]<po.dot(node.direction[i]);
	}
	 
	mySolution findsolutionInTheNode(myNode node, myPoint p)//�����ṹ�Ľڵ���Ѱ�����Ž�
	{
		mySolution solution;
		for(int i=0;i<node.N;i++)//�������еķ�֧
		{
			if(isInBranch(p,node,i))//�ڵ�i����֧��
			{
				myVec po(p, node.pos);
				double kcosbeta = po.dot(node.direction[i]);
				double ksquare = po.length()*po.length();
				double a = node.gamma[i] - kcosbeta;
				double b = node.radius*node.radius - ksquare;
				double c = node.radius*node.radius*kcosbeta - ksquare*node.gamma[i]; //���η���
					//(gamma-k*(cos(theta))lambda^2 + (r^2-k^2)lambda + (kcos(theta)*r^2 - gamma*k^2) = 0 
				double d = -b/(2*a) + sqrt((b*b-4*a*c)/(4*a*a)); //ȡ������Ľ�
				//if(d<0) //��Щ����£�d��ֵΪ��
				//{
				//	cout<<i<<' '<<p.getx()<<' '<<p.gety()<<' '<<a<<' '<<b<<' '<<c<<' '<<d<<endl;
				//}				
				if(node.delta[i]==-1 || node.delta[i]>d) //���ڵ�ǰ��֧
				{
					solution.flag = false;
					solution.pos.setx(node.pos.getx()+node.direction[i].getx()*d);
					solution.pos.sety(node.pos.gety()+node.direction[i].gety()*d);
					solution.val = solution.pos.distance(p)/sqrt(d*d + node.radius*node.radius + 2*d*node.gamma[i]);
					return solution;
				}else //�ⳬ����ǰ��֧
				{
					solution = findsolutionInTheNode(*node.next[i],p); //������һ��֧�����Ž�
					return solution;
				}
			}
		}
		//֮ǰ�ı�����û�н������ô��ǰ�ڵ�λ��Ϊ���Ž�
		solution.flag = false;
		solution.pos = node.pos;
		solution.val = solution.pos.distance(p)/node.radius;		
		return solution;
	}

	bool isInBranch(myPoint p,myNode node, int i,bool flag) //�жϵ�p�Ƿ��ڸ��ڵ�ĵ�i����֧
		// �����������������Լ��Բ֮�� & ͹��֮��
	{
		if(flag==false)
		{
			double radius = node.radius*node.radius/2/node.gamma[i]; // r'=r*r/(2*gamma)
			myPoint center;
			center.setx(node.pos.getx()-radius*node.direction[i].getx());
			center.sety(node.pos.gety()-radius*node.direction[i].gety());
			//��һ���ж�����������֮�⣨Ҫ֮Ҫ֮��
			myVec po(p, node.pos);
			return center.distance(p)<radius && node.gamma[i]<po.dot(node.direction[i]);
		}else
		{
			myVec po(p, node.pos);
			return node.gamma[i]<po.dot(node.direction[i]);
		}
	}

	mySolution findsolutionInTheRootNode(myNode node, myPoint p, vector<bool> flag)//�����ṹ�Ľڵ���Ѱ�����Ž�
	{
		mySolution solution;
		for(int i=0;i<node.N;i++)//�������еķ�֧
		{
			if(isInBranch(p,node,i,flag[i]))//�ڵ�i����֧��
			{
				myVec po(p, node.pos);
				double kcosbeta = po.dot(node.direction[i]);
				double ksquare = po.length()*po.length();
				double a = node.gamma[i] - kcosbeta;
				double b = node.radius*node.radius - ksquare;
				double c = node.radius*node.radius*kcosbeta - ksquare*node.gamma[i]; //���η���
				//(gamma-k*(cos(theta))lambda^2 + (r^2-k^2)lambda + (kcos(theta)*r^2 - gamma*k^2) = 0 
				double d = -b/(2*a) + sqrt((b*b-4*a*c)/(4*a*a)); //ȡ������Ľ�
				//if(d<0) //��Щ����£�d��ֵΪ��
				//{
				//	cout<<i<<' '<<p.getx()<<' '<<p.gety()<<' '<<a<<' '<<b<<' '<<c<<' '<<d<<endl;
				//}				
				if(node.delta[i]==-1 || node.delta[i]>d) //���ڵ�ǰ��֧
				{
					solution.flag = false;
					solution.pos.setx(node.pos.getx()+node.direction[i].getx()*d);
					solution.pos.sety(node.pos.gety()+node.direction[i].gety()*d);
					solution.val = solution.pos.distance(p)/sqrt(d*d + node.radius*node.radius + 2*d*node.gamma[i]);
					return solution;
				}else //�ⳬ����ǰ��֧
				{
					solution = findsolutionInTheNode(*node.next[i],p); //������һ��֧�����Ž�
					return solution;
				}
			}
		}
		//֮ǰ�ı�����û�н������ô��ǰ�ڵ�λ��Ϊ���Ž�
		solution.flag = false;
		solution.pos = node.pos;
		solution.val = solution.pos.distance(p)/node.radius;		
		return solution;
	}

	mySolution findsolution(myPoint p)//�������ṹ��ȷ�����Ž�����Ž��ֵ
	{
		mySolution solution;
		if(isInConvexHull(p))//p��͹���ڣ����Ž�������Զ������ֵΪ1
		{
			solution.flag = true;
			solution.pos.setx(0);
			solution.pos.sety(0);
			solution.val = 1;			
		}else//p��͹���⣬�ӽ�����ṹ��Ѱ�����Ž�
		{
			solution = findsolutionInTheRootNode(*root,p,flag);
		}
		return solution;
	}

	void getbound(myPoint &p1, myPoint &p2) //����root��Ӧ�ķ�֧ȷ���߽�
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