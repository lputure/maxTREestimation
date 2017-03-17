The codes in this file is for the maximum TRE estimation problem in Image Guided System. 

	Please see https://arxiv.org/pdf/1703.00112.pdf for definition of LBWMEC.

(1) myBoundLBWMEC.m: the matlab function for estimating the maximum TRE.
	

(2) myLBWMEC.mexw64: Mex file which implements the core function myLBWMEC.
	
	(2.1) mytest.h: the class to handle with the estimation
	(2.2) myLBWMEC.cpp: the interface to connect C++ to Matlab
	(2.3) test.m: A test for myLBWMEC function.

(3) other .h .hpp .cpp file: files used to implement myLBWMEC using C++

	(3.1) myVec.h: 		A self defined class to handle with 2D vectors
	(3.2) myPoint.h:	A self defined class to handle with 2D points
	(3.3) mysolution.h:	A self defined class to store the solution of LBWMEC
	(3.3) myNode.h & myMultiTree.h:	The class to store the tree shape feasible solution set of LBWMEC problem
	(3.4) myconvexhull.h:	Functions to compute the convex hull of the given point set
	(3.5) myminiball.h & myminball.cpp: Functions to compute MEC problem

	

