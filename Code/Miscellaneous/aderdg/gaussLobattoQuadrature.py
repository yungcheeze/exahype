#!/usr/bin/python

import math

def gaulob(xL,xR,degree):
	eps = math.pow(10,-15)
	maxiter = 1000
	stop = False
	iter = 0
	diff = 0

	nodes = [0.] * degree
	weight = [0.] * degree
	nodes_old = [0.] * degree
	
	P = [0.] * (degree * degree)

	for i in range(degree):
		nodes[i] = math.cos(math.pi*i/(degree-1))
	
	while(not stop and iter < maxiter):
		iter = iter + 1
		nodes_old = nodes[:]
		for i in range(degree):
			P[i] = 1
			P[i+degree] = nodes[i]
		for k in range(2,degree):
			for i in range(degree):
				P[i+k*degree] = ((2*k-1)*nodes[i]*P[i+(k-1)*degree] - (k-1)*P[i+(k-2)*degree]) / k
		for i in range(degree):
			nodes[i] = nodes_old[i] - (nodes[i] * P[i+(degree-1)*degree] - P[i+(degree-2)*degree]) / (degree*P[i+(degree-1)*degree])
		if(max([abs(x-y) for x, y in zip(nodes,nodes_old)]) < eps):
			break
		
	for i in range(degree):
		weight[i] = 2. / (degree*(degree-1)*P[i+(degree-1)*degree]*P[i+(degree-1)*degree])
		weight[i] = 0.5*weight[i]*(xR-xL)
		nodes[i] = xL+0.5*(nodes[i]+1)*(xR-xL)
	
	return nodes, weight