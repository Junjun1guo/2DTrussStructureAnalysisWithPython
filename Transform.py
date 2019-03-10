#-*-coding: UTF-8-*-
##############################################################################
#  Author: Junjun Guo
#  E-mail: guojj@tongji.edu.cn/guojj_ce@163.com
#    Date: 20/02/2019
#  Environemet: Successfully excucted in python 2.7
##############################################################################
import abc
import math
import numpy as np

#转换矩阵基类
class Transform(object):
	@abc.abstractmethod
	def transformMatrix(self):
		return

#2维坐标转换矩阵
class Transform2D(Transform):
	#类初始化，pointI=(xi,yi),pointJ=(xj,yj)
	def __init__ (self,pointI,pointJ):

		self.pointI=pointI
		self.pointJ=pointJ
		self.T=np.zeros((4,4))
	#2维平面坐标转换矩阵
	# T=[[C,S,0,0],[-S,C,0,0],[0,0,C,S],[0,0,-S,C]]
	def transformMatrix(self):
		
		xI=self.pointI[0]
		yI=self.pointI[1]
		xJ=self.pointJ[0]
		yJ=self.pointJ[1]
		vectorxy=(xJ-xI,yJ-yI)
		if vectorxy[0]>0 and vectorxy[1]==0:
			alpha=0
		elif vectorxy[0]>0 and vectorxy[1]>0:
			alpha=math.atan(vectorxy[1]/float(vectorxy[0]))
		elif vectorxy[0]==0 and vectorxy[1]>0:
			alpha=math.pi/2.0
		elif vectorxy[0]<0 and vectorxy[1]>0:
			alpha=math.atan(vectorxy[1]/float(vectorxy[0]))+math.pi
		elif vectorxy[0]<0 and vectorxy[1]==0:
			alpha=math.pi
		elif vectorxy[0]<0 and vectorxy[1]<0:
			alpha=math.atan(vectorxy[1]/float(vectorxy[0]))+math.pi
		elif vectorxy[0]==0 and vectorxy[1]<0:
			alpha=math.pi*1.5
		elif vectorxy[0]>0 and vectorxy[1]<0:
			alpha=math.atan(vectorxy[1]/float(vectorxy[0]))
		cos=math.cos(alpha)
		sin=math.sin(alpha)
		
		self.T[0,0]=cos
		self.T[0,1]=sin
		self.T[1,0]=-sin
		self.T[1,1]=cos
		self.T[2,2]=cos
		self.T[2,3]=sin
		self.T[3,2]=-sin
		self.T[3,3]=cos
		
		return self.T


if __name__=="__main__":
	pointI=(0,0)
	pointJ=(1,-1)
	transInstance=Transform2D(pointI,pointJ)
	T=transInstance.transformMatrix()
	print np.linalg.inv(T)