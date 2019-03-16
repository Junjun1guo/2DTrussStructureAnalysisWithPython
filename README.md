2DTrussStructureAnalysisWithPython
==================================
force analysis for 2d truss structure based on the object-oriented programming of Python language


ModelTruss2D.py
---------------
```python
#-*-coding: UTF-8-*-
################################################################################
#  Author: Junjun Guo
#  E-mail: guojj@tongji.edu.cn/guojj_ce@163.com
#    Date: 20/02/2019
#  Environemet: Successfully excucted in python 2.7
################################################################################
import numpy as np
import abc
import math
import copy
import shelve
from scipy.integrate import quad
from sympy import *
from compiler.ast import flatten
from Transform import Transform2D
#####################################################
class BaseFEM(object):
	
	#Process raw data
	def dataPreProcess(self,data,column):

		rows=len(flatten(data.tolist()))/column
		return data.reshape((rows,column))
	
	#change array data to dictionary data
	def arrayToDict (self,data):
		
		dictReturn={}
		m,n=data.shape
		for i in range(m):
			dictReturn[data[i,0]]=(data[i,1:])
		return dictReturn
	
	#return shape function 
	@abc.abstractmethod
	def shapeFunction (self):
		return
#################################################################################		
#################################################################################		

class Truss2D(BaseFEM):

	def __init__ (self,points,elements,loads,constraints):
		self.largeNumber=10e10
		proPoints=self.dataPreProcess (points,3)
		self.nodesDict=self.arrayToDict(proPoints)

		proElements=self.dataPreProcess(elements,6)
		self.elementsDict=self.arrayToDict(proElements)
		
		proLoads=self.dataPreProcess(loads,3)
		self.loadsDict=self.arrayToDict(proLoads)
		
		proConstraints=self.dataPreProcess(constraints,3)
		self.constrainsDict=self.arrayToDict(proConstraints)
			
		self.eleTransMatDict=self.eleTransformMat ()
		self.eleStiffGloDict=self.elementStiffnessGlobal ()
		self.nodeLinkDict=self.nodeNumLinkGlobalStiff()
		self.globalStruStiffMat=self.structureStiffness()

		self.nodeWeightMat=self.selfWeight()
		self.totalLoad()

		self.finalStruStiffMat=self.processedStruStiff()
		self.finalLoadMat=self.processedLoad()
		self.globalNodeDisp=self.nodeDispGlobal()
		self.aixalForce=self.axialForceN()
		self.axialStress=self.axialNStress()
		self.resultDB=self.resultSave()
######################################################
	# shape function for 2D truss element
	# [(L-x)/L,x/L]
	def shapeFunction (self,eleNum):
		
		x=symbols("x")
		L=self.eleLength (eleNum)
		shapeFunMat=np.mat([[(L-x)/float(L),0,x/float(L),0],\
			[0,(L-x)/float(L),0,x/float(L)]])
		return shapeFunMat,x
######################################################
	# load density in local x, y direction
	def weightLoadDense (self,eleNum):

		nodeI=self.elementsDict[eleNum][0]
		nodeJ=self.elementsDict[eleNum][1]
		xI=self.nodesDict[nodeI][0]
		yI=self.nodesDict[nodeI][1]
		xJ=self.nodesDict[nodeJ][0]
		yJ=self.nodesDict[nodeJ][1]
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
		qx=self.elementsDict[eleNum][4]*sin
		qy=self.elementsDict[eleNum][4]*cos
		return np.mat([qx,qy]).T
######################################################
	# calculate element length
	def eleLength (self,eleNum):

		nodeI=self.elementsDict[eleNum][0]
		nodeJ=self.elementsDict[eleNum][1]
		nodeIx=self.nodesDict[nodeI][0]
		nodeIy=self.nodesDict[nodeI][1]
		nodeJx=self.nodesDict[nodeJ][0]
		nodeJy=self.nodesDict[nodeJ][1]
		L=np.sqrt((nodeJx-nodeIx)**2+(nodeJy-nodeIy)**2)
		return L
######################################################
	#local stiffness matrix for each element
	def __elementStiffnessMat (self,elementDict):
		
		eleNumber=elementDict.keys()[0]
		dictValue=elementDict.values()[0]
		A=dictValue[2]
		E=dictValue[3]
		L=self.eleLength (eleNumber)
		eleK=np.mat([[1,0,-1,0],[0,0,0,0],[-1,0,1,0],[0,0,0,0]])\
			*(A*E/float(L))
		return eleK
#####################################################	
	# global stiffness matrix for each element
	# kglobal=T.T*kloc*T
	def elementStiffnessGlobal (self):

		eleStiffGloDict={}
		eleKeys=self.elementsDict.keys()
		for i in range(len(eleKeys)):
			T=self.eleTransMatDict[eleKeys[i]]
			eleK=self.__elementStiffnessMat ({eleKeys[i]:self.elementsDict[eleKeys[i]]})
			temp=np.dot(T.T,eleK)
			eleKGlobal=np.dot(temp,T)
			eleStiffGloDict[eleKeys[i]]=eleKGlobal
		return eleStiffGloDict
#####################################################	
	# element transform matrix T
	def eleTransformMat (self):
		
		eleTransMat={}
		eleKeys=self.elementsDict.keys()
		for i1 in range(len(eleKeys)):
			nodeI=self.elementsDict[eleKeys[i1]][0]
			nodeJ=self.elementsDict[eleKeys[i1]][1]
			transInstance=Transform2D(self.nodesDict[nodeI],self.nodesDict[nodeJ])
			T=np.mat(transInstance.transformMatrix())
			eleTransMat[eleKeys[i1]]=T
		return eleTransMat
####################################################
	def nodeNumLinkGlobalStiff (self):
		nodeLink={}
		nodeKeys=self.nodesDict.keys()
		for i1 in range(len(nodeKeys)):
			nodeLink[nodeKeys[i1]]=i1*2
		return nodeLink
####################################################
	# Structureal stiffness matrix
	def structureStiffness(self):
		length=len(self.nodesDict)
		struStiffMat=np.mat(np.zeros((length*2,length*2)))
		eleList=self.elementsDict.keys()
		for i1 in range(len(eleList)):
			nodeI=self.elementsDict[eleList[i1]][0]
			nodeJ=self.elementsDict[eleList[i1]][1]
			nodeLinkI=self.nodeLinkDict[nodeI]
			nodeLinkJ=self.nodeLinkDict[nodeJ]
			eleGlo=self.eleStiffGloDict[eleList[i1]]
			struStiffMat[nodeLinkI:(nodeLinkI+2),nodeLinkI:(nodeLinkI+2)]=eleGlo[0:2,0:2]+\
				struStiffMat[nodeLinkI:(nodeLinkI+2),nodeLinkI:(nodeLinkI+2)]
			struStiffMat[nodeLinkI:(nodeLinkI+2),nodeLinkJ:(nodeLinkJ+2)]=eleGlo[0:2,2:4]+\
				struStiffMat[nodeLinkI:(nodeLinkI+2),nodeLinkJ:(nodeLinkJ+2)]
			struStiffMat[nodeLinkJ:(nodeLinkJ+2),nodeLinkI:(nodeLinkI+2)]=eleGlo[2:4,0:2]+\
				struStiffMat[nodeLinkJ:(nodeLinkJ+2),nodeLinkI:(nodeLinkI+2)]
			struStiffMat[nodeLinkJ:(nodeLinkJ+2),nodeLinkJ:(nodeLinkJ+2)]=eleGlo[2:4,2:4]+\
				struStiffMat[nodeLinkJ:(nodeLinkJ+2),nodeLinkJ:(nodeLinkJ+2)]
		return struStiffMat
###################################################
	#convert selfweight of each element to its nodal force
	def selfWeight (self):
		
		length=len(self.nodesDict)
		nodeForceMat=np.mat(np.zeros((length*2,1)))
		eleList=self.elementsDict.keys()
		for i1 in range(len(eleList)):
			nodeI=self.elementsDict[eleList[i1]][0]
			nodeJ=self.elementsDict[eleList[i1]][1]
			nodeLinkI=self.nodeLinkDict[nodeI]
			nodeLinkJ=self.nodeLinkDict[nodeJ]
			N,x=self.shapeFunction (eleList[i1])
			q=self.weightLoadDense (eleList[i1])
			L=self.eleLength (eleList[i1])
			inteFun=np.dot(N.T,q)
			eleLoadLocalIx=integrate(inteFun[0,0],(x,0,L))
			eleLoadLocalIy=integrate(inteFun[1,0],(x,0,L))
			eleLoadLocalJx=integrate(inteFun[2,0],(x,0,L))
			eleLoadLocalJy=integrate(inteFun[3,0],(x,0,L))
			T=self.eleTransMatDict[eleList[i1]]
			localWeightMat=np.mat([eleLoadLocalIx,eleLoadLocalIy,eleLoadLocalJx,eleLoadLocalJy]).T
			globalWeightMat=np.dot(T.T,localWeightMat)
			nodeForceMat[nodeLinkI:(nodeLinkI+2),0]=nodeForceMat[nodeLinkI:(nodeLinkI+2),0]+\
				globalWeightMat[0:2,0]
			nodeForceMat[nodeLinkJ:(nodeLinkJ+2),0]=nodeForceMat[nodeLinkJ:(nodeLinkJ+2),0]+\
				globalWeightMat[2:4,0]
		return nodeForceMat*-1
####################################################
	# add node loads to global selfweight
	def totalLoad (self):

		loadNodes=self.loadsDict.keys()
		for i1 in range(len(loadNodes)):
			nodeLink=self.nodeLinkDict[loadNodes[i1]]
			self.nodeWeightMat[nodeLink:(nodeLink+2),0]=self.nodeWeightMat[nodeLink:(nodeLink+2),0]+\
				np.mat(self.loadsDict[loadNodes[i1]]).T
		return 
####################################################
	# boundary conditon process for structural stiffness matrix
	def processedStruStiff (self):
		
		constraintKeys=self.constrainsDict.keys()
		globalStructCopy=copy.deepcopy(self.globalStruStiffMat)
		for i in range(len(constraintKeys)):
			indexI=self.nodeLinkDict[constraintKeys[i]]
			globalStructCopy[indexI,indexI]=globalStructCopy[indexI,indexI]*self.largeNumber
			globalStructCopy[(indexI+1),(indexI+1)]=globalStructCopy[(indexI+1),(indexI+1)]\
				*self.largeNumber
		return globalStructCopy
####################################################
	# boundary condition process for total load matrix
	def processedLoad (self):
		constraintKeys=self.constrainsDict.keys()
		totalLoadCopy=copy.deepcopy(self.nodeWeightMat)
		for i in range(len(constraintKeys)):
			indexI=self.nodeLinkDict[constraintKeys[i]]
			totalLoadCopy[indexI,0]=self.globalStruStiffMat[indexI,indexI]*self.largeNumber\
				*self.constrainsDict[constraintKeys[i]][0]
			totalLoadCopy[(indexI+1),0]=self.globalStruStiffMat[(indexI+1),(indexI+1)]*\
				self.largeNumber*self.constrainsDict[constraintKeys[i]][1]
		return totalLoadCopy
####################################################
	# calculate node displacement in global coordinate
	def nodeDispGlobal (self):

		constraintKeys=self.constrainsDict.keys()
		wholeStruStiffness=copy.deepcopy(self.finalStruStiffMat)
		wholeLoad=copy.deepcopy(self.finalLoadMat)
		boundaryList=[]
		for i in range(len(constraintKeys)):
			nodeLinkI=self.nodeLinkDict[constraintKeys[i]]
			if self.constrainsDict[constraintKeys[i]][0]==0:
				boundaryList.append(nodeLinkI)
			else:
				continue
			if self.constrainsDict[constraintKeys[i]][1]==0:
				boundaryList.append(nodeLinkI+1)
			else:
				continue
		bWholeStruStiff=np.delete(wholeStruStiffness,boundaryList,axis=0)
		cWholeStruStiff=np.delete(bWholeStruStiff,boundaryList,axis=1)
		bwholeLoad=np.delete(wholeLoad,boundaryList,axis=0)
		cwholeLoad=np.delete(bwholeLoad,boundaryList,axis=1)
		invStruStiff=np.linalg.inv(cWholeStruStiff)
		nodeDisp=np.dot(invStruStiff,cwholeLoad)
		allDimenList=[x for x in range(len(self.nodesDict.keys())*2)]
		nonZeroNodeDisp=[x for x in allDimenList if x not in boundaryList]
		resultDispMat=np.mat(np.zeros((len(allDimenList),1)))
		for i1 in range(len(nodeDisp)):
			resultDispMat[nonZeroNodeDisp[i1],0]=nodeDisp[i1]
		return resultDispMat
####################################################
	# results display
	def axialForceN (self):
		self.globalNodeForce=np.dot(self.globalStruStiffMat,self.globalNodeDisp)
		eleKeys=self.eleTransMatDict.keys()
		axialForce={}
		localDisp={}

		for i in range(len(eleKeys)):
			partialNodeDisp=[]
			partialForce=[]
			nodeI=self.elementsDict[eleKeys[i]][0]
			nodeJ=self.elementsDict[eleKeys[i]][1]
			LinkI=self.nodeLinkDict[nodeI]
			LinkJ=self.nodeLinkDict[nodeJ]
			partialNodeDisp.append(self.globalNodeDisp[LinkI,0])
			partialNodeDisp.append(self.globalNodeDisp[(LinkI+1),0])
			partialNodeDisp.append(self.globalNodeDisp[LinkJ,0])
			partialNodeDisp.append(self.globalNodeDisp[(LinkJ+1),0])
			localEleDisp=np.dot(self.eleTransMatDict[eleKeys[i]],np.mat(partialNodeDisp).T)
			L=self.eleLength (eleKeys[i])
			N=-1*localEleDisp[0]*self.elementsDict[eleKeys[i]][2]*self.elementsDict[eleKeys[i]][3]\
			/float(L)
			axialForce[eleKeys[i]]=N
		return axialForce
####################################################
	# axial stress of truss element
	def axialNStress (self):
		eleKeys=self.elementsDict.keys()
		axialStressDict={}
		for i in range(len(eleKeys)):
			axialStressDict[eleKeys[i]]=self.aixalForce[eleKeys[i]]\
				/float(self.elementsDict[eleKeys[i]][2])
		return axialStressDict
####################################################
	# save result to shelve database library 
	def resultSave (self):
		resultDB = shelve.open('2DTrussDB.db', flag='n')
		resultDB["TransMatDict"]=self.eleTransMatDict
		resultDB["globalStiffMat"]=self.globalStruStiffMat
		resultDB["nodeLinkDict"]=self.nodeLinkDict
		resultDB["globalNodeDispMat"]=self.globalNodeDisp
		resultDB["eleDict"]=self.elementsDict
		resultDB["nodeDict"]=self.nodesDict
		resultDB["axialForceDict"]=self.aixalForce
		resultDB["axialStressDict"]=self.axialStress
		resultDB.close()
####################################################			
if __name__=="__main__":
	#nodes=(nodeNumber,x,y)
	nodes=np.loadtxt("2DTrussNodes.txt")
	#elements=(elementNumber,I,J,A,E,r)
	elements=np.loadtxt("2DTrussElements.txt")
	#loads=(nodeNumber,globalXValue,globalYValue)
	loads=np.loadtxt("2DTrussLoads.txt")
	#constraints=(nodeNumber,xDisp,yDisp)
	# xDisp=0 means fix
	constraints=np.loadtxt("2DTrussConstraints.txt")
	instance=Truss2D(nodes,elements,loads,constraints)
	resultDB = shelve.open('2DTrussDB.db', flag='n')
