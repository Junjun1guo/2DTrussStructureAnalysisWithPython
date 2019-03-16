2DTrussStructureAnalysisWithPython
==================================
force analysis for 2d truss structure based on the object-oriented programming of Python language



* [ViewTruss2D.py](#viewtruss2dpy)

![](https://github.com/junjun1guo/2DTrussStructureAnalysisWithPython/raw/Truss2DModelAndView/view2.png)  

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
```
ViewTruss2D.py
---------------
```python
#-*-coding: UTF-8-*-
# -*- coding: utf-8 -*- 
import wx
import wx.xrc
import shelve
import os
import numpy as np
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from wx import glcanvas
import sys
import math
from modelTruss2D import BaseFEM,Truss2D


################################################################################
################################################################################
class MyCanvasBase(glcanvas.GLCanvas):
    def __init__(self, parent):
        glcanvas.GLCanvas.__init__(self, parent, -1)
        self.init = True
        self.context = glcanvas.GLContext(self)
        
        # initial mouse position
        self.lastx = self.x = 30
        self.lasty = self.y = 30
        self.size = None
        self.Bind(wx.EVT_SIZE, self.OnSize)
        self.Bind(wx.EVT_PAINT, self.OnPaint)

        self.Bind(wx.EVT_LEFT_DOWN, self.OnMouseDown)
        self.Bind(wx.EVT_LEFT_UP, self.OnMouseUp)
        self.Bind(wx.EVT_MOTION, self.OnMouseMotion)

        self.arrayLenth=0.06
        self.arrayDeg=30
        self.arrayRatio=0.5
        self.circleRadius=0.015
        self.circleSegNum=100

#################################################
    def OnSize(self,event):
        wx.CallAfter(self.DoSetViewport)
        event.Skip()
#################################################
    def DoSetViewport(self):
        size = self.GetClientSize()
        self.Size=size
        self.SetCurrent(self.context)
        glViewport(0, 0, size.width, size.height)
#################################################
    def OnPaint(self, event):

        dc = wx.PaintDC(self)
        self.SetCurrent(self.context)
        if not self.init:
            self.InitGL()
            self.init = True
        self.OnDraw()
#################################################
    def OnMouseDown(self, evt):
        self.CaptureMouse()
        self.x, self.y = self.lastx, self.lasty = evt.GetPosition()
#################################################
    def OnMouseUp(self, evt):
        self.ReleaseMouse()
#################################################
    def OnMouseMotion(self, evt):
		pass
#        if evt.Dragging() and evt.LeftIsDown():
#            self.lastx, self.lasty = self.x, self.y
#            self.x, self.y = evt.GetPosition()
#            self.Refresh(False)
################################################################################
################################################################################
class CubeCanvas(MyCanvasBase):
    def InitGL(self):
        # set viewing projection
        glMatrixMode(GL_PROJECTION)
        glFrustum(-0.5, 0.5, -0.5, 0.5, 1.0, 3.0)
        # position viewer
        glMatrixMode(GL_MODELVIEW)
        glTranslatef(0.0, 0.0, -2.0)
        # position object
        glRotatef(self.y, 1.0, 0.0, 0.0)
        glRotatef(self.x, 0.0, 1.0, 0.0)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
################################################
    def dataProject (self):
        print loadDB.get("nodeDict")
################################################
    def drawPoints (self,loadDB):
        nodeDict=loadDB.get("nodeDict")
        keyList=nodeDict.keys()
        xValueList=[nodeDict[key][0] for key in keyList]
        yValueList=[nodeDict[key][1] for key in keyList]
        xmin=min(xValueList)
        xmax=max(xValueList)
        ymin=min(yValueList)
        ymax=max(yValueList)
        self.normNodeDict={key:((nodeDict[key][0]*2-(xmin+xmax))/float\
			(xmax-xmin)*0.7,(nodeDict[key][1]*2-(ymin+ymax))/float(ymax-ymin)\
			*0.7) for key in keyList}

        glPointSize(5)
        glBegin(GL_POINTS)
        glColor(1, 0, 0)
        for key in keyList:
            glVertex2f(self.normNodeDict[key][0], self.normNodeDict[key][1])
        glEnd()
        glFlush()
        self.SwapBuffers()
        self.SwapBuffers()
################################################
    def OnDraw(self,*args, **kwargs):
        # clear color and depth buffers
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glFinish()
        self.SwapBuffers()
#################################################
    def textDraw (self,color,position,textString):
		glColor(color[0], color[1], color[2])
		glRasterPos2f(position[0],position[1])
		for eachChar in textString:
			glutBitmapCharacter(GLUT_BITMAP_9_BY_15,ord(eachChar))
#################################################
    def nodeNumerDraw (self):
		nodeList=self.normNodeDict.keys()
		textcolor=(0, 1, 0)
		for node in nodeList:
			textposition=(self.normNodeDict[node][0]+0.02,\
				self.normNodeDict[node][1]+0.02)
			textString=str(int(node))
			self.textDraw(textcolor,textposition,textString)
		glFlush()
		self.SwapBuffers()
#################################################
    def eleDraw (self,loadDB):
		eleDict=loadDB.get("eleDict")
		keyList=eleDict.keys()
		self.SwapBuffers()
		glBegin( GL_LINES)
		glColor(0, 0, 1)
		for key in keyList:
			I=eleDict[key][0]
			J=eleDict[key][1]
			glVertex2f(self.normNodeDict[I][0],self.normNodeDict[I][1])
			glVertex2f(self.normNodeDict[J][0],self.normNodeDict[J][1])
		glEnd()
		glFlush()
		self.SwapBuffers()
#################################################
    def eleNumberDraw (self,loadDB):
		eleDict=loadDB.get("eleDict")
		nodeDict=self.normNodeDict
		eleNumList=eleDict.keys()
		textColor=(1,0,0)
		self.SwapBuffers()
		for eleNum in eleNumList:
			nodeI=eleDict[eleNum][0]
			nodeJ=eleDict[eleNum][1]
			nodeIx=nodeDict[nodeI][0]
			nodeIy=nodeDict[nodeI][1]
			nodeJx=nodeDict[nodeJ][0]
			nodeJy=nodeDict[nodeJ][1]
			textPosition=((nodeIx+nodeJx)*0.5,(nodeIy+nodeJy)*0.5)
			textString=str(int(eleNum))
			self.textDraw(textColor,textPosition,textString)
		
		glFlush()
		self.SwapBuffers()	
#################################################
    def lineDrawPoints (self,textColor,pointsList):
		glBegin( GL_LINES)
		glColor(textColor[0], textColor[1], textColor[2])
		for i1 in range(len(pointsList)):
			glVertex2f(pointsList[i1][0],pointsList[i1][1])
		glEnd()
#################################################
    def xDirArrayPositive (self,textColor,nodeX,nodeY):
		arrayEndX=nodeX-self.arrayLenth
		arrayEndY=nodeY
		beta=self.arrayDeg*math.pi/float(180)
		r=self.arrayLenth*self.arrayRatio
		upX=nodeX-r*math.cos(beta)
		upY=nodeY+r*math.sin(beta)
		downX=upX
		downY=nodeY-r*math.sin(beta)
		pointList=[(nodeX,nodeY),(arrayEndX,arrayEndY),(nodeX,nodeY),\
			(upX,upY),(nodeX,nodeY),(downX,downY)]
		self.lineDrawPoints(textColor,pointList)
		return (arrayEndX,arrayEndY)
#################################################
    def xDirArrayNegtive (self,textColor,nodeX,nodeY):
		arrayEndX=nodeX+self.arrayLenth
		arrayEndY=nodeY
		beta=self.arrayDeg*math.pi/float(180)
		r=self.arrayLenth*self.arrayRatio
		upX=nodeX+r*math.cos(beta)
		upY=nodeY+r*math.sin(beta)
		downX=upX
		downY=nodeY-r*math.sin(beta)
		pointList=[(nodeX,nodeY),(arrayEndX,arrayEndY),(nodeX,nodeY),\
			(upX,upY),(nodeX,nodeY),(downX,downY)]
		self.lineDrawPoints(textColor,pointList)
		return (arrayEndX,arrayEndY)
#################################################
    def yDirArrayPostive (self,textColor,nodeX,nodeY):
		arrayEndX=nodeX
		arrayEndY=nodeY-self.arrayLenth
		beta=self.arrayDeg*math.pi/float(180)
		r=self.arrayLenth*self.arrayRatio
		leftX=nodeX-r*math.sin(beta)
		leftY=nodeY-r*math.cos(beta)
		rightX=nodeX+r*math.sin(beta)
		rightY=leftY
		pointList=[(nodeX,nodeY),(arrayEndX,arrayEndY),(nodeX,nodeY),\
			(leftX,leftY),(nodeX,nodeY),(rightX,rightY)]
		self.lineDrawPoints(textColor,pointList)
		return (arrayEndX,arrayEndY)
#################################################
    def yDirArrayNegtive (self,textColor,nodeX,nodeY):
		arrayEndX=nodeX
		arrayEndY=nodeY+self.arrayLenth
		beta=self.arrayDeg*math.pi/float(180)
		r=self.arrayLenth*self.arrayRatio
		leftX=nodeX-r*math.sin(beta)
		leftY=nodeY+r*math.cos(beta)
		rightX=nodeX+r*math.sin(beta)
		rightY=leftY
		pointList=[(nodeX,nodeY),(arrayEndX,arrayEndY),(nodeX,nodeY),\
			(leftX,leftY),(nodeX,nodeY),(rightX,rightY)]
		self.lineDrawPoints(textColor,pointList)
		return (arrayEndX,arrayEndY)
#################################################
    def loadDraw (self,loadDB):
		loadDict=loadDB.get("loadDict")
		nodeDict=self.normNodeDict
		loadNumList=loadDict.keys()
		textColory=[x/float(255) for x in (255,106,106)]
		textColorx=[x/float(255) for x in (255,255,0)]
		self.SwapBuffers()
		for loadNode in loadNumList:
			if loadDict[loadNode][0]>0:
				position=self.xDirArrayPositive(textColorx,\
					nodeDict[loadNode][0],nodeDict[loadNode][1])
				textString=str(np.abs(loadDict[loadNode][0]))
				self.textDraw (textColorx,position,textString)
			if loadDict[loadNode][0]<0:
				position=self.xDirArrayNegtive(textColorx,\
					nodeDict[loadNode][0],nodeDict[loadNode][1])
				textString=str(np.abs(loadDict[loadNode][0]))
				self.textDraw (textColorx,position,textString)
			if loadDict[loadNode][1]>0:
				position=self.yDirArrayPostive(textColory,\
					nodeDict[loadNode][0],nodeDict[loadNode][1])
				textString=str(np.abs(loadDict[loadNode][1]))
				self.textDraw (textColory,position,textString)
			if loadDict[loadNode][1]<0:
				position=self.yDirArrayNegtive(textColory,\
					nodeDict[loadNode][0],nodeDict[loadNode][1])
				textString=str(np.abs(loadDict[loadNode][1]))
				self.textDraw (textColory,position,textString)
		glFlush()
		self.SwapBuffers()
#################################################
    def circleDraw (self,color,nodeX,nodeY):
		R=self.circleRadius
		n=self.circleSegNum
		glBegin(GL_LINE_LOOP)
		glColor(color[0], color[1], color[2])
		for i in range(n):
			glVertex2f(nodeX+R*math.cos(2*math.pi/float(n)*i),\
				nodeY+R*math.sin(2*math.pi/float(n)*i))		
		glEnd()
#################################################
    def constraintDraw(self,loadDB):
		constraintDict=loadDB.get("ConstraintDict")
		nodeDict=self.normNodeDict
		loadNumList=constraintDict.keys()
		textColor=[0,1,0]
		self.SwapBuffers()
		for nodeCons in loadNumList:
			self.circleDraw (textColor,nodeDict[nodeCons][0],\
				nodeDict[nodeCons][1])
		glFlush()
		self.SwapBuffers()
			
		
################################################################################
################################################################################
class baseViewFEM(wx.Frame):
	def __init__ (self):
		self.screenSize=wx.DisplaySize()
		title="Finite element analysis for structure"
		sizeScreen=(self.screenSize[0],self.screenSize[1])
		wx.Frame.__init__(self,None,-1,title,size=sizeScreen)
		
		self.UI=self.initUI()
		self.statusBarSet()

		self.layoutFrame()
		self.displayPanel()

		self.processPanel()
		self.layoutProcessPanel()
		self.preprocessBoxSizer()
		self.runButton()
		self.postprocessBoxSizer()

		self.processPanel.SetSizer(self.layoutProcessPanel)
		self.processPanel.Layout()

		self.SetSizer( self.layoutFrame )
		self.Layout()
		
		self.sizer = wx.BoxSizer(wx.HORIZONTAL)
		self.canvas1 = CubeCanvas(self.displayPanel)
		self.sizer.Add(self.canvas1, 1, wx.LEFT | wx.TOP | wx.GROW)
		self.displayPanel.SetSizer(self.sizer)
#################################################
	def initUI (self):
		
		self.icon = wx.Icon('FrameIco.ico', wx.BITMAP_TYPE_ICO)
		self.SetIcon(self.icon) 
#################################################
	def statusBarSet (self):

		self.statusBar=self.CreateStatusBar(number=3)
		statusBar1=self.statusBar.SetStatusText("Version: 1.0.0",0)
		statusBar2=self.statusBar.SetStatusText("Author: Junjun Guo",1)
		statusBar3=self.statusBar.SetStatusText("Email: guojj@tongji.edu.cn",2)
#################################################
	def layoutFrame (self):
		
		self.layoutFrame = wx.FlexGridSizer( 1, 2, 10, 10 )
		self.layoutFrame.AddGrowableCol( 0 )
		self.layoutFrame.AddGrowableRow( 0 )
		self.layoutFrame.SetFlexibleDirection( wx.VERTICAL )
#################################################
	def displayPanel (self):

		self.displayPanel = wx.Panel( self)
		self.layoutFrame.Add( self.displayPanel, 1, wx.EXPAND |wx.ALL, 5 )
#################################################
	def processPanel(self):

		self.processPanel = wx.Panel( self, wx.ID_ANY, wx.DefaultPosition, \
			wx.Size( 200,-1 ),wx.TAB_TRAVERSAL )
		self.layoutFrame.Add( self.processPanel, 1, wx.ALL|wx.EXPAND, 5 )
#################################################
	def layoutProcessPanel(self):

		self.layoutProcessPanel = wx.FlexGridSizer( 3, 1, 10, 10 )
		self.layoutProcessPanel.SetFlexibleDirection( wx.VERTICAL )
		self.layoutProcessPanel.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )
#################################################
	def preprocessBoxSizer (self):

		self.preProcess = wx.StaticBoxSizer( wx.StaticBox( self.processPanel\
			,wx.ID_ANY,"preProcess" ), wx.VERTICAL )
		self.preProcess.SetMinSize( wx.Size( -1,100)) 
		self.layoutProcessPanel.Add( self.preProcess, 1, wx.EXPAND, 5 )
#################################################
	def runButton (self):

		self.runButton = wx.Button( self.processPanel, wx.ID_ANY, "Run",\
			wx.DefaultPosition, wx.DefaultSize, 0 )
		self.runButton.Enable( False )
		self.runButton.SetBackgroundColour( wx.Colour( 0, 255, 255 ) )
		self.layoutProcessPanel.Add( self.runButton, 0, wx.ALIGN_CENTER|\
			wx.ALL, 5 )
#################################################
	def postprocessBoxSizer (self):

		self.postProcess = wx.StaticBoxSizer( wx.StaticBox( self.processPanel\
			, wx.ID_ANY,"postProcess" ), wx.VERTICAL )
		self.layoutProcessPanel.Add( self.postProcess, 1, wx.EXPAND, 5 )
################################################################################
################################################################################
class FEM2DTrussView(baseViewFEM):
	def __init__ (self):
		super(FEM2DTrussView,self).__init__()
		self.copyCanvas1=self.canvas1
		self.loadProcess=self.loadData()
		
		self.eventHandle()
###################################################################
	def buttonCreate (self,label,tipString,backColor,sizerS):
		buttonCreat = wx.Button( self.preProcess.GetStaticBox(), wx.ID_ANY, label,\
			wx.DefaultPosition, wx.DefaultSize, 0 )
		buttonCreat.SetForegroundColour( wx.SystemSettings.GetColour(\
			wx.SYS_COLOUR_BTNTEXT ) )
		buttonCreat.SetBackgroundColour( wx.Colour( backColor[0], \
			backColor[1], backColor[2] ) )
		buttonCreat.SetToolTipString( tipString )
		sizerS.Add( buttonCreat, 0, wx.ALIGN_CENTER|wx.ALL, 5 )
		return buttonCreat
###################################################################		
	def loadData (self):
		bSizer1 = wx.BoxSizer( wx.VERTICAL )
		self.Node=self.buttonCreate("Node","node(nodeNumber,x,y)",\
			( 0, 255, 255 ),bSizer1)
		self.Element=self.buttonCreate("element",\
			"element(elementNum,nodeI,nodeJ,area,elasticStiffness,unitWeight)",\
			( 0, 255, 255 ),bSizer1)
		self.Element.Enable( False )
		self.Load=self.buttonCreate("Load","load(nodeNumber,loadValueInGlobalX"\
			+"loadValueInGlobalY)", ( 0, 255, 255 ),bSizer1)
		self.Load.Enable( False )
		self.Constraint=self.buttonCreate("Constraint",\
			"constraint(nodeNumber,xDispInGlobalCoord,yDispInGlobalCoord)",\
			( 0, 255, 255 ),bSizer1)
		self.Constraint.Enable( False )
		self.preProcess.Add( bSizer1, 1, wx.EXPAND, 5 )

		nodeDispSizer1 = wx.FlexGridSizer( 1, 3, 0, 0 )
		nodeDispSizer1.AddGrowableCol( 1 )
		nodeDispSizer1.SetFlexibleDirection( wx.BOTH )
		nodeDispSizer1.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )

		self.NodeStaticText = wx.StaticText( self.postProcess.GetStaticBox(), wx.ID_ANY,\
			"Node:", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.NodeStaticText.Wrap( -1 )
		self.NodeStaticText.Enable( False )
		nodeDispSizer1.Add( self.NodeStaticText, 0, wx.ALIGN_CENTER|wx.ALL, 1 )

		self.NodeCtrText= wx.TextCtrl( self.postProcess.GetStaticBox(), \
			wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.Size( 80,-1 ), 0 )
		self.NodeCtrText.Enable( False )
		nodeDispSizer1.Add( self.NodeCtrText, 0, wx.ALIGN_CENTER|wx.ALL, 1 )

		self.NodeDispButton = wx.Button( self.postProcess.GetStaticBox(), wx.ID_ANY,"Ok",\
			wx.DefaultPosition, wx.Size( 40,-1 ), 0 )
		self.NodeDispButton.SetBackgroundColour( wx.Colour( 0, 255, 255 ) )
		self.NodeDispButton.Enable( False )
		nodeDispSizer1.Add( self.NodeDispButton, 0, wx.ALIGN_CENTER|wx.ALL, 0 )
		self.postProcess.Add( nodeDispSizer1, 1, wx.EXPAND, 5 )

		nodeDispSizer2 = wx.FlexGridSizer( 1, 2, 0, 0 )
		nodeDispSizer2.SetFlexibleDirection( wx.BOTH )
		nodeDispSizer2.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )

		self.NodeDispStaticText = wx.StaticText(self.postProcess.GetStaticBox(),\
			wx.ID_ANY,"Disp:", wx.DefaultPosition, wx.Size( -1,-1 ), 0 )
		self.NodeDispStaticText.Wrap( -1 )
		self.NodeDispStaticText.Enable( False )
		nodeDispSizer2.Add(self.NodeDispStaticText, 0, wx.ALIGN_CENTER|wx.ALL, 1)

		self.NodeDispCtrlText = wx.TextCtrl( self.postProcess.GetStaticBox(), \
			wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.Size( 150,-1 ), 0 )
		self.NodeDispCtrlText.Enable( False )

		nodeDispSizer2.Add(self.NodeDispCtrlText, 0, wx.ALIGN_CENTER|wx.ALL, 5 )
		self.postProcess.Add(nodeDispSizer2, 1, wx.EXPAND, 5 )

		self.StaticLine= wx.StaticLine( self.postProcess.GetStaticBox(),\
			wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LI_HORIZONTAL )
		self.postProcess.Add( self.StaticLine, 0, wx.EXPAND |wx.ALL, 5 )

		eleForceSizer1 = wx.FlexGridSizer( 1, 3, 0, 0 )
		eleForceSizer1.SetFlexibleDirection( wx.BOTH )
		eleForceSizer1.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )

		self.eleStaticText = wx.StaticText( self.postProcess.GetStaticBox(),\
			wx.ID_ANY,"Element:", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.eleStaticText.Wrap( -1 )
		self.eleStaticText.Enable( False )
		eleForceSizer1.Add( self.eleStaticText, 0, wx.ALIGN_CENTER|wx.ALL, 1 )

		self.eleCtrText = wx.TextCtrl( self.postProcess.GetStaticBox(),\
			wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.Size( 80,-1 ), 0 )
		self.eleCtrText.Enable( False )
		eleForceSizer1.Add( self.eleCtrText, 0, wx.ALIGN_CENTER|wx.ALL, 1 )

		self.eleForceButton = wx.Button(self.postProcess.GetStaticBox(),\
			wx.ID_ANY, "Ok", wx.DefaultPosition, wx.Size( 40,-1 ), 0 )
		self.eleForceButton.SetBackgroundColour( wx.Colour( 0, 255, 255 ) )
		self.eleForceButton.Enable( False )
		eleForceSizer1.Add(self.eleForceButton, 0, wx.ALIGN_CENTER|wx.ALL, 5 )
		self.postProcess.Add( eleForceSizer1, 1, wx.EXPAND, 5 )

		eleForceSizer2 = wx.FlexGridSizer( 1, 2, 0, 0 )
		eleForceSizer2.SetFlexibleDirection( wx.BOTH )
		eleForceSizer2.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )

		self.eleForceStaticText = wx.StaticText( self.postProcess.GetStaticBox(),\
			wx.ID_ANY,"AxialForce:", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.eleForceStaticText.Wrap( -1 )
		self.eleForceStaticText.Enable( False )
		eleForceSizer2.Add( self.eleForceStaticText, 0, wx.ALIGN_CENTER|wx.ALL, 5 )

		self.eleForceCtrText = wx.TextCtrl(self.postProcess.GetStaticBox(),\
			wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		self.eleForceCtrText.Enable( False )
		eleForceSizer2.Add( self.eleForceCtrText, 0, wx.ALIGN_CENTER|wx.ALL, 5 )
		self.postProcess.Add( eleForceSizer2, 1, wx.EXPAND, 5 )

		imagePanel= wx.Panel( self.postProcess.GetStaticBox(), wx.ID_ANY,\
			wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		self.postProcess.Add( imagePanel, 1, wx.EXPAND |wx.ALL, 5 )
		

		image = wx.Image('FrameIco.ico', wx.BITMAP_TYPE_ICO)
		temp = image.ConvertToBitmap()
		self.bmp = wx.StaticBitmap(parent=imagePanel, bitmap=temp)



###################################################################
	def eventHandle (self):

		self.loadDB = shelve.open('dataStore.db', flag='n')
		self.Node.Bind(wx.EVT_BUTTON, self.nodeButtonClicked)
		self.Element.Bind(wx.EVT_BUTTON, self.eleButtonClicked)
		self.Load.Bind(wx.EVT_BUTTON, self.LoadButtonClicked)
		self.Constraint.Bind(wx.EVT_BUTTON, self.constraintButtonClicked)
		self.runButton.Bind(wx.EVT_BUTTON, self.runButtonClicked)
		self.NodeDispButton.Bind(wx.EVT_BUTTON, self.nodeDispButtonClicked)
		self.eleForceButton.Bind(wx.EVT_BUTTON, self.eleForceButtonClicked)
###################################################################
	def fileLoad (self):
		wildcard = "Text Files (*.txt)|*.txt" 
		dlg = wx.FileDialog(self, "Choose a file", os.getcwd(), "", wildcard,wx.FD_OPEN) 
		if dlg.ShowModal() == wx.ID_OK:
			f = open(dlg.GetPath(), 'r') 
		return f
###################################################################		
	def nodeButtonClicked (self,event):
		f=self.fileLoad()
		self.nodeData=np.loadtxt(f)
		nodeDict=self.arrayToDict(self.nodeData)
		self.loadDB["nodeDict"]=nodeDict
		self.copyCanvas1.drawPoints(self.loadDB)
		self.copyCanvas1.nodeNumerDraw ()
		self.Node.Enable(False)
		self.Element.Enable(True)
###################################################################
	def eleButtonClicked (self,event):
		f=self.fileLoad()
		self.eleData=np.loadtxt(f)
		eleDict=self.arrayToDict(self.eleData)
		self.loadDB["eleDict"]=eleDict
		self.copyCanvas1.eleDraw(self.loadDB)
		self.copyCanvas1.eleNumberDraw(self.loadDB)
		self.Element.Enable(False)
		self.Load.Enable(True)
###################################################################
	def LoadButtonClicked (self,event):
		f=self.fileLoad()
		self.loadData=np.loadtxt(f)
		loadDict=self.arrayToDict(self.loadData)
		self.loadDB["loadDict"]=loadDict
		self.copyCanvas1.loadDraw(self.loadDB)
		self.Load.Enable(False)
		self.Constraint.Enable(True)
###################################################################
	def constraintButtonClicked (self,event):
		f=self.fileLoad()
		self.constraintData=np.loadtxt(f)
		constraintDict=self.arrayToDict(self.constraintData)
		self.loadDB["ConstraintDict"]=constraintDict
		self.copyCanvas1.constraintDraw(self.loadDB)
		self.Constraint.Enable(False)
		self.runButton.Enable( True )
###################################################################
	def runButtonClicked (self,event):
		instance=Truss2D(self.nodeData,self.eleData,\
			self.loadData,self.constraintData)
		self.NodeStaticText.Enable(True)
		self.NodeCtrText.Enable(True)
		self.NodeDispButton.Enable(True)
		self.NodeDispStaticText.Enable(True)
		self.NodeDispCtrlText.Enable(True)
		self.eleStaticText.Enable(True)
		self.eleCtrText.Enable(True)
		self.eleForceButton.Enable(True)
		self.eleForceStaticText.Enable(True)
		self.eleForceCtrText.Enable(True)
		self.resultDB1 = shelve.open('2DTrussDB.db', flag='r')
###################################################################
	def nodeDispButtonClicked (self,event):
		node=int(self.NodeCtrText.GetValue())
		nodeDispMat=self.resultDB1.get("globalNodeDispMat")
		nodeLinkDict=self.resultDB1.get("nodeLinkDict")
		try:
			x=nodeDispMat[nodeLinkDict[node],0]
			y=nodeDispMat[nodeLinkDict[node]+1,0]
			floatX=float('%.6f' % x)
			floatY=float('%.6f' % y)
			nodeValue=(floatX,floatY)
			self.NodeDispCtrlText.SetValue(str(nodeValue))
		except KeyError:
			self.NodeDispCtrlText.SetValue("Node not exists!")
###################################################################
	def eleForceButtonClicked (self,event):
		element=int(self.eleCtrText.GetValue())
		eleForceDict=self.resultDB1.get("axialForceDict")
		try:
			axialForce=eleForceDict[element]
			floatAxialForce=float('%.3f' % axialForce)
			self.eleForceCtrText.SetValue(str(floatAxialForce))
		except KeyError:
			self.eleForceCtrText.SetValue("Ele not exists!")
###################################################################
	def arrayToDict (self,data):
		dictReturn={}
		try:
			m,n=data.shape
			for i in range(m):
				dictReturn[data[i,0]]=(data[i,1:])
		except ValueError:
			dictReturn[data[0]]=(data[1:])
		return dictReturn
################################################################################
################################################################################
if __name__=="__main__":
	app=wx.App()
	frame=FEM2DTrussView()
	frame.Show(True)
	app.MainLoop()
