# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 6.13-4 replay file
# Internal Version: 2014_01_04-02.03.49 126873
# Run by yang on Thu Mar 30 11:07:51 2017
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup

def LocateComma (line):
    CommaList = []
    for i in range (0,len(line)):
        if line[i] == ',':
            CommaList.append(i)
    return CommaList

def LocateLeftBracket (line):
    LBList = []
    for i in range (0,len(line)):
        if line[i] == '[':
            LBList.append(i)
    return LBList  

def LocateRightBracket (line):
    RBList = []
    for i in range (0,len(line)):
        if line[i] == ']':
            RBList.append(i)
    return RBList      

def DecomposeData(line):    
    DataList = []    
    temp_cm_List = LocateComma(line)
    temp_lb_List = LocateLeftBracket(line)
    temp_rb_List = LocateRightBracket(line)
                
    for n in range(0, len(temp_cm_List)):
        temp_X = float(line[temp_lb_List[n]+1:temp_cm_List[n]])
        temp_Y = float(line[temp_cm_List[n]+1:temp_rb_List[n]])
                    
        tempVertex = [temp_X,temp_Y]
        DataList.append(tempVertex)
    return [temp_X,temp_Y]
#==============================================================================
# Here please indicate the exact seedlist file name
#==============================================================================
#InputFileName = "Final_Random_SeedList_2017-04-05_14-11-29"
#InputFileName = "Final_Random_SeedList_2017-04-13_09-39-02"

InputFileList = ['Geo_Lenticular_0.68_Configuration#1',\
                 'Geo_Lenticular_0.68_Configuration#2',\
                 'Geo_Real_WL1_0.68',\
                 'Geo_Real_WL2_0.68_Configuration#1',\
                 'Geo_Real_WL2_0.68_Configuration#2',\
                 'Geo_Rectangular_0.68',\
                 'Gradient_WL1_0.73',\
                 'Gradient_WL2_0.659',\
                 'Square_0.5',\
                 'Square_0.6',\
                 'Square_0.7',\
                 'Square_0.752']



for InputFileName in InputFileList:
    executeOnCaeStartup()
    session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(referenceRepresentation=ON)
    session.viewports['Viewport: 1'].setValues(displayedObject=None)
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',sheetSize=1000.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    
    #s.rectangle(point1=(0.0, 0.0), point2=(7.0, 7.0))
    session.viewports['Viewport: 1'].view.setValues(nearPlane=15.3441, 
        farPlane=24.2539, width=35.6365, height=15.5069, cameraPosition=(-3.42912, 
        4.68285, 19.799), cameraTarget=(-3.42912, 4.68285, 0))
    mdb.models['Model-1'].sketches['__profile__'].sketchOptions.setValues(sheetSize=1000.0, gridSpacing=5, gridAuto=OFF)
    session.viewports['Viewport: 1'].view.fitView()
    
    f=open('%s.dat'%InputFileName,'r')
    
    TLines = []
    TLine = []
    
    line = f.readline()
    while line[0:1] == '#':
        line = f.readline()
    line=f.readline()
    TLines.append(line)
    datalist = []
    while line[0:1] != '*': 
        line = f.readline()
        TLines.append(line)
      
    for j in range (0,len(TLines)-1):
        TLine.append(TLines[j])
    
    for tl in TLine:
        commalist = []
        templist = []
        
        for i in range (0,len(tl)):
            if tl[i] == ",":
                commalist.append(i)
        #print commalist
        
        #print tl[1:commalist[0]]
        templist.append(int(tl[1:commalist[0]]))
        templist.append(float(tl[commalist[0]+2:commalist[1]]))
        templist.append(float(tl[commalist[1]+2:commalist[2]]))
        templist.append(float(tl[commalist[4]+2:commalist[5]]))    
        
        datalist.append(templist)
        #print templist
        
    for dl in datalist:
        #print dl
        s.CircleByCenterPerimeter(center=(dl[1], dl[2]), point1=(dl[1]+dl[3], dl[2]))
    
    mdb.models['Model-1'].sketches.changeKey(fromName='__profile__', 
        toName='%s'%InputFileName)
    line = f.readline()
    #print line
    OuterVertexList = []
    InnerVertexList = []
    while line[0:1] != '*':    
        OuterVertexList.append(DecomposeData(line))
        line = f.readline()
    
    firstpoint = OuterVertexList[0]
    OuterVertexList.append(firstpoint)
    for i in range(0,len(OuterVertexList)-1):
        s.Line(point1=(OuterVertexList[i][0], OuterVertexList[i][1]), point2=(OuterVertexList[i+1][0], OuterVertexList[i+1][1]))
    s.unsetPrimaryObject()
    f.close()