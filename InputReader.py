# -*- coding: utf-8 -*-
def LocateSemicolon (line):
    SemicolonList = []
    for i in range (0,len(line)):
        if line[i] == ';':
            SemicolonList.append(i)
    return SemicolonList
    
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
    return DataList

def DecomposeEdge(line):
    EdgeList = []
    semicolonList = LocateSemicolon(line)
    if len(semicolonList) == 1:
        templine = line[semicolonList[0]+2:-2]
        EdgeList.append(DecomposeData(templine))
    else:
        for n in range(0,len(semicolonList)-1):
            templine = line[semicolonList[n]+2:semicolonList[n+1]-1]
            
            EdgeList.append(DecomposeData(templine))
        templine = line[semicolonList[-1]+2:-2]
        #print templine
        EdgeList.append(DecomposeData(templine))
    return EdgeList
    

def ReadFile(FileName):
    p = open ("%s.idt"%FileName,'r')
    line = p.readline()
# Initialize the default parameters    
    OuterVertexList = []    
    InnerVertexList = []
    SubRegionsList = []
    RigidEdges = []
    FlexibleEdges = []
    Radius = 0.0
    DistributionType = 'Gaussian'
    std_dev = 0.0
    Vf = 0.0
    rgdGamma = 1.05
    flbGamma = 0.5
    clsGamma = 1.05
    cTime = 0.5
    maxIncrement = 20000
    displayInterval = 500
    
    while len(line) != 0:
        if line[0:1] == '*':
            if line[1:-1] == 'OuterVertexList':                
                line = p.readline()
                OuterVertexList = DecomposeData(line)
                #print OuterVertexList
                
            if line[1:-1] == 'InnerVertexList':
                line = p.readline()
                if "Void" in line or "void" in line:
                    InnerVertexList = [[]]
                else:
                    InnerVertexList = DecomposeData(line)                
                #print InnerVertexList
                
            if line[1:-1] == 'Subregion':
                line = p.readline()
                SubRegionsList.append(DecomposeData(line))
                SubRegionsList.append(float(line[0:LocateSemicolon(line)[0]]))
                
                
            if line[1:-1] == 'RigidEdges':
                line = p.readline()
                Symbol = line[0:LocateSemicolon(line)[0]]
                if Symbol == 'A' or Symbol == 'a':
                    RigidEdges.append('A')
                    RigidEdges.append(line[LocateSemicolon(line)[0]+1:-1])
                    #print RigidEdges
                
                    #print RigidEdges
                    #print line
                    line  = p.readline()
                    if line[0:1] != '#':
                        Symbol = line[0:LocateSemicolon(line)[0]]                
                        if Symbol == 'M' or Symbol == 'm':
                            RigidEdges.append('M')
                            if "Void" in line or "void" in line:
                                RigidEdges.append('Void')
                                #print RigidEdges
                            else:
                                tempEdge = DecomposeEdge(line)
                                for te in tempEdge:
                                    RigidEdges.append(te)
                            continue
                #print RigidEdges
                
            if line[1:-1] == 'FlexibleEdges':
                line = p.readline()
                Symbol = line[0:LocateSemicolon(line)[0]]
                if Symbol == 'A' or Symbol == 'a':
                    FlexibleEdges.append('A')
                    if "Void" in line or "void" in line:
                        FlexibleEdges.append('Void')
                        #print RigidEdges
                    else:
                        FlexibleEdges.append(line[LocateSemicolon(line)[0]+1:-1])
                    line  = p.readline()
                    if line[0:1] != '#':
                        Symbol = line[0:LocateSemicolon(line)[0]]
                        if Symbol == 'M' or Symbol == 'm':
                            FlexibleEdges.append('M')
                            if "Void" in line or "void" in line:
                                FlexibleEdges.append('Void')
                                #print FlexibleEdges
                            else:
                                tempEdge = DecomposeEdge(line)
                                #print tempEdge                                
                                FlexibleEdges.append(tempEdge)
                                #print FlexibleEdges
                            continue
                
                #print FlexibleEdges                   
                
            if line[1:-1] == 'Radius':
                line = p.readline()
                Radius = float(line)
                #print Radius
                
            if line[1:-1] == 'DistributionType':
                line = p.readline()
                DistributionType = line[0:LocateSemicolon(line)[0]]
                std_dev  = float(line[LocateSemicolon(line)[0]+1:-1])
                #print Radius
                
            if line[1:-1] == 'Vf':
                line = p.readline()
                Vf = float(line)
                #print Vf
                
            if line[1:-1] == 'rgdGamma':
                line = p.readline()
                rgdGamma = float(line)
                #print rgdGamma
                
            if line[1:-1] == 'flbGamma':
                line = p.readline()
                flbGamma = float(line)
                #print flbGamma
                
            if line[1:-1] == 'clsGamma':
                line = p.readline()                
                clsGamma = float(line)
                #print clsGamma
                
            if line[1:-1] == 'cTime':
                line = p.readline()
                cTime = float(line)
                #print cTime
                
            if line[1:-1] == 'maxIncrement':
                line = p.readline()
                maxIncrement = int(line)
                #print maxIncrement
                
            if line[1:-1] == 'displayInterval':
                line = p.readline()
                displayInterval = int(line)
                #print displayInterval 
                        
        line = p.readline()
    
    SumList = [OuterVertexList,InnerVertexList,SubRegionsList,RigidEdges,FlexibleEdges,Radius,\
               DistributionType, std_dev,Vf,rgdGamma,flbGamma,clsGamma,cTime,maxIncrement,displayInterval]
    
    return SumList        

    
#print ReadFile("Test")
#ReadFile("Irregular_three_zones_vf0.6")


