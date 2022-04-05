#Created by Keerthi Krishna PARVATHANENI 10/01/2018
# This is code will generate the random distribution of the circular particles in an n sided polygon.
# Here we use a collision based event model described in wrok of 
# " https://torquatocpanel.deptcpanel.princeton.edu/wp-content/uploads/files/publications/papers//paper-229.pdf  "
# Aleksandar Donev a,b, Salvatore Torquato a,b,c,*, Frank H. Stillinger
# Here we improve this method such that we can control the local packing fraction in the polygon.
# -*- coding: utf-8 -*-
import numpy as np
import scipy.stats as stats
import scipy.optimize as opt
import pylab as pl
import matplotlib.gridspec as gridspec
import datetime
import time
from InputReader import *

# Local constants
PI=np.pi
SQRT=np.sqrt
INF = float('inf')

# Compute the cross product of the vextor
def Cross_Product(Point_1,Point_2):
    Cross_Product_value = Point_1[0]*Point_2[1]-Point_1[1]*Point_2[0]
    #print Cross_Product_value
    return Cross_Product_value    

# Compute dot product of the vector
def Dot_Product(Vector_1,Vector_2):
    Dot_Product_Value = Vector_1[0]*Vector_2[0]+Vector_1[1]*Vector_2[1] 
    return Dot_Product_Value
 
# sort the vertices based on the collision time  based on the spartial distribution from the fixed wall and flexible walls  
def AssortBdrNds(VertexList):
    AssortedList = []   
    Temp = sorted(VertexList)
    #print Temp    
    Temp_Copy = [tp for tp in Temp]
    #print Temp_Copy        
    AssortedList.append(Temp[0])
    Temp.remove(AssortedList[-1])
    #print Temp 
    
    for i in Temp:        
        Temp_Values = []
        for j in Temp:
            Vector_1 = [i[0]-AssortedList[-1][0],i[1]-AssortedList[-1][1]]
            Vector_2 = [j[0]-AssortedList[-1][0],j[1]-AssortedList[-1][1]]            
            Temp_Values.append(Cross_Product(Vector_1,Vector_2)) 
            #print i,AssortedList[-1],j,Temp_Values[-1]
        if min(Temp_Values) >= 0:
            AssortedList.append(i)
            #print i            
            Temp_Copy.remove(AssortedList[-1])                        
        else:
            #print "Odd Point",i,j,"Please input a convex DOWN hull"
            continue
    #print  AssortedList   
    UpHull = [up for up in Temp_Copy]
    Temp_U = []
    Temp_U.append(Temp_Copy[0])
    
    for i in range(1,len(Temp_Copy)):
        Temp_Values = []
        for j in range(1,len(Temp_Copy)):
            Vector_1 = [Temp_Copy[i][0]-Temp_U[-1][0],Temp_Copy[i][1]-Temp_U[-1][1]]
            Vector_2 = [Temp_Copy[j][0]-Temp_U[-1][0],Temp_Copy[j][1]-Temp_U[-1][1]]
            Temp_Values.append(Cross_Product(Vector_1,Vector_2)) 
            
        if max(Temp_Values) <= 0:
            Temp_U.append(Temp_Copy[i])            
            UpHull.remove(Temp_U[-1])             
        else:
            print ("Odd Point",Temp_Copy[i],Temp_Copy[j],"Please input a convex UP hull")
            continue        
    Temp_U.reverse()
    #print Temp_U
    for i in Temp_U:
        if i not in AssortedList:
            AssortedList.append(i)
        else:
            continue
    #print AssortedList
    return AssortedList    
    
#Test1 = [[-4,2],[-2,-2],[-4,-1],[3,2],[3,-2],[5,1],[0,-2],[-3,2],[-2.3,2],[-3.5,2],[-1,2],[0,-3]]

#Test = AssortBdrNds(Test1)

def Map(Qdg_Nodes,IsoparametricCS_Point):
    Ksi = IsoparametricCS_Point[0]
    Eta = IsoparametricCS_Point[1]
    Astd_Qdg_Nodes = AssortBdrNds(Qdg_Nodes)
    #print Astd_Qdg_Nodes
    Node_0 =  Astd_Qdg_Nodes[0]
    Node_1 =  Astd_Qdg_Nodes[1]
    Node_2 =  Astd_Qdg_Nodes[2]
    Node_3 =  Astd_Qdg_Nodes[3]
    N0 = (1-Ksi)*(1-Eta)/4
    N1 = (1+Ksi)*(1-Eta)/4
    N2 = (1+Ksi)*(1+Eta)/4
    N3 = (1-Ksi)*(1+Eta)/4
    Glb_X = N0*Node_0[0]+N1*Node_1[0]+N2*Node_2[0]+N3*Node_3[0]
    Glb_Y = N0*Node_0[1]+N1*Node_1[1]+N2*Node_2[1]+N3*Node_3[1]
    
    return([Glb_X,Glb_Y])

Test1=[[-4,-2],[3,2],[-6,3],[2,-3]]
Test2=[[3,2],[2,-3],[8,4],[9,-3]]
Test3=[[11,1],[12,-2],[8,4],[9,-3]]
def Plot_Interpolation(VertexLists,N):
    Glb_X=[]
    Glb_Y=[]
    GX=[]
    GY=[]    
    for vtx in VertexLists:
        X = np.random.uniform(-1,1,N)
        Y = np.random.uniform(-1,1,N)
        #pl.figure(figsize=(6,6))
        #pl.xlim(-1,1)
        #pl.ylim(-1,1)
        #pl.scatter(X,Y)
        for i in range (0,N):
            Temp = Map(vtx,[X[i],Y[i]])
            Glb_X.append(Temp[0])
            Glb_Y.append(Temp[1])
        
        for v in vtx:
            GX.append(v[0])
            GY.append(v[1])
            
    W = max(GX)-min(GX)
    H = max(GY)-min(GY)
    c = 0.5
    c_lim = 1.1
    pl.figure(figsize=(W*c,H*c))
    pl.xlim(c_lim*min(GX),c_lim*max(GX))
    pl.ylim(c_lim*min(GY),c_lim*max(GY))
    pl.scatter(Glb_X,Glb_Y)
    pl.scatter(GX,GY,color="red")
    pl.show()
    return
#Plot_Interpolation([Test1,Test2,Test3],N)
#Test_Triangle=[[0,0],[2,2],[1,1],[4,0]]
#Plot_Interpolation([Test_Triangle],N)
def Calculate_Tri_Area(VertexList):
    if VertexList != []: 
        a = np.sqrt((VertexList[1][0]-VertexList[2][0])**2+(VertexList[1][1]-VertexList[2][1])**2)
        b = np.sqrt((VertexList[2][0]-VertexList[0][0])**2+(VertexList[2][1]-VertexList[0][1])**2)
        c = np.sqrt((VertexList[1][0]-VertexList[0][0])**2+(VertexList[1][1]-VertexList[0][1])**2)
        p = (a+b+c)/2
        S = np.sqrt(p*(p-a)*(p-b)*(p-c))
    else:
        S = 0.0
    return S

def Calculate_Polygon_Area(VertexList):
    Astd_Nodes = AssortBdrNds(VertexList)
    Polygon_Area = 0.
    if VertexList != []:
        for i in range(1,len(Astd_Nodes)-1):
            Polygon_Area = Polygon_Area + Calculate_Tri_Area([Astd_Nodes[0],Astd_Nodes[i],Astd_Nodes[i+1]])
    else:
        Polygon_Area = 0.0
    return Polygon_Area
    
#print Calculate_Polygon_Area([[0,0],[0,4],[4,0]])


def Draw_Circle(center,radius,NOP=180):
    theta = np.linspace(0,2*PI,NOP)
    x=radius * np.cos(theta) + center[0]
    y=radius * np.sin(theta) + center[1]
    pl.plot(x,y,linewidth=0.050,color='r')
    pl.scatter(center[0],center[1], 0.08, 'black' )
    return 
      
    
Test = [[-4,2],[-2,-2],[-4,-1],[3,2],[3,-2],[5,1],[0,-2],[-3,2],[0,-3]]
VL = [[-3,-3],[4,-2],[6,2],[4,5],[-6,4]]
sd1 =[[-6,4],[-3,-3],[-1,-1],[0,2]]
sd2 =[[-3,-3],[-1,-1],[4,-2],[1,0]]
sd3 =[[1,0],[4,-2],[6,2],[0,2]]
sd4 =[[6,2],[4,5],[0,2],[-6,4]]
Test5 = [[[-4,2],[-2,-2]],[[-4,-1],[3,2]],[[3,-2],[5,1]],[[0,-2],[-3,2]],[0,-3]]
#Plot_Interpolation([sd1,sd2,sd3,sd4],N)

def Add_Rgd_Edge(Signal,EndPoints):
    Rgd_Edges = []
    # FOR AUTO ADD ENDPOINTS = [[POINT1],[POINT2],[POINT3]...[POINTN]]
    if Signal == 'Automatic' or Signal == 'Auto' or Signal == 'A':
        Astd_EPs = AssortBdrNds(EndPoints)
        for ep in range (0,len(Astd_EPs)-1):
            Rgd_Edges.append([Astd_EPs[ep],Astd_EPs[ep+1]])
        Rgd_Edges.append([Astd_EPs[-1],Astd_EPs[0]])
        #print Rgd_Edges
        
    # FOR MANUALLY ADD ENDPOINTS = [[[POINT1],[POINT2]]...[[POINTN-1],[POINTN]]]
    elif Signal == 'Manual' or Signal == 'M':
        #print "PLEASE MAKE SURE THE EDGES ARE ASSORTED ANTICLOCKWISE"
        #print "The END POINTS YOU HAVE INPUT ARE:"
        #print EndPoints
        for ep in EndPoints:
            Rgd_Edges.append(ep)
        #print Rgd_Edges
    return Rgd_Edges

    #Add_Rgd_Edge('M',[[[6,2],[4,5]],[[-1,-1],[4,-2]]])

def Add_Flb_Edge(Signal,EndPoints):
    Flb_Edges = []
    # FOR AUTO ADD ENDPOINTS = [[POINT1],[POINT2],[POINT3]...[POINTN]]
    if Signal == 'Automatic' or Signal == 'Auto' or Signal == 'A':
        Astd_EPs = AssortBdrNds(EndPoints)
        for ep in range (0,len(Astd_EPs)-1):
            Flb_Edges.append([Astd_EPs[ep],Astd_EPs[ep+1]])
        Flb_Edges.append([Astd_EPs[-1],Astd_EPs[0]])
        #print Flb_Edges
        
    # FOR MANUALLY ADD ENDPOINTS = [[[POINT1],[POINT2]]...[[POINTN-1],[POINTN]]]
    elif Signal == 'Manual' or Signal == 'M':
        #print "PLEASE MAKE SURE THE EDGES ARE ASSORTED ANTICLOCKWISE"
        #print "The END POINTS YOU HAVE INPUT ARE:"
        #print EndPoints
        for ep in EndPoints:
            Flb_Edges.append(ep)
        #print Flb_Edges
    return Flb_Edges

def nDomain(VertexList,Vf,Radius):
    S = Calculate_Polygon_Area(VertexList)
    n = int(round(S*Vf/PI/Radius**2)) # +1 is modified on 11-01-2017
    return n

# The generated seedlist is in format [X-coord,Y-coord]
def Scatter_Subregion(VertexList,Vf,R):
    Astd_VL = AssortBdrNds(VertexList)
    N = nDomain(VertexList,Vf,R)
    print ('The scattered seed number is'),N
    Seed_List = [] 
    for i in range(0,N):
        X = np.random.uniform(-0.90,0.90,1)
        Y = np.random.uniform(-0.90,0.90,1)
        Glb_XY = Map(Astd_VL,[X[0],Y[0]])    
        Seed_List.append(Glb_XY)    
    return Seed_List

# COMBINE SUBREGION IS NECESSARY FOR ALL CASE  
# Seed_Lists = [[Seed_List1],[Seed_List2],...,[Seed_ListN]]
# Seed_List = [[No_Seed0,Xcoord,Ycoord],[No_Seed2,Xcoord,Ycoord],...,[No_SeedN,Xcoord,Ycoord]]  
def Combine_Subregions(Seed_Lists):
    Seed_List = []
    No_Seed = 0
    for SL in Seed_Lists:
        for sl in SL:
            sl.insert(0,No_Seed) #Assign a unique No. to the seed
            Seed_List.append(sl)
            No_Seed = No_Seed+1
    return Seed_List

#print Combine_Subregions([Scatter_Subregion(sd1,0.5,1.2),Scatter_Subregion(sd3,0.5,1.2)])


def Calculate_Vf(Seed_List,OuterVertexList,InnerVertexLists):
    CircleArea = 0.
    Pure_OVL = []
    for ovl in OuterVertexList:
        if ovl != []:
            Pure_OVL.append(ovl)
    for sl in Seed_List:
        CircleArea = CircleArea + PI*sl[5]**2
    OuterArea = Calculate_Polygon_Area(Pure_OVL)
    InnerArea = 0.0
    
    if InnerVertexLists != [[]]:
        for ivl in InnerVertexLists:
            InnerArea = InnerArea + Calculate_Polygon_Area(ivl)
    else:
        InnerArea = 0.0
    Vf = CircleArea/(OuterArea - InnerArea)
    #print 'Calculate_Vf',Vf
    return Vf         


# The fuction below projects a vector[Vx,Vy] to another vector
def Calculate_Vproj(Vx,Vy,Vector):    
    m = Vector
    m_mdl = SQRT (m[0]**2 + m[1]**2)
    u = [m[0]/m_mdl,m[1]/m_mdl]
    v = [-m[1]/m_mdl,m[0]/m_mdl]
    
    Velocity_collision = Vx*u[0]+Vy*u[1] # coord value in local x direction
    Velocity_parallel = Vx*v[0]+Vy*v[1] # coord value in local y direction
    #Velocity_noncollision = SQRT(Vx**2+Vy**2-Vcollision**2)
    #Vector_collision = [Velocity_collision*u[0],Velocity_collision*u[1]]
    #Vector_noncollision = [Vx-Velocity_collision*u[0],Vy-Velocity_collision*u[1]]
    return [Velocity_collision,Velocity_parallel]


def Transform_from_Local_to_Global(localx,localy,SIN,COS):
    GlobalX = COS*localx - SIN*localy
    GlobalY = SIN*localx + COS*localy
    return [GlobalX,GlobalY]
    
def Transform_from_Global_to_Local(GlobalX,GlobalY,SIN,COS):
    localx = COS*GlobalX + SIN*GlobalY
    localy = -SIN*GlobalX + COS*GlobalY
    return [localx,localy]
#print Transform_from_Global_to_Local(2,3,1/SQRT(5),2/SQRT(5))
#print Transform_from_Global_to_Local(3.1304951684997055,1.7888543819998319,-1/SQRT(5),2/SQRT(5))

#dxy=[-9,-3]
#mdl=SQRT(dxy[0]**2+dxy[1]**2)
#sin=dxy[1]/mdl
#cos=dxy[0]/mdl
#l = Calculate_Vproj(-3,0,dxy)
#gx = Transform_from_Local_to_Global(l[0],l[1],sin,cos)[0]
#print l,gx

#After every collision the direction matrix must be updated
# this function will compute the updated direction matrix
def Calculate_vtclOrient(GivenPoint,EndPoint1,EndPoint2):
    m = GivenPoint[0]
    n = GivenPoint[1]
    a = EndPoint1[0]
    b = EndPoint1[1]
    c = EndPoint2[0]
    d = EndPoint2[1]
    t = c-a
    s = d-b
    e = (a*s**2+m*t**2-b*s*t+n*s*t)/(s**2+t**2)
    f = (b*t**2+n*s**2-a*s*t+m*s*t)/(s**2+t**2)
    COS = (e-m)/SQRT((e-m)**2+(f-n)**2)
    SIN = (f-n)/SQRT((e-m)**2+(f-n)**2)
    module = SQRT((e-m)**2+(f-n)**2)
    M1 = (c-a)*(m-a) + (d-b)*(n-b)
    M2 = (a-c)*(m-c) + (b-d)*(n-d)
    if min(M1,M2) <= 0:
        clsPossibility = False
    else:
        clsPossibility = True
    return [[e-m,f-n],[SIN,COS],module,clsPossibility]
    
#print Calculate_vtclOrient([3,3],[0,-0],[0,-2])
#print Calculate_vtclOrient([3,3],[0,-2],[0,-0])


#Calculation of the collision time in generlized co-ordiantes
def Calculate_gnrlzdTime(cTime,Gamma,Seed_List,AllEdges):
    Edges = []
    for edge in AllEdges:
        if edge != []:
            Edges.append(edge)
        else:
            continue
    
    gnrlzdTime = []
    No_Seed = int(Seed_List[0])
    XY = [Seed_List[1],Seed_List[2]]             
    TempList = []    
    minTime = 0.0    
            
    while AllEdges != [] and minTime <= 0:
        Vx = Seed_List[3]
        Vy = Seed_List[4]
        Radius = Seed_List[5]
        Alpha = Seed_List[6]#Growth_Rate
        R = Gamma*Radius
        
        for edge in Edges:
            TotalList = Calculate_vtclOrient(XY,edge[0],edge[1])
            Vector_Perpendicular = TotalList[0]
            Vector_Orientation  = TotalList[1]
            SIN = Vector_Orientation[0]
            COS = Vector_Orientation[1]
            Vector_Module = TotalList[2]
            clsPossibility = TotalList[3]
            #calculate v_vtcl and v_prll
            VelocityList = Calculate_Vproj(Vx,Vy,Vector_Perpendicular)
            Velocity_collision = VelocityList[0]
            Velocity_parallel = VelocityList[1]            
            
            #CASE 1
            if Velocity_collision <= 0 and Vector_Module >= R and abs(Velocity_collision) > Alpha and clsPossibility == True:  
                tempTime = INF
                Numberlist = [No_Seed]
                NLV = Transform_from_Local_to_Global(Velocity_collision,Velocity_parallel,SIN,COS)
                New_Vx = NLV[0]
                New_Vy = NLV[1]
                NewVelocity = [New_Vx,New_Vy]
                #print 'CASE1',No_Seed,XY,Vector_Module,R,edge,Velocity_collision,Alpha
            
            #CASE 2
            elif Velocity_collision <= 0 and Vector_Module >= R and abs(Velocity_collision) < Alpha and clsPossibility == True:
                tempTime = (Vector_Module - R)/(Velocity_collision + Alpha)
                Numberlist = [No_Seed]
                NLV = Transform_from_Local_to_Global(Velocity_collision - Alpha,Velocity_parallel,SIN,COS)
                New_Vx = NLV[0]
                New_Vy = NLV[1]
                NewVelocity = [New_Vx,New_Vy]
                #print 'CASE2',No_Seed,XY,Vector_Module,R,edge,Velocity_collision,Alpha
                
            #CASE 3
            elif Velocity_collision > 0 and Vector_Module >= R and abs(Velocity_collision) > Alpha and clsPossibility == True:
                tempTime = (Vector_Module - R)/(Velocity_collision + Alpha)
                Numberlist = [No_Seed]
                NLV = Transform_from_Local_to_Global(-Velocity_collision - Alpha,Velocity_parallel,SIN,COS)
                New_Vx = NLV[0]
                New_Vy = NLV[1]
                NewVelocity = [New_Vx,New_Vy]
                #print 'CASE3',No_Seed,XY,Vector_Module,R,edge,Velocity_collision,Alpha
                
            #CASE 4
            elif Velocity_collision > 0 and Vector_Module >= R and abs(Velocity_collision) < Alpha and clsPossibility == True:
                tempTime = (Vector_Module - R)/(Velocity_collision + Alpha)
                Numberlist = [No_Seed]
                NLV = Transform_from_Local_to_Global(-Velocity_collision - Alpha,Velocity_parallel,SIN,COS)
                New_Vx = NLV[0]
                New_Vy = NLV[1]
                NewVelocity = [New_Vx,New_Vy]
                #print 'CASE4',No_Seed,XY,Vector_Module,R,edge,Velocity_collision,Alpha
                
            #CASE 5
            elif Velocity_collision <= 0 and Vector_Module < R and abs(Velocity_collision) > Alpha and clsPossibility == True:
                tempTime = (-Vector_Module + R)/(-Velocity_collision - Alpha)/cTime/0.99
                Numberlist = [No_Seed]
                NLV = Transform_from_Local_to_Global(Velocity_collision - Alpha,Velocity_parallel,SIN,COS)
                New_Vx = NLV[0]
                New_Vy = NLV[1]
                NewVelocity = [New_Vx,New_Vy]
                print ('CASE5'),No_Seed,XY,Vector_Module,R,edge,Velocity_collision,Alpha
                
            #CASE 6
            elif Velocity_collision <= 0 and Vector_Module < R and abs(Velocity_collision) < Alpha and clsPossibility == True:
                tempTime = -INF
                Numberlist = [No_Seed]
                NLV = Transform_from_Local_to_Global(Velocity_collision - Alpha,Velocity_parallel,SIN,COS)
                New_Vx = NLV[0]
                New_Vy = NLV[1]
                Seed_List[3] = New_Vx
                Seed_List[4] = New_Vy
                NewVelocity = [New_Vx,New_Vy]
                #print 'CASE6',No_Seed,XY,Vector_Module,R,edge,Velocity_collision,Alpha
            
            #CASE 7
            elif Velocity_collision > 0 and Vector_Module < R and abs(Velocity_collision) > Alpha and clsPossibility == True:
                tempTime = -INF
                Numberlist = [No_Seed]
                NLV = Transform_from_Local_to_Global(-Velocity_collision - Alpha,Velocity_parallel,SIN,COS)
                New_Vx = NLV[0]
                New_Vy = NLV[1]
                Seed_List[3] = New_Vx
                Seed_List[4] = New_Vy
                NewVelocity = [New_Vx,New_Vy]
                #print 'CASE7',No_Seed,XY,Vector_Module,R,edge,Velocity_collision,Alpha
    
            #CASE 8
            elif Velocity_collision > 0 and Vector_Module < R and abs(Velocity_collision) < Alpha and clsPossibility == True:
                tempTime = -INF
                Numberlist = [No_Seed]
                NLV = Transform_from_Local_to_Global(-Velocity_collision - Alpha,Velocity_parallel,SIN,COS)
                New_Vx = NLV[0]
                New_Vy = NLV[1]
                Seed_List[3] = New_Vx
                Seed_List[4] = New_Vy
                NewVelocity = [New_Vx,New_Vy]
                #print 'CASE8',No_Seed,XY,Vector_Module,R,edge,Velocity_collision,Alpha
            
            #Case 9
            elif Alpha != abs(Velocity_collision) and clsPossibility == False:
                tempTime = INF
                Numberlist = [No_Seed]
                NLV = Transform_from_Local_to_Global(Velocity_collision,Velocity_parallel,SIN,COS)
                New_Vx = NLV[0]
                New_Vy = NLV[1]
                NewVelocity = [New_Vx,New_Vy]
                #print 'CASE9',No_Seed,XY,Vector_Module,R,edge,Velocity_collision,Alpha
    
            #CASE10 Alpha = abs(Velocity_collision)
            else:
                tempTime = -INF
                Seed_List[6] = 0.9*Alpha
                Numberlist = [No_Seed]                
                #print 'CASE10',No_Seed,XY,Vector_Module,R,edge,Velocity_collision,Alpha
    
            
            if tempTime >= 0 and gnrlzdTime == []:
                gnrlzdTime.append(tempTime)
                TempList = [tempTime,Numberlist,NewVelocity,Seed_List[3],Seed_List[4],Seed_List[6],[Vector_Module,R],edge]
                minTime = tempTime
            elif tempTime >= 0 and gnrlzdTime != [] and tempTime < gnrlzdTime[0]:
                gnrlzdTime[0] = tempTime
                TempList = [tempTime,Numberlist,NewVelocity,Seed_List[3],Seed_List[4],Seed_List[6],[Vector_Module,R],edge]
                minTime = tempTime
            else:
                #print 'Error is here'
                continue 

    if AllEdges == []:
        NewVelocity = [Seed_List[3],Seed_List[4]]
        TempList = [INF,No_Seed,NewVelocity,Seed_List[3],Seed_List[4],Seed_List[6],[SQRT(Vx*Vx+Vy*Vy),R]]
            
        #print minTime
        #N = N+1
        #print N
    return TempList #[gnrlzdTime,Numberlist,NewVelocity]
        

#test_seedlist = [1,4,3,0,-0.03,2.8,0.01]
#test_rgdedge = [[[0,0],[-8,0]]]
#print Calculate_rgdTime(test_seedlist,test_rgdedge)
#==============================================================================
#  ⇣⇣⇣ SeedList format [Seed.NO,X-coord,Y-coord,Vx,Vy,Radius,Growth_rate(Alpha)]
#  ⇣⇣⇣ Rgd_Edges = [[Edge1],[Edge2],...,[EdgeN]]
#==============================================================================
def Calculate_rgdTime(cTime,rgdGamma,Seed_List,Rgd_Edges):
    rgdList = Calculate_gnrlzdTime(cTime,rgdGamma,Seed_List,Rgd_Edges)        
    return rgdList #[rgdTime,Numberlist,NewVelocity]

#==============================================================================
#  ⇣⇣⇣ SeedList format [Seed.NO,X-coord,Y-coord,Vx,Vy,Radius,Growth_rate(Alpha)]
#  ⇣⇣⇣ Flb_Edges = [[Edge1],[Edge2],...,[EdgeN]]    
#==============================================================================
def Calculate_flbTime(cTime,flbGamma,Seed_List,Flb_Edges):
    flbList = Calculate_gnrlzdTime(cTime,flbGamma,Seed_List,Flb_Edges)            
    return flbList #[flbTime,Numberlist,NewVelocity]     
             
    

#==============================================================================
# ⇣⇣⇣ SeedLists format [..,[Seed.NO,X-coord,Y-coord,Vx,Vy,Radius,Growth_rate],..]    
#==============================================================================
def Calculate_clsTime(clsGamma,Seed_List):
    clsTime = 0.
    clsTimeList = []
    TempList1 = []
    TempList2 = []
    #Beta = 1.05 # Scale factor of center distance
    for i in range (0,len(Seed_List)):
        SLi = Seed_List[i]
        No_Seedi = int(SLi[0])
        XYi = [SLi[1],SLi[2]]
        Vxi = SLi[3]
        Vyi = SLi[4]
        Ri = SLi[5]
        Alphai = SLi[6] #Growth_Rate
                
            
        for j in range (i+1,len(Seed_List)):    
            SLj = Seed_List[j]           
            No_Seedj = int(SLj[0])
            XYj = [SLj[1],SLj[2]]
            Vxj = SLj[3]
            Vyj = SLj[4]
            Rj = SLj[5]
            Alphaj = SLj[6] #Growth_Rate
    
            dXY = [XYj[0]-XYi[0],XYj[1]-XYi[1]]
            dXY_Module = SQRT(dXY[0]**2 + dXY[1]**2)            

    
            dV = [Vxi-Vxj,Vyi-Vyj]
            dV_Module = SQRT(dV[0]**2 + dV[1]**2)
            SigmaR = (Ri + Rj)*clsGamma
            SigmaAlpha = Alphai + Alphaj
    
            A = dV_Module**2 - SigmaAlpha**2
            B = (SLi[1]-SLj[1])*dV[0] + (SLi[2]-SLj[2])*dV[1] - SigmaR*SigmaAlpha
            C = dXY_Module**2 - SigmaR**2
            
            if B**2 - A*C >= 0.:
                T1 = (-B + SQRT(B**2-A*C))/A
                T2 = (-B - SQRT(B**2-A*C))/A
                if T2 > 0.:
                    Temptime = T2
                elif T1 > 0. and T2 < 0.:
                    Temptime = T1
                else:
                    Temptime = INF
            else:
                Temptime = INF
            
            if clsTimeList != [] and Temptime < clsTime and Temptime > 1e-6:
                clsTime = Temptime
                Numberlist = [No_Seedi,No_Seedj]
            elif clsTimeList == [] and Temptime > 1e-6:
                clsTime = Temptime
                clsTimeList.append(clsTime)
                Numberlist = [No_Seedi,No_Seedj]
            elif clsTimeList == [] and Temptime <= 1e-6:
                clsTime = 0.0
                clsTimeList.append(clsTime)
                TempList1.append([No_Seedi,No_Seedj])
                Numberlist = TempList1
            elif clsTimeList != [] and Temptime < clsTime and Temptime <= 1e-6:
               #clsTimeList != [] and Temptime <= 1e-5:
                clsTime = 0.0
                clsTimeList = [0.0]
                TempList2.append([No_Seedi,No_Seedj])
                Numberlist = TempList2
        
    #clsTime = (dXY_Module - SigmaR)/(SigmaAlpha - dV)
            
    # print [New_Vxi,New_Vyi,New_Vxj,New_Vyj]  
        
    #Numberlist = [No_Seedi,No_Seedj]
    #NewVelocity = [[New_Vxi,New_Vyi],[New_Vxj,New_Vyj]]
            
    return [clsTime,Numberlist]

def Scale_Vf(Seed_List,Vf_target,Vf_present):
    if Vf_present >= Vf_target:
        scale_coefficient = SQRT(Vf_target/Vf_present)
    else:
        scale_coefficient = 1.0
        
    R_Hist = []
    R_oldHist = []
    for sl in Seed_List:
        Radius_old = sl[5]
        R_oldHist.append(Radius_old)        
        Radius_new = scale_coefficient * Radius_old
        sl[5] = Radius_new
        R_Hist.append(Radius_new)
    
    return [Seed_List,R_Hist]

def Plot_Seed(SeedList,OuterVertexList,InnerVertexList,Rgd_Edges,Flb_Edges):
    #X_coord=[]
    #Y_coord=[]

    Outer_X=[]
    Outer_Y=[]    
    Inner_X=[]
    Inner_Y=[]
    green_1 = (210/256.,240/256.,28/256.)
    black = (0.,0.,0.)
    white = (1.,1.,1.)
    for out in OuterVertexList:
        Outer_X.append(out[0])
        Outer_Y.append(out[1])    
    
    maxW = max(Outer_X)-min(Outer_X)
    maxH = max(Outer_Y)-min(Outer_Y)
    Scalefactor = maxW/.5
    W = maxW/Scalefactor
    H = maxH/Scalefactor
    
    for ox in range (0,len(Outer_X)):
        Outer_X[ox] = Outer_X[ox]/Scalefactor
    for oy in range (0,len(Outer_Y)):
        Outer_Y[oy] = Outer_Y[oy]/Scalefactor
    for ix in range (0,len(Inner_X)):
        Inner_X[ix] = Inner_X[ix]/Scalefactor
    for iy in range (0,len(Inner_Y)):
        Inner_Y[iy] = Inner_Y[iy]/Scalefactor
    
    c = 1.576
    c_limax = 1.01
    c_limin = 1.0
    
    fig = pl.figure(figsize=(W*1.573,H*c),dpi=1280)
    fig.patch.set_facecolor('none')
    fig.patch.set_alpha(0.0)
    
    ax = pl.axes()
    #ax.patch.set_alpha(0.0)
    ax.patch.set_facecolor('none')
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
 
    pl.tick_params(axis="both", which="both", bottom="off", top="off",\
                labelbottom="on", left="off", right="off", labelleft="on")
    
    pl.xlim(c_limin*min(Outer_X)-0.02*W,c_limax*max(Outer_X))
    pl.ylim(c_limin*min(Outer_Y)-0.02*H,c_limax*max(Outer_Y))
    
    pl.tick_params(axis="both", which="both", bottom="off", top="off",\
                labelbottom="off", left="off", right="off", labelleft="off")
    for sl in SeedList:
        #X_coord.append(sl[1])
        #Y_coord.append(sl[2])
        #m=SQRT(sl[3]**2 + sl[4]**2)/Scalefactor
        
        sl_x = sl[1]/Scalefactor
        sl_y = sl[2]/Scalefactor
        sl_vx = sl[3]/Scalefactor
        sl_vy = sl[4]/Scalefactor
        sl_R = sl[5]/Scalefactor 
        m = SQRT(sl_vx**2+sl_vy**2)
        
        Draw_Circle([sl_x,sl_y],sl_R,)
        
        headwidth = 0.005
        headlength = 0.01
        ax.arrow(sl_x,sl_y,0.025*sl_vx/m-headwidth,0.025*sl_vy/m-headlength, head_width=headwidth, linewidth=0.02, head_length=headlength, fc='k', ec='k')
        # P.arrow( x, y, dx, dy, **kwargs )
        pl.text(sl_x + 0.01*abs(sl_x),sl_y + 0.01*abs(sl_y),sl[0],color='red',fontsize=0.4)
    
    if len(InnerVertexList) > 1:    
        for inn in InnerVertexList:
            Inner_X.append(inn[0]/Scalefactor)
            Inner_Y.append(inn[1]/Scalefactor) 
            
    for rgd in Rgd_Edges:
        if len(rgd) > 0:
            XList = [rgd[0][0]/Scalefactor,rgd[1][0]/Scalefactor]
            YList = [rgd[0][1]/Scalefactor,rgd[1][1]/Scalefactor]
            pl.plot(XList,YList,'b-',linewidth=0.25)
        
    for flb in Flb_Edges:
        if len(flb) > 0:
            XList = [flb[0][0]/Scalefactor,flb[1][0]/Scalefactor]
            YList = [flb[0][1]/Scalefactor,flb[1][1]/Scalefactor]
            pl.plot(XList,YList,'r--',linewidth=0.2)
    #pl.scatter(X_coord, Y_coord)
    pl.scatter(Outer_X,Outer_Y, 0.5, color="red")
    pl.scatter(Inner_X,Inner_Y, color="green")    
    pl.show()
    #date_time = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
    #pl.savefig("C:\Users\yang\Pictures\python.plot\%s.png"%date_time,format='png', bbox_inches='tight', transparent=True, dpi=1280)

#==============================================================================
#         ⇣⇣⇣ SeedList format [Seed.NO,X-coord,Y-coord,Vx,Vy,Radius,Growth_rate]    
#==============================================================================
# Here starts the main program
#===============================================================================
#                  PART 1 Generation of original seed list.
#===============================================================================
# Step.1 Initialize the required vertex lists (Outer/Inner),Rgd/Flb Boundaries,



#InputFileName = 'Test'
#InputFileName = 'Irregular_three_zones_vf0.6'
#InputFileName = 'Void_three_zones_vf0.65'
#InputFileName = 'Void_three_zones_vf0.58_New'
#InputFileName = 'Void_three_zones_vf0.58'
#InputFileName = 'WeftL1_Six_zones_vf0.76'
#InputFileName = 'WeftL2_Four_zones_vf0.659'
#InputFileName = 'Square_0.5'
#InputFileName = 'Square_0.6'
#InputFileName = 'Square_0.7'
#InputFileName = 'Square_0.8'
#InputFileName = 'Rectangular_0.68'
#InputFileName = 'Lenticular_0.68'
#InputFileName = 'Real_Weft_L1_0.68'
#InputFileName = 'Real_Weft_L2_0.68'
#InputFileName = 'Test_leakage'
#InputFileName = 'Time_Square_0.45'
#InputFileName = 'Time_Square_0.55'
#InputFileName = 'Time_Square_0.65'
#InputFileName = 'X30_0.65'
InputFileName = 'X20-U'
#InputFileName = 'X20-G'



SumList = ReadFile(InputFileName)
#OuterVertexList = [[1,1],[8.5,0.9],[13.2,1.5],[13.0,4.0],[7.0,4.5],[0,4.0],[0.0,2.0]]
OuterVertexList = SumList[0]
InnerVertexList = SumList[1]

#Subregion_1 =[[1.0,1.0],[7,4.5],[0.,4.0],[0,2.0]]
#Subregion_2 =[[1.0,1.0],[8.5,0.9],[13.0,4.0],[7.0,4.5]]
#Subregion_3 =[[8.5,0.9],[13.2,1.5],[13.0,4.0],[13.1,2.75]]

RigidEdgeList = SumList[3]
Rgd_Edges = []
for i in range (0,len(RigidEdgeList)/2):
    if RigidEdgeList[2*i] == 'A':
        if 'void' not in RigidEdgeList[2*i+1] and 'Void' not in RigidEdgeList[2*i+1]:
            tempEdges = []
            tempEdges =  Add_Rgd_Edge('A',OuterVertexList)
            for te in tempEdges:
                Rgd_Edges.append(te)            
        else:
            Rgd_Edges.append([])
    elif RigidEdgeList[2*i] == 'M':
        if 'void' not in RigidEdgeList[2*i+1] and 'Void' not in RigidEdgeList[2*i+1]:
            tempEdges = []
            tempEdges = Add_Rgd_Edge('M',[RigidEdgeList[2*i+1]])
            for te in tempEdges:
                Rgd_Edges.append(te)
        else:
            Rgd_Edges.append([])
Temp_Rgd = []
for rgd in Rgd_Edges:
    if rgd != []:
        Temp_Rgd.append(rgd)
Rgd_Edges = Temp_Rgd


FlexibleEdgeList = SumList[4]
#print SumList[4]
Flb_Edges = []
for i in range (0,len(FlexibleEdgeList)/2):
    if FlexibleEdgeList[2*i] == 'A':
        if 'void' not in FlexibleEdgeList[2*i+1] and 'Void' not in FlexibleEdgeList[2*i+1]:
            tempEdges = []
            tempEdges = Add_Rgd_Edge('A',OuterVertexList)
            for te in tempEdges:
                Flb_Edges.append(te)            
        else:
            Flb_Edges.append([])
    elif FlexibleEdgeList[2*i] == 'M':
        if 'void' not in FlexibleEdgeList[2*i+1] and 'Void' not in FlexibleEdgeList[2*i+1]:
            tempEdges = []        
            tempEdges = Add_Flb_Edge('M',FlexibleEdgeList[2*i+1])
            #print FlexibleEdgeList[2*i+1]
            for te in tempEdges:
                Flb_Edges.append(te)
        else:
            Flb_Edges.append([])
#print '****',Flb_Edges
#Flb_Edges = Add_Flb_Edge('M',[[[1.0,1.0],[7,4.5]],[[8.5,0.9],[13,4.0]]])
R_tgt = SumList[5]
DistributionType = SumList[6]
std_dev = SumList[7]


# Step.2 Scatter subregions and Combine the generated subregions

#SeedList = Combine_Subregions([Scatter_Subregion(Subregion_1,0.5,R_tgt),Scatter_Subregion(Subregion_2,0.6,R_tgt)\
#                              ,Scatter_Subregion(Subregion_3,0.7,R_tgt)])

SubregionsList = SumList[2]
Seed_Lists = []
for i in range (0,len(SubregionsList)/2):
    Seed_Lists.append(Scatter_Subregion(SubregionsList[2*i],SubregionsList[2*i+1],R_tgt))
SeedList = Combine_Subregions(Seed_Lists)
#ST_SL = Combine_Subregions([Scatter_Subregion(ST_OVL,0.8,R_tgt)])

#Plot_Seed(SeedList,OuterVertexList,InnerVertexList)

# Step.3 Assign random Velocity and Growth rate
for sl in SeedList:
    Vx = np.random.uniform(-R_tgt/2.,R_tgt/2.,1)
    Vy = np.random.uniform(-R_tgt/2.,R_tgt/2.,1)
    R0 = 0.
    if DistributionType == 'Gaussian' or DistributionType == 'G' or DistributionType == 'g':
        grth = np.random.normal(1.0*R_tgt,std_dev,1)
    elif DistributionType == 'Uniform' or DistributionType == 'U' or DistributionType == 'u':
        grth = np.random.uniform((1.0-std_dev)*R_tgt/5.,(1.0+std_dev)*R_tgt/5.,1)
    elif DistributionType == 'Constant' or DistributionType == 'C' or DistributionType == 'c':
        grth = [R_tgt/5.0]
    else:
        grth = [1.0]
        
    for data in [Vx[0],Vy[0],R0,grth[0]]:
        sl.append(data)

Plot_Seed(SeedList,OuterVertexList,InnerVertexList,Rgd_Edges,Flb_Edges)
#Plot_Seed(ST_SL,ST_OVL,ST_IVL,ST_RgdEdges,ST_FlbEdges)

#===============================================================================
#PART 2 Growth and motions inside the domain until the specified Vf is satisfied
#===============================================================================
# Step.4 Calculate dTime
Ninc = 1
vf = 0.0
Vf = SumList[8]
Vf_target = Vf
rgdGamma = SumList[9]
flbGamma = SumList[10]
clsGamma = SumList[11]
cTime = SumList[12]
N_total = SumList[13]
N_interval = SumList[14]
dTime = 0.0
Time = 0.0
dVfdTime = []
date_time = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
#f = open("D:\Caba/Random_SeedList_%s.dat"% (date_time),"w+")
f = open("Random_SeedList_%s.dat"% (date_time),"w+")

T_Origin = time.time()
while Ninc <= N_total and vf < Vf:    
    
# Step.5 Move all the circles to Time + dTime (Update position, velocity, Radius)
    T_Start = time.time()
    
    rgdTimeList = []
    flbTimeList = []
    clsTimeList = []
    rgdTime = 0.
    flbTime = 0.
    UpdatedSeed = []
    
    clsTimeList = Calculate_clsTime(clsGamma,SeedList)
    clsTime = clsTimeList[0]
    for i in range (0,len(SeedList)):
        rgdTempList = Calculate_rgdTime(cTime,rgdGamma,SeedList[i],Rgd_Edges)
        rgdTemptime = rgdTempList[0]
        No = rgdTempList[1][0]
        SeedList[No][3] = rgdTempList[3]
        SeedList[No][4] = rgdTempList[4]
        SeedList[No][6] = rgdTempList[5]        
        if rgdTimeList != [] and rgdTemptime < rgdTime:
            rgdTime = rgdTemptime
            rgdTimeList = rgdTempList
        elif rgdTimeList == []:
            rgdTime = rgdTemptime
            rgdTimeList = rgdTempList
            
        flbTempList = Calculate_flbTime(cTime,flbGamma,SeedList[i],Flb_Edges)
        flbTemptime = flbTempList[0]
        No = rgdTempList[1][0]
        SeedList[No][3] = rgdTempList[3]
        SeedList[No][4] = rgdTempList[4]
        SeedList[No][6] = rgdTempList[5]
        if flbTimeList != [] and flbTemptime < flbTime:
            flbTime = flbTemptime
            flbTimeList = flbTempList
        elif flbTimeList == []:
            flbTime = flbTemptime
            flbTimeList = flbTempList
        
    dTime = min(clsTime,rgdTime,flbTime)
    #Plot_Seed(SeedList,OuterVertexList,InnerVertexList,Rgd_Edges,Flb_Edges)    
    for sl in SeedList:
        No_Seed = sl[0]
        X_coord = sl[1]
        Y_coord = sl[2]
        Vx = sl[3]
        Vy = sl[4]
        Radius = sl[5]
        Alpha = sl[6]
         
        New_X = X_coord + Vx * dTime
        New_Y = Y_coord + Vy * dTime
        
        New_Radius = Radius + cTime*Alpha*dTime
        sl[1] = New_X
        sl[2] = New_Y
        sl[5] = New_Radius
        
    
    if dTime == clsTime and clsTime > 0.:
        i = clsTimeList[1][0]
        j = clsTimeList[1][1]
        #print i,j,'collision'
        XYi = [SeedList[i][1],SeedList[i][2]]
        XYj = [SeedList[j][1],SeedList[j][2]]  
        Ri = SeedList[i][5]
        Rj = SeedList[j][5]
        dXY = [XYj[0]-XYi[0],XYj[1]-XYi[1]]
        dXY_Module = SQRT(dXY[0]**2 + dXY[1]**2)            

        COS = dXY[0]/dXY_Module
        SIN = dXY[1]/dXY_Module
            
        Vxi = SeedList[i][3]
        Vyi = SeedList[i][4]        
        Vxj = SeedList[j][3]
        Vyj = SeedList[j][4]
        
        Alphai = SeedList[i][6]
        Alphaj = SeedList[j][6]
        SigmaAlpha = Alphai + Alphaj        
            
        VHi = Calculate_Vproj(Vxi,Vyi,dXY)[0]
        VLi = Calculate_Vproj(Vxi,Vyi,dXY)[1]
        VHj = Calculate_Vproj(Vxj,Vyj,dXY)[0]
        VLj = Calculate_Vproj(Vxj,Vyj,dXY)[1]
                       
        UpdatedSeed = [i,j,[VHi,VHj],dXY_Module,Ri+Rj]
            
        New_VHi_local = -abs(VHj) - 1.0*SigmaAlpha
        #New_VLi_local = VLi
        New_VLi_local = np.random.uniform(-R_tgt/2.,R_tgt/2.,1)[0]
        New_VHj_local = abs(VHi) + 1.0*SigmaAlpha
        #New_VLj_local = VLj*VLi/abs(VLi)*2.0
        New_VLj_local = np.random.uniform(-R_tgt/2.,R_tgt/2.,1)[0]
            
        New_Vxi = Transform_from_Local_to_Global(New_VHi_local,New_VLi_local,SIN,COS)[0]
        New_Vyi = Transform_from_Local_to_Global(New_VHi_local,New_VLi_local,SIN,COS)[1]
        New_Vxj = Transform_from_Local_to_Global(New_VHj_local,New_VLj_local,SIN,COS)[0]
        New_Vyj = Transform_from_Local_to_Global(New_VHj_local,New_VLj_local,SIN,COS)[1]
            
        #New_Vxi = np.random.uniform(-R_tgt/2.,R_tgt/2.,1)[0]
        #New_Vyi = np.random.uniform(-R_tgt/2.,R_tgt/2.,1)[0]
        #New_Vxj = np.random.uniform(-R_tgt/2.,R_tgt/2.,1)[0]
        #New_Vyj = np.random.uniform(-R_tgt/2.,R_tgt/2.,1)[0]
            
        SeedList[i][3] = New_Vxi
        SeedList[i][4] = New_Vyi
        SeedList[j][3] = New_Vxj
        SeedList[j][4] = New_Vyj
    
    elif dTime == clsTime and clsTime == 0.:
        NumberList = clsTimeList[1]
        
        for sl in NumberList:
            i = sl[0]
            j = sl[1]
            
            XYi = [SeedList[i][1],SeedList[i][2]]
            XYj = [SeedList[j][1],SeedList[j][2]]  

            dXY = [XYj[0]-XYi[0],XYj[1]-XYi[1]]
            dXY_Module = SQRT(dXY[0]**2 + dXY[1]**2)            
    
            COS = dXY[0]/dXY_Module
            SIN = dXY[1]/dXY_Module
                
            New_VHi_local = np.random.uniform(-R_tgt/2.,R_tgt/2.,1)[0]
            New_VLi_local = np.random.uniform(-R_tgt/2.,R_tgt/2.,1)[0]
            New_VHj_local = np.random.uniform(-R_tgt/2.,R_tgt/2.,1)[0]
            New_VLj_local = np.random.uniform(-R_tgt/2.,R_tgt/2.,1)[0]
                
            New_Vxi = Transform_from_Local_to_Global(New_VHi_local,New_VLi_local,SIN,COS)[0]
            New_Vyi = Transform_from_Local_to_Global(New_VHi_local,New_VLi_local,SIN,COS)[1]
            New_Vxj = Transform_from_Local_to_Global(New_VHj_local,New_VLj_local,SIN,COS)[0]
            New_Vyj = Transform_from_Local_to_Global(New_VHj_local,New_VLj_local,SIN,COS)[1]
            
            SeedList[i][3] = New_Vxi
            SeedList[i][4] = New_Vyi
            SeedList[j][3] = New_Vxj
            SeedList[j][4] = New_Vyj
        
    elif dTime == rgdTime:
        i = rgdTimeList[1][0]
        edge = rgdTimeList[-1]
        UpdatedSeed = [i,edge]
        #print i,'rgdedge',rgdTimeList[-1]
        New_Vx = rgdTimeList[2][0]
        New_Vy = rgdTimeList[2][1]
        SeedList[i][3] = New_Vx
        SeedList[i][4] = New_Vy
    
    elif dTime == flbTime:
        i = flbTimeList[1][0]
        edge = flbTimeList[-1]
        UpdatedSeed = [i,edge]
        
        #print i,'flbedge'        
        New_Vx = flbTimeList[2][0]
        New_Vy = flbTimeList[2][1]
        SeedList[i][3] = New_Vx
        SeedList[i][4] = New_Vy
        

# Step.6 Caculate current Vf and update the Ninc
    vf = Calculate_Vf(SeedList,OuterVertexList,InnerVertexList)
    #print 'vf is',vf
    Ninc = Ninc + 1
    Time = Time + dTime
    T_End = time.time()
    deltaT = T_End - T_Start
    totalT = T_End - T_Origin
    dVfdTime.append([Ninc,vf,Time,dTime,deltaT,totalT])
    for n in range (0,int(N_total/N_interval)):
        if Ninc in [n*N_interval + 1]:
            Plot_Seed(SeedList,OuterVertexList,InnerVertexList,Rgd_Edges,Flb_Edges)
            print('parameters'), vf,Ninc,dTime,UpdatedSeed
            f.write('Increment N%i'%Ninc)
            f.write('\n')
            for sl in SeedList:
                f.write(str(sl))
                f.write('\n')

print ('BEFORE SCALE'),vf
Plot_Seed(SeedList,OuterVertexList,InnerVertexList,Rgd_Edges,Flb_Edges)
        
# PART 3 Scale the circles to get a precise Vf value


FinalList = Scale_Vf(SeedList,Vf,vf)
SeedList = FinalList[0]
R_Hist =  FinalList[1]
# PART 4 Plot pictures
vf = Calculate_Vf(SeedList,OuterVertexList,InnerVertexList) 
print ('AFTER SCALE'),vf
Plot_Seed(SeedList,OuterVertexList,InnerVertexList,Rgd_Edges,Flb_Edges)
pl.hist(R_Hist,bins = 15)
pl.show()

g = open("dVfdTime_Random_SeedList_%s.dat"% (date_time),"w+")
g.write('# The input file name is %s.idt\n'%InputFileName)
g.write('Ninc\tvf\tTime\tdTime\tdeltaT [s]\ttotalT [s]\n')
g.write('1\t0.0\t0.0\t0.0\t0.0\t0.0\n')
for d in dVfdTime:
    g.write(str(d[0]))
    g.write('\t')
    g.write(str(d[1]))
    g.write('\t')
    g.write(str(d[2]))
    g.write('\t')
    g.write(str(d[3]))
    g.write('\t')
    g.write(str(d[4]))
    g.write('\t')
    g.write(str(d[5]))
    g.write('\n')


Time_Now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
h = open("Final_Random_SeedList_%s.dat"% (date_time),"w+")
h.write('# The input file name is %s.idt\n'%InputFileName)
h.write('# Generated on %s\n'%Time_Now)
h.write('# The total number of seed is %i\n'%len(SeedList))
h.write('# The final fiber volume fraction is %f\n'%vf)
h.write('*SeedList\n')
for sl in SeedList:
    h.write(str(sl))
    h.write('\n')
h.write('*OuterVertexList\n')
for ovl in (AssortBdrNds(OuterVertexList)):
    h.write(str(ovl))
    h.write('\n')
h.write('*InnerVertexList\n')
for ivl in (AssortBdrNds(InnerVertexList)):
    h.write(str(ivl))
    h.write('\n')

#GammaList = [rgdGamma,flbGamma,clsGamma]
#New_State = DriveRgdMothion(Vf,SeedList,1000,50,OuterVertexList,InnerVertexList,Rgd_Edges,Flb_Edges,GammaList,cTime)
#New_SeedList = New_State[0]
#DriveRgdMothion(Vf,New_SeedList,1500,50,OuterVertexList,InnerVertexList,Rgd_Edges,Flb_Edges,GammaList,cTime)
g.close()    
f.close()    
h.close()



