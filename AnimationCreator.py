import numpy as np
import scipy.stats as stats
import scipy.optimize as opt
import pylab as pl
#import matplotlib.gridspec as gridspec
import datetime
import time

PI=np.pi
SQRT=np.sqrt
INF = float('inf')
def Cross_Product(Point_1,Point_2):
    Cross_Product_value = Point_1[0]*Point_2[1]-Point_1[1]*Point_2[0]
    #print Cross_Product_value
    return Cross_Product_value    

def Dot_Product(Vector_1,Vector_2):
    Dot_Product_Value = Vector_1[0]*Vector_2[0]+Vector_1[1]*Vector_2[1] 
    return Dot_Product_Value
    
def AssortBdrNds(VertexList):
    AssortedList = []   
    Temp = sorted(VertexList)    
    Temp_Copy = [tp for tp in Temp]        
    AssortedList.append(Temp[0])
    Temp.remove(AssortedList[-1])
    
    for i in Temp:        
        Temp_Values = []
        for j in Temp:
            Vector_1 = [i[0]-AssortedList[-1][0],i[1]-AssortedList[-1][1]]
            Vector_2 = [j[0]-AssortedList[-1][0],j[1]-AssortedList[-1][1]]            
            Temp_Values.append(Cross_Product(Vector_1,Vector_2))            
        if min(Temp_Values) >= 0:
            AssortedList.append(i)            
            Temp_Copy.remove(AssortedList[-1])                        
        else:
            #print "Odd Point",i,j,"Please input a convex DOWN hull"
            continue
        
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
            print ("Odd Point"),Temp_Copy[i],Temp_Copy[j],"Please input a convex UP hull"
            continue        
    Temp_U.reverse()
    for i in Temp_U:
        if i not in AssortedList:
            AssortedList.append(i)
        else:
            continue
    return AssortedList    
    
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
def Draw_Circle(center,radius,NOP=180):
    theta = pl.linspace(0,2*PI,NOP)
    x=radius * np.cos(theta) + center[0]
    y=radius * np.sin(theta) + center[1]
    pl.plot(x,y,linewidth=0.8,color='red')
    pl.scatter(center[0],center[1],linewidth=0.5, )
    return 

def Plot_Seed(SeedList,OuterVertexList,InnerVertexList,Rgd_Edges,Flb_Edges):
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
    Scalefactor = maxW/10.0
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
    
    c = 1.5
    c_limax = 1.01
    c_limin = 1.0
    
    pl.figure(figsize=(W*c,H*c),dpi=1280)
    ax = pl.axes()
    ax.patch.set_facecolor(white)
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
        
        headwidth = 0.05
        headlength = 0.1
        #ax.arrow(sl_x,sl_y,0.25*sl_vx/m-headwidth,0.25*sl_vy/m-headlength, head_width=headwidth, head_length=headlength, fc='k', ec='k')
        # P.arrow( x, y, dx, dy, **kwargs )
        pl.text(sl_x + 0.01*abs(sl_x),sl_y + 0.01*abs(sl_y),sl[0],color='red',fontsize=14)
    
    if len(InnerVertexList) > 1:    
        for inn in InnerVertexList:
            Inner_X.append(inn[0]/Scalefactor)
            Inner_Y.append(inn[1]/Scalefactor) 
            
    for rgd in Rgd_Edges:
        if len(rgd) > 0:
            XList = [rgd[0][0]/Scalefactor,rgd[1][0]/Scalefactor]
            YList = [rgd[0][1]/Scalefactor,rgd[1][1]/Scalefactor]
            pl.plot(XList,YList,'b-',linewidth=2.5)
        
    for flb in Flb_Edges:
        if len(flb) > 0:
            XList = [flb[0][0]/Scalefactor,flb[1][0]/Scalefactor]
            YList = [flb[0][1]/Scalefactor,flb[1][1]/Scalefactor]
            pl.plot(XList,YList,'r--',linewidth=2)
    pl.scatter(X_coord, Y_coord)
    pl.scatter(Outer_X,Outer_Y,color="red")
    pl.scatter(Inner_X,Inner_Y,color="green")    
    pl.show()
    
#OuterVertexList = [[0,163.7],[56.5,183],[132,151],[200,110],[190,30],[156,0],[0,0],[180,129],[0,86]]
#OuterVertexList = [[0,0],[130,0],[0,130],[130,130]]
OuterVertexList = [[-6,-6],[6.0,-6],[6,6.0],[-6.0,6.0]]

InnerVertexList = [[]]
Rgd_Edges = Add_Rgd_Edge('A',OuterVertexList)
#Flb_Edges = Add_Flb_Edge('M',[[[132,151],[100,86]],[[100,86],[112,80]],[[112,80],[180,129]],\
            #[[0,0],[100,86]],[[112,80],[156,0]],[[56.5,183],[65,154]],[[65,154],[95,144]],[[95,144],[132,151]]])
#Flb_Edges = Add_Flb_Edge('M',[[[0,0],[0.1,0.1]]])
Flb_Edges = Add_Flb_Edge('M',[[[-6,5],[6,5]],[[-6,-5],[6,-5]],[[-5,5],[-5,-5]],[[5,5],[5,-5]]])
#f=open('Random_SeedList_2017-04-14_12-08-19.dat','r')
#f=open('Random_SeedList_2017-04-14_10-55-49.dat','r')
#f=open('Random_SeedList_2017-04-12_17-34-15.dat','r')
#f=open('Final_Random_SeedList_2017-04-27_01-03-13.dat','r')
#f=open('Final_Random_SeedList_2017-04-12_21-58-22.dat','r')
f=open('Random_SeedList_2017-08-27_16-35-54.dat','r')

#f=open('.dat','r')
#f=open('SeedList_Test.dat','r')

line = f.readline()


while len(line) != 0:
    #print line[0:-1],len(line)
    seedlist = []
    while line[0:1] == "[":
        seedlist.append(line)
        line = f.readline()
        
    line = f.readline()

    SeedList = []
    
    for tl in seedlist:
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
        templist.append(float(tl[commalist[2]+2:commalist[3]]))
        templist.append(float(tl[commalist[3]+2:commalist[4]]))
        templist.append(float(tl[commalist[4]+2:commalist[5]])) 
        templist.append(float(tl[commalist[5]+2:-2]))
        
        #datalist.append(templist)
        SeedList.append(templist)        
    #print SeedList
    #Plot_Seed(SeedList,OuterVertexList,InnerVertexList,Rgd_Edges,Flb_Edges)
    #print seedlist
Plot_Seed(SeedList,OuterVertexList,InnerVertexList,Rgd_Edges,Flb_Edges)