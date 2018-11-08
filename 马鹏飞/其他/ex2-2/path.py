import numpy as np
from scipy.interpolate import Rbf
from scipy.spatial     import Voronoi
from ctypes import *
import math

def initialProbability(points):
    # a list of tuples. 
    # each tuple indicates a point.
    u1 = -2.0; u2 = -1.8;# u1 = -2			
    s1 = s2 = np.sqrt(1);# u2 = -1.8		
    r12 = 0.0;# s1 = s2 = 1.0		# 
    # r12 = 0.0
    value=[]        # np.zeros(len(points))
    for i in range(len(points)):
        x = points[i][0]
        y = points[i][1]
        temp = s2**2*(x - u1)**2 + s1**2*(y - u2)**2 - 2 * r12*s1*s2*(y - u2)*(x - u1)
        temp = np.exp(-0.5*temp / ((1.0 - r12*r12)*s1*s1*s2*s2))
        temp = temp / (2.0*np.pi*s1*s2*np.sqrt(1.0 - r12*r12))
        value.append(temp)
    return value

def initialMoment(points):
    moments=[]
    n = len(points)
    for i in range(n):
        moment=[]
        moment.append(points[i][0])
        moment.append(points[i][1])
        moment.append(points[i][0]*points[i][1])
        moment.append(points[i][0]*points[i][0])
        moment.append(points[i][1]*points[i][1])
        moments.append(moment)
    return moments

def oderk4(moments,t,dt):
    def RHS(M10, M01, M11, M20, M02, TM, MARK):
        GAMA = 0.2; BETA = 0.0; EPSON = 0.3; OMIGA1 = 1.0; OMIGA = 1.6; SIGMA1 = 1.0; SIGMA2 = 0.4
        if MARK == 1:
            ret = M01
        if MARK == 2:
            ret = -(GAMA*BETA / OMIGA1 / OMIGA1)*(3.0*M01*M02 - 2.0*M01*M01*M01) - OMIGA1*OMIGA1*EPSON*(3.0*M10*M20 - 2.0*M10*M10*M10) - GAMA*M01 - OMIGA1*OMIGA1*M10 + OMIGA1*SIGMA1*np.sin(OMIGA*TM)
        if MARK == 3:
            ret = M02 + (GAMA*BETA / OMIGA1 / OMIGA1)*(2.0*M01*M01*M01*M10 - 3.0*M02*M11) + OMIGA1*OMIGA1*EPSON*(2.0*M10*M10*M10*M10 - 3.0*M20*M20) - GAMA*M11 - OMIGA1*OMIGA1*M20 + OMIGA1*SIGMA1*M10*np.sin(OMIGA*TM)
        if MARK == 4:
            ret = 2.0*M11
        if MARK == 5:
            ret = (GAMA*BETA / OMIGA1 / OMIGA1)*(4.0*M01*M01*M01*M01 - 6.0*M02*M02) + OMIGA1*OMIGA1*EPSON*(4.0*M10*M10*M10*M01 - 6.0*M20*M11) - 2.0*GAMA*M02 - 2.0*OMIGA1*OMIGA1*M11 + 2.0*OMIGA1*SIGMA1*M01*np.sin(OMIGA*TM)+OMIGA1*OMIGA1*SIGMA2*SIGMA2
        return ret
    # it start now !
    n=len(moments)
    for i in range(n):
        M10 = moments[i][0]
        M01 = moments[i][1]
        M11 = moments[i][2]
        M20 = moments[i][3]
        M02 = moments[i][4]
        DT0 = dt / 100
        TM0 = t - dt
        while np.abs(TM0 - t) > 1.0e-3:
            M10_1 = M10 + DT0*RHS(M10, M01, M11, M20, M02, TM0, 1)
            M01_1 = M01 + DT0*RHS(M10, M01, M11, M20, M02, TM0, 2)
            M11_1 = M11 + DT0*RHS(M10, M01, M11, M20, M02, TM0, 3)
            M20_1 = M20 + DT0*RHS(M10, M01, M11, M20, M02, TM0, 4)
            M02_1 = M02 + DT0*RHS(M10, M01, M11, M20, M02, TM0, 5)
            M10_2 = 0.75*M10 + 0.25*M10_1 + 0.25*DT0*RHS(M10_1, M01_1, M11_1, M20_1, M02_1, TM0, 1)
            M01_2 = 0.75*M01 + 0.25*M01_1 + 0.25*DT0*RHS(M10_1, M01_1, M11_1, M20_1, M02_1, TM0, 2)
            M11_2 = 0.75*M11 + 0.25*M11_1 + 0.25*DT0*RHS(M10_1, M01_1, M11_1, M20_1, M02_1, TM0, 3)
            M20_2 = 0.75*M20 + 0.25*M20_1 + 0.25*DT0*RHS(M10_1, M01_1, M11_1, M20_1, M02_1, TM0, 4)
            M02_2 = 0.75*M02 + 0.25*M02_1 + 0.25*DT0*RHS(M10_1, M01_1, M11_1, M20_1, M02_1, TM0, 5)
            M10_3 = M10 / 3.0 + 2.0*M10_2 / 3.0 + 2.0*DT0*RHS(M10_2, M01_2, M11_2, M20_2, M02_2, TM0, 1) / 3.0
            M01_3 = M01 / 3.0 + 2.0*M01_2 / 3.0 + 2.0*DT0*RHS(M10_2, M01_2, M11_2, M20_2, M02_2, TM0, 2) / 3.0
            M11_3 = M11 / 3.0 + 2.0*M11_2 / 3.0 + 2.0*DT0*RHS(M10_2, M01_2, M11_2, M20_2, M02_2, TM0, 3) / 3.0
            M20_3 = M20 / 3.0 + 2.0*M20_2 / 3.0 + 2.0*DT0*RHS(M10_2, M01_2, M11_2, M20_2, M02_2, TM0, 4) / 3.0
            M02_3 = M02 / 3.0 + 2.0*M02_2 / 3.0 + 2.0*DT0*RHS(M10_2, M01_2, M11_2, M20_2, M02_2, TM0, 5) / 3.0
            M10 = M10_3
            M01 = M01_3
            M11 = M11_3
            M20 = M20_3
            M02 = M02_3
            TM0 = TM0 + DT0
        moments[i][0] = M10
        moments[i][1] = M01
        moments[i][2] = M11
        moments[i][3] = M20
        moments[i][4] = M02
    return moments

def pathIntegration(points,moments,area,probability_old):
    # the probability function
    def p(point, moment):
        try:
            x = point[0]
            y = point[1]
            u1 = moment[0]; u2 = moment[1]
            s1 = np.sqrt(moment[3] - moment[0]*moment[0])
            s2 = np.sqrt(moment[4] - moment[1]*moment[1])
            r12 = (moment[2] - moment[0]*moment[1]) / (s1*s2)
            temp = s2*s2*(x - u1)*(x - u1) + s1*s1*(y - u2)*(y - u2) - 2.0 * r12*s1*s2*(y - u2)*(x - u1)
            temp = np.exp(-0.5*temp / ((1.0 - r12*r12)*s1*s1*s2*s2))
            ret = temp / (2.0*np.pi*s1*s2*np.sqrt(1.0 - r12*r12))
        except:
            ret = 0.0
        if math.isnan(ret):
            ret = 0.0
        return ret
    # transform lists to arrays
    n=len(points)
    moments = np.array(moments)
    points  = np.array(points)
    area    = np.array(area)
    # now, it start !
    probability_new=[]
    for i in range(n):
        temp=0
        for j in range(n):
            # call the probability function
            if probability_old[j] > 0.000001:
                pp = p(points[i], moments[j])
                temp += area[j]*pp*probability_old[j]
        probability_new.append(temp)
    return probability_new


def changePointsAndProbability(points,probability):
    # unwrap points
    x=[];y=[]
    for point in points:
        x.append(point[0])
        y.append(point[1])
    # get new points
    z=[]
    for p in probability:
        if p<0:
            z.append(0.15)
        else:
            if p>0.16:
                z.append(0.05)
            else:
                z.append(0.15-0.625*p)
    # bubble point
    points_new = bubblePoints(x, y, z, 5.0)

    x_new=[];y_new=[]
    for point in points_new:
        x_new.append(point[0])
        y_new.append(point[1])
    
    rbf_p = Rbf(x,y,probability)
    probability_new = list(rbf_p(x_new,y_new))

    return points_new,probability_new



def bubblePoints(x,y,z,limits):
    if len(x) != len(y) | len(x) != len(y):
        print("x y z should be lists with the same length")
        return 0
    functype_ideal = CFUNCTYPE(c_float,c_float,c_float)
    bubble = cdll.LoadLibrary("bubble.dll").bubblePoint
    bubble.restype = POINTER(c_double)
    bubble.argtypes = (c_double,functype_ideal,POINTER(c_int))
    rbf=Rbf(x,y,z)
    def distance(a,b):
        return rbf([a],[b])
    points=[]
    n = c_int()
    point=bubble(limits,functype_ideal(distance),n)
    for i in range(n.value):
        points.append((point[i],point[i+n.value]))
    return points


def areas(points):
	# calculate the area of a poly
    def poly_area(vertices):
        if not isinstance(vertices,np.ndarray):
            print("the type of input should be array when calculate the poly area")
            return 0
        area = 0
        for i in range(vertices.shape[0]-2):
            area += tri_area(vertices[[-1,i,i+1]])
        return np.abs(area)
    # calculate the area of a triangle 
    def tri_area(vertices):
        if not isinstance(vertices,np.ndarray):
            print("the type of input should be array when calculate the tri area")
            return 0
        if vertices.shape[0] != 3:
            print("this is not a triangle")
            return 0
        x1=vertices[0][0]
        y1=vertices[0][1]
        x2=vertices[1][0]
        y2=vertices[1][1]
        x3=vertices[2][0]
        y3=vertices[2][1]
        return 0.5*(x1*y2 + x2*y3 + x3*y1 - x1*y3 - x2*y1 - x3*y2)
    # calculate voronoi diagram
    vor=Voronoi(points)
    area_list=[]
    for i in range(len(points)):
        region = vor.regions[vor.point_region[i]]
        # the area of polys on the edge are zero.
        area=0
        if -1 not in region:
            area=poly_area(vor.vertices[region])
        area_list.append(area)
    return area_list


# time variables
dt  = 2.0*np.pi/1.6
t   = dt

# initial mesh grid
y   = np.linspace(-1,1,50)*5.0
x   = np.linspace(-1,1,50)*5.0
x,y = np.meshgrid(x,y)
x,y = x.ravel(),y.ravel()

# these are initial variables
points      = list(zip(x,y))
area        = areas(points)
moments     = initialMoment(points)
moments     = oderk4(moments, t, dt)
probability = initialProbability(points)
probability = pathIntegration(points,moments,area,probability)  

# these const variables are used to interpolate new moments
# const_moments = moments
# const_points  = points


# start the loop
while t<160:
	
    t = t + dt

    points, probability = changePointsAndProbability(points, probability)
    
    moments  =  initialMoment(points)
    moments  =  oderk4(moments, dt, dt)
    area     =  areas(points)
    
    probability  =  pathIntegration(points,moments,area,probability)
    
    volume   =  0.0
    
    for i in range(len(probability)):
        volume += area[i]*probability[i]
    
    for i in range(len(probability)):
        probability[i] /= volume
    
    # output the outcomes in every step
    file_name = '%.2f.txt' % t
    f = open(file_name, "w") 
    for i in range(len(points)):
        print("%.4f %.4f %.4f" % (points[i][0], points[i][1],probability[i]), file = f)
    f.close()