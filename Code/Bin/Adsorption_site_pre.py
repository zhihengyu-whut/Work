# -*- coding: utf-8 -*-
import math
def PointOnLine(Points):
    '''
    通过海伦公式判断空间中三点是否共线
    :param Points:点的坐标<numpy.array>
    :return:不共线返回0,反之返回1
    '''
    edge_A = Distance(Points[0], Points[1])
    edge_B = Distance(Points[1], Points[2])
    edge_C = Distance(Points[2], Points[0])

    area = 0.5*(edge_A + edge_B + edge_C);
    if (area*(area - edge_A)*(area - edge_B)*(area - edge_C) >0):
        return 0; #面积大于零 就是一个三角形 三点不共线
    else:
        return 1;			

def Distance(Acoord,Bcoord):
    '''
    计算任意两点间距（笛卡尔坐标系）
    :param Acoord,Bcoord:A、B两点的坐标
    :return:两点间距
    '''
    return math.sqrt((Acoord[0]-Bcoord[0])**2+(Acoord[1]-Bcoord[1])**2+(Acoord[2]-Bcoord[2])**2)

def NormalVectorSurface(Points):
    '''
    计算平面法向量,并归一化
    :param Points:平面上点的坐标
    :return:返回平面法向量
    '''
    if PointOnLine(Points) == 0: #判断点是否共线
        Acoord = Points[0]
        Bcoord = Points[1]
        Ccoord = Points[2]
        NVS_a = (Bcoord[1]-Acoord[1])*(Ccoord[2]-Acoord[2])-(Ccoord[1]-Acoord[1])*(Bcoord[2]-Acoord[2]) #a=(y2-y1)*(z3-z1)-(y3-y1)*(z2-z1)
        NVS_b = (Bcoord[2]-Acoord[2])*(Ccoord[0]-Acoord[0])-(Ccoord[2]-Acoord[2])*(Bcoord[0]-Acoord[0]) #b=(z2-z1)*(x3-x1)-(z3-z1)*(x2-x1)
        NVS_c = (Bcoord[0]-Acoord[0])*(Ccoord[1]-Acoord[1])-(Ccoord[0]-Acoord[0])*(Bcoord[1]-Acoord[1]) #c=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
        a = NVS_a/math.sqrt(NVS_a**2+NVS_b**2+NVS_c**2)
        b = NVS_b/math.sqrt(NVS_a**2+NVS_b**2+NVS_c**2)
        c = NVS_c/math.sqrt(NVS_a**2+NVS_b**2+NVS_c**2)
        normalvector = (a,b,c)
        return normalvector
    else:
        print('Points Coord mistake')

def Top_site(Points):
    top_site = Points
    return top_site

def Bridge_site(Points):
    point1 = Points[0]
    point2 = Points[1]
    bridge_site_x = 0.5*(point1[0]+point2[0])
    bridge_site_y = 0.5*(point1[1]+point2[1])
    bridge_site_z = 0.5*(point1[2]+point2[2])
    bridge_site = (bridge_site_x,bridge_site_y,bridge_site_z)
    return bridge_site

def Hollow_site(Points):
    PointA = Points[0]
    PointB = Points[1]
    PointC = Points[2]
    hollow_site_x = (PointA[0]+PointB[0]+PointC[0])/3
    hollow_site_y = (PointA[1]+PointB[1]+PointC[1])/3
    hollow_site_z = (PointA[2]+PointB[2]+PointC[2])/3
    hollow_site = (hollow_site_x,hollow_site_y,hollow_site_z)
    return hollow_site