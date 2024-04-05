# -*- coding: utf-8 -*-
import math
import random

class Sphere:
    def __init__(self, center=None, radius=0):
        self.center = center
        self.radius = radius

class Point:
    """
    创建Point类,保存三维空间中点的x,y,z坐标
    """
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
    
    def __str__(self):
        return "(" + str(self.x) + ", " + str(self.y) + ", " + str(self.z) + ")"
    
def dist(p1, p2):
    return math.sqrt((p1.x - p2.x)**2 + (p1.y - p2.y)**2 + (p1.z - p2.z)**2)

def is_close(value1, value2, tolerance=1e-9):
    return math.isclose(value1, value2, abs_tol=tolerance)

def createSphere0(sphere):
    sphere.center = Point(0, 0, 0)
    sphere.radius = 0

def createSphere1(p1, sphere):
    sphere.center = p1
    sphere.radius = 0

def createSphere2(p0, p1, sphere):
    sphere.center = Point((p0.x + p1.x) / 2, (p0.y + p1.y) / 2, (p0.z + p1.z) / 2)
    sphere.radius = dist(p0, p1) / 2

def createSphere3(p0, p1, p2, sphere):
    a = dist(p0, p1)
    b = dist(p1, p2)
    c = dist(p2, p0)
    if a >= b and a >= c:
        sphere.center = Point((p0.x + p1.x) / 2, (p0.y + p1.y) / 2, (p0.z + p1.z) / 2)
        sphere.radius = a / 2
    elif b >= c:
        sphere.center = Point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2, (p1.z + p2.z) / 2)
        sphere.radius = b / 2
    else:
        sphere.center = Point((p2.x + p0.x) / 2, (p2.y + p0.y) / 2, (p2.z + p0.z) / 2)
        sphere.radius = c / 2

def createSphere4(p0, p1, p2, p3, sphere):
    a = dist(p0, p1)
    b = dist(p0, p2)
    c = dist(p0, p3)
    d = dist(p1, p2)
    e = dist(p1, p3)
    f = dist(p2, p3)

    if a >= b and a >= c and a >= d and a >= e and a >= f:
        createSphere2(p0, p1, sphere)
    elif b >= c and b >= d and b >= e and b >= f:
        createSphere2(p0, p2, sphere)
    elif c >= d and c >= e and c >= f:
        createSphere2(p0, p3, sphere)
    elif d >= e and d >= f:
        createSphere2(p1, p2, sphere)
    elif e >= f:
        createSphere2(p1, p3, sphere)
    else:
        createSphere2(p2, p3, sphere)

def createSphere(support, numPoint, sphere, visited):
    if numPoint == 0:
        createSphere0(sphere)
    elif numPoint == 1:
        createSphere1(support[0], sphere)
    elif numPoint == 2:
        createSphere2(support[0], support[1], sphere)
    elif numPoint == 3:
        createSphere3(support[0], support[1], support[2], sphere)
    else:
        createSphere4(support[0], support[1], support[2], support[3], sphere)

    for i in range(numPoint):
        #if dist(support[i], sphere.center) > sphere.radius:
        if not visited[i] and not is_close(dist(support[i], sphere.center), sphere.radius):
            visited[i] = True
            newSupport = []
            newNumPoint = 0
            for j in range(i + 1):
                if not visited[j]:
                    newSupport.append(support[j])
                    newNumPoint += 1
            createSphere(newSupport, newNumPoint, sphere, visited)
            visited[i] = False

def smallestEnclosingSphere(points):
    numPoint = len(points)
    support = []
    for i in range(numPoint):
        support.append(points[i])
    
    sphere = Sphere()
    visited = [False] * numPoint
    createSphere(support, numPoint, sphere, visited)
    return sphere.center, sphere.radius

def output_SmallestEnclosingSphere(points):
    center, radius = smallestEnclosingSphere(points)
    return center, radius

# 清除缓存
#import gc
#gc.collect()