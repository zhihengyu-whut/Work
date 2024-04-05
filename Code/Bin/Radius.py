# -*- coding: utf-8 -*-
import numpy as np
from scipy.spatial.distance import euclidean
from scipy.spatial.distance import cdist
 
def get_SES(points):
    coordinates = np.array([[point.x, point.y, point.z] for point in points[:-1]])  # 提取其他点的坐标
    center = np.mean(coordinates, axis=0)  # 计算除了最后一个点之外的中心点
    last_point_coordinate = np.array([points[-1].x, points[-1].y, points[-1].z])  # 提取最后一个点的坐标
    radius = euclidean(center, last_point_coordinate)  # 计算中心点和最后一个点之间的欧氏距离
    
    distances = cdist([last_point_coordinate], coordinates)  # 计算最后一个点到其他点的距离
    min_distance = np.min(distances)  # 获取最小距离
    
    return center, radius, min_distance


def output_SES(points):
    center, radius, min_distancer = get_SES(points)
    return center, radius, min_distancer
