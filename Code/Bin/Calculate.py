# -*- coding: utf-8 -*-
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
#数据提取

def extract_data(data_path, E0, directory):
    '''
    提取计算数据存于system中
    :param data_path:待分析数据文件的保存路径
    :return:
    system{
        ['atoms']:              list,
        ['energy_atoms']:       list
        }
    '''
    import os
    import numpy as np
    from ase.io import vasp
    atoms_ptcluster = list()
    energy_ptcluster = list()
    system = {}
    for dir in directory:
        target_dir = data_path + dir
        data_dirs = os.listdir(target_dir)
        for data_dir in data_dirs:
            if data_dir.isdecimal():
                contcar = target_dir + "/" + data_dir + "/" + "CONTCAR"
                oszicar = target_dir + "/" + data_dir + "/" + "OSZICAR"
                oszicar_data = extract_OSZICAR(oszicar)
                # 吸附能：Adopation_energy = E_ptH - E_pt - 0.5E_H2
                energy_pt = oszicar_data['energy_atom'] - E0[dir] - 0.5 * E0["H2"]
                num_elec_steps = oszicar_data['num_elec_steps']
                # 若吸附能大于4eV,则可认为该吸附体系的结构不合理，略去
                # 若电子步数大于60，则可认为该吸附体系的计算未收敛，略去
                if energy_pt < 4 and energy_pt > -1 and num_elec_steps <= 60:
                    atoms_ptcluster.append(vasp.read_vasp(contcar))
                    energy_ptcluster.append(energy_pt)
    
    system['atoms'] = np.array(atoms_ptcluster,dtype=object)
    system['energy_atoms'] = np.array(energy_ptcluster)
    return system

def extract_data_pro(data_path, E0, directory):
    '''
    提取计算数据存于system中
    :param data_path:待分析数据文件的保存路径
    :return:
    system{
        ['atoms']:              list,
        ['energy_atoms']:       list
        }
    '''
    import os
    import numpy as np
    from ase.io import vasp
    atoms_ptcluster = list()
    energy_ptcluster = list()
    system = {}
    for dir in directory:
        target_dir = data_path + dir
        data_dirs = os.listdir(target_dir)
        for data_dir in data_dirs:
            if data_dir.isdecimal():
                contcar = target_dir + "/" + data_dir + "/" + "CONTCAR"
                oszicar = target_dir + "/" + data_dir + "/" + "OSZICAR"
                oszicar_data = extract_OSZICAR(oszicar)
                # 吸附能：Adopation_energy = E_ptH - E_pt - 0.5E_H2
                energy_pt = oszicar_data['energy_atom'] - E0[dir] - 0.5 * E0["H2"]
                num_elec_steps = oszicar_data['num_elec_steps']
                # 若吸附能大于0.5eV,则可认为该吸附体系的排斥太强，略去
                # 若电子步数大于60，则可认为该吸附体系的计算未收敛，略去
                if energy_pt < 0.5 and energy_pt > -1 and num_elec_steps <= 60:
                    atoms_ptcluster.append(vasp.read_vasp(contcar))
                    energy_ptcluster.append(energy_pt)
    
    system['atoms'] = np.array(atoms_ptcluster,dtype=object)
    system['energy_atoms'] = np.array(energy_ptcluster)
    #system['atoms'] = atoms_ptcluster
    #system['energy_atoms'] = energy_ptcluster
    return system
########################################################################################################################################
# 提取POSCAR文件信息【注释、缩放系数、基矢、元素、元素原子个数、原子总数、坐标类型、笛卡尔坐标（原子位置）、分数坐标（原子位置）】
def extract_POSCAR(file_name):
    '''
    提取POSCAR文件的内容存于system中
    :param file_name: POSCAR文件的保存路径
    :return:
    system{
        ['comment']:            string,
        ['scale_coeff']:        float,
        ['box_coord']:          list,
        ['atom']:               list,
        ['atom_num']:           list,
        ['all_item']:           int,
        ['coord_type']:         string,
        ['coord_Cartesian']:    list,
        ['coord_Direct']:       list 
        }
    '''
    import numpy as np
    file = open(file_name,mode='r')
    content = file.readlines()
    system = {}
    system['comment'] = content[0].strip()
    system['scale_coeff'] = float(content[1])
    for i in range(2,5):
        system.setdefault('box_coord',[]).append(
            [float(x)*system['scale_coeff'] for x in content[i].split()])
    system['atom'] = content[5].strip().split()
    system['atom_num'] = [int(x) for x in content[6].strip().split()]
    system['all_atom'] = sum(system['atom_num'])
    system['coord_type'] = 'Cartesian' if content[7].strip().startswith('C') \
                                       or content[7].strip().startswith('c') \
                                     else 'Direct'
    if system['coord_type'] == 'Cartesian':
        for i in range(8,8+system['all_atom']):
            system.setdefault('coord_Cartesian',[]).append(
                [float(x) for x in content[i].split()[:3]])
            system['coord_Direct'] = \
                np.dot(np.linalg.inv(np.array(system['box_coord']).T),np.array(system['coord_Cartesian']).T).T.tolist()
    else:
        for i in range(8,8+system['all_atom']):
            system.setdefault('coord_Direct',[]).append(
                [float(x) for x in content[i].split()[:3]])
            system['coord_Cartesian'] = \
                np.dot(np.array(system['box_coord']).T,np.array(system['coord_Direct']).T).T.tolist()
    file.close()
    return system 

#提取OSZICAR文件信息【能量、电子步数】
def extract_OSZICAR(file_name):
    '''
    提取OSZICAR文件中信息并保存在system中
    :param file_name:VASP计算得到的OSZICAR文件的保存路径
    :return:
    system{
		['energy_atom']:        list,
        ['num_elec_steps']:     list		
	    }
    '''
    file = open(file_name,mode='r')
    content = file.readlines()
    system = {}
    system['energy_atom'] = float(content[-1].split()[4])
    system['num_elec_steps'] = int(content[-2].split()[1])
    file.close()
    return system

####################################################################################################################################
def min_enclosing_sphere(atoms):#,symbol_list):
    '''
    给定吸附体系的atoms对象,返回最小包围球的半径与球心坐标
    :param atoms:atoms对象:atoms
    :param symbol_list:所需确定的元素:list
    :return:
    system{
        ['radius']
        ['circle_center']
        } 返回最小包围球的半径与球心坐标
    '''
    import sys
    sys.path.append('/home/yzh/Mywork/Mycode/Bin/')
    import Radius
    system = {}
    list_of_point = list()
    for atom in atoms:
        #if atom.symbol in symbol_list:
        #    continue
        point = Radius.Point(atom.position[0], atom.position[1], atom.position[2])
        list_of_point.append(point)
    sphere = Radius.output_Radius(list_of_point)
    system['center'] = sphere[0]
    system['radius'] = sphere[1]
    return system #radius, circle_center

def MES(atoms):
    '''
    给定吸附体系的atoms对象,返回最后一个点到除了最后一个点之外的中心点的距离radius,球心center,最后一个点到其余点距离的最小值md
    :param atoms:atoms对象:atoms
    :return:
    system{
        ['radius']
        ['circle_center']
        ['min_distance']
        } 返回最小包围球的半径与球心坐标
    '''
    import sys
    sys.path.append('/home/zhihengyu/WORK/Code/Bin/')
    import Radius
    system = {}
    list_of_point = list()
    for atom in atoms:
        point = Point(atom.position[0], atom.position[1], atom.position[2])
        list_of_point.append(point)
    sphere = Radius.output_SES(list_of_point)
    system['center'] = sphere[0]
    system['radius'] = sphere[1]
    system['min_distance'] = sphere[2]
    return system 

#计算并返回两点（三维空间）之间的距离
def distance_3D_dot(initial_dot,terminal_dot):
    import math
    return math.sqrt((terminal_dot[0]-initial_dot[0])**2+(terminal_dot[1]-initial_dot[1])**2+(terminal_dot[2]-initial_dot[2])**2)


def CreateWulffstructure(num):
    from ase.cluster.wulff import wulff_construction
    atoms = wulff_construction('Pt', 
                           surfaces=[(1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 1, 1)], 
                           energies=[1.856,1.871,1.488,1.763], 
                           size=num, 
                           structure='fcc', 
                           latticeconstant=3.92)
    return atoms

##################################################################################################################################

def get_basis_sets_by_KMeans(eigenvector,n_clusters):
    '''
	:param eigenvector: 二维numpy数组,每一行(axis=0)表示一个ase.Atoms对象的特征向量
	:param n_clusters: 需要筛选得到的吸附构型的目标个数
	:return: 返回筛选得到作为基组的吸附构型的特征向量，及其索引信息
	'''
    import numpy as np
    from sklearn.cluster import KMeans
    from sklearn.metrics import pairwise_distances_argmin_min
    is_basis_vector = [False for _ in range(len(eigenvector))]
    basis_sets = list()
    k_means = KMeans(n_clusters=n_clusters, random_state=1,n_init='auto').fit(eigenvector)
    basis_index, _ = pairwise_distances_argmin_min(k_means.cluster_centers_, eigenvector)
    for index in basis_index:
        is_basis_vector[index] = True
        basis_sets.append(eigenvector[index])
    basis_sets = np.array(basis_sets)
    return basis_sets, basis_index

################################################################################################################################

def get_represent_vector(list_of_feature_vector, basis_sets, power=3):
    '''
	:param list_of_feature_vector: 二维numpy数组,每一行表示一个待求ase.Atoms对象的频率谱
	:param basis_sets: 二维numpy数组,每一行(axis=0)表示一个经过筛选的ase.Atoms对象的频率谱
	:return: 一维numpy数组,待求ase.Atoms对象与基组中的全体基矢进行相似度评估得到的相似度向量,用于机器学习模型输入
    '''
    import numpy as np
    list_of_feature_vector = np.divide(list_of_feature_vector.T, np.sqrt(np.diag(np.dot(list_of_feature_vector, list_of_feature_vector.T)))).T
    basis_sets = np.divide(basis_sets.T, np.sqrt(np.diag(np.dot(basis_sets, basis_sets.T)))).T
    pruduct = np.dot(list_of_feature_vector, basis_sets.T)
    return np.power(pruduct, power)

# 得到相似度矩阵的改良方法，充分使用numpy内置的并行运算
def get_similarity_matrix(eigenvector, power=3):
    '''
	:param eigenvector: 二维numpy数组,每一行(axis=0)表示一个ase.Atoms对象的特征向量
	:param power: 计算相似度过程中取的阶数,为一个正数,默认为3
	:return: 二维numpy数组,matrix[row, column]表示第row个与第column个ase.Atoms对象的相似度
	'''
    import numpy as np
    pruduct = np.dot(eigenvector, eigenvector.T)
    diag = np.sqrt(np.diag(pruduct))
    similarity_matrix = np.power(np.divide(np.divide(pruduct, diag).T, diag), power)
    return similarity_matrix


def get_surface_relations(points):
    '''
    获取点集构成的表面关系。
    :param points: (numpy.ndarray) 点集，形状为(N, 3)。
    :return:
        Tuple[numpy.ndarray, List[List[Tuple[int, int]]], List[List[Tuple[int, int]]]]:
            返回一个元组，包含三个元素：
            - 表示点集的数组；
            - 表示每个表面所对应的点的列表，其中每个元素表示一个表面对应的三个点的索引；
            - 表示每个表面对应的边的列表，其中每个元素表示一个边对应的两个点的索引。
    '''
    from scipy.spatial import Delaunay
    # 利用Delaunay三角化计算表面
    tri = Delaunay(points)
    faces = tri.convex_hull

    # 对每个点寻找所在表面并进行汇总
    point2face = {}
    for i, face in enumerate(faces):
        for p in face:
            if p not in point2face:
                point2face[p] = [i]
            else:
                break

    # 整理得到每个点所在的第一个表面
    point2face_first = {p: f[0] for p, f in point2face.items()}

    # 整理得到每个表面包含的点的编号列表
    face_points = []
    for i, face_idx in enumerate(faces):
        non_empty_face = [p for p in face_idx if point2face_first[p] == i]
        if non_empty_face:
            face_points.append(non_empty_face)


    # 构建每条边到面的映射关系，只保留每条边对应的第一组表面
    edge_to_face = {}
    for i, triangle in enumerate(faces):
        for j in range(3):
            edge = tuple(sorted([triangle[j], triangle[(j+1)%3]]))
            if edge not in edge_to_face:
                edge_to_face[edge] = i

    # 构建面到边的映射关系
    face_to_edges = {}
    for edge, face in edge_to_face.items():
        if face in face_to_edges:
            face_to_edges[face].append(edge)
        else:
            face_to_edges[face] = [edge]

    # 将边按照所在面汇总
    face_line = []
    for edges in face_to_edges.values():
        face_line.append(edges)

    return points, face_points, face_line



#####################################################################################################################################
def farthest_point_sample(xyz, npoint):
    '''
    Input:
        xyz: pointcloud data, [B, N, 3]
        npoint: number of samples
    Return:
        centroids: sampled pointcloud index, [B, npoint]
    '''
    import numpy as np
    xyz = xyz.transpose(0,2,1)
    B, N, C = xyz.shape
    
    centroids = np.zeros((B, npoint))    # 采样点矩阵（B, npoint）
    distance = np.ones((B, N)) * 1e10                       # 采样点到所有点距离（B, N）

    batch_indices = np.arange(B)        # batch_size 数组
    
    barycenter = np.sum((xyz), 1)                                    #计算重心坐标 及 距离重心最远的点
    barycenter = barycenter/xyz.shape[1]
    barycenter = barycenter.reshape(B, 1, C)   #numpy中的reshape相当于torch中的view

    dist = np.sum((xyz - barycenter) ** 2, -1)
    farthest = np.argmax(dist,1)                                     #将距离重心最远的点作为第一个点,这里跟torch.max不一样

    for i in range(npoint):
        #print("-------------------------------------------------------")
        #print("The %d farthest pts %s " % (i, farthest))
        centroids[:, i] = farthest                                      # 更新第i个最远点
        centroid = xyz[batch_indices, farthest, :].reshape(B, 1, C)        # 取出这个最远点的xyz坐标
        dist = np.sum((xyz - centroid) ** 2, -1)                     # 计算点集中的所有点到这个最远点的欧式距离,-1消掉了xyz那个维度
        #print("dist    : ", dist)
        mask = dist < distance
        #print("mask %i : %s" % (i,mask))
        distance[mask] = dist[mask]                                     # 更新distance,记录样本中每个点距离所有已出现的采样点（已采样集合中的点）的最小距离
        #print("distance: ", distance)

        farthest = np.argmax(distance, -1)                           # 返回最远点索引
 
    return centroids

#########################################################################################################################