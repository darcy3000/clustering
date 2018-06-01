import numpy as np
# import matplotlib.pyplot as plt
from collections import OrderedDict as OD
import pickle
import networkx as nx
import random
#import seaborn as sns
import scipy.cluster.hierarchy as sch
import heapq

# Reading the file
file = open("data.txt", "r")
data = file.read()
data.replace('\n', "")
x_data = data.split('>')

# Special case
del x_data[0]

points = OD()  # Storing all points where key is chr_*** and values is the entire string
name_to_id = {}  # Hashing from name to id
id_to_name = {}  # Hashing from id to name

def distance(cluster_1, cluster_2):
    len1 = len(cluster_1)
    len2 = len(cluster_2)
    dist = 0
    for i in cluster_1:
        for j in cluster_2:
            dist += proximity_matrix[i, j]

    dist /= (len1 * len2)
    return dist

for i in x_data:
    sequence = i.split('\n')
    val = ""
    key = sequence[0]
    for j in range(1, len(sequence)):
        val = val + sequence[j]
    points[key] = val

ids = np.arange(0, len(x_data))
dendro_tree = np.arange(0, len(x_data))  # Will be used to populate the z matrix
cluster_no = (2 * len(x_data)) - 3  # will be used to keep track while renaming of clusters
Z = np.zeros((len(x_data) - 1, 4))
hash={}

for ind, i in enumerate(points):
    name_to_id[i] = ind
    id_to_name[ind] = i

proximity_matrix = pickle.load(open('proximity.pkl', 'rb'))


def split(C):
    A = C
    B = []
    a = np.zeros((len(C)))
    index_to_id = []
    for ind, i in enumerate(C):
        a[ind] = np.sum(proximity_matrix[i, C])
        index_to_id.append(i)

    a /= (len(C) - 1)
    max_ind = np.argmax(a)
    B.append(A[max_ind])
    del A[max_ind]
    while (len(A) > 1):
        a_dis = np.zeros((len(A)))
        b_dis = np.zeros((len(A)))
        for ind, i in enumerate(A):
            a_dis[ind] = np.sum(proximity_matrix[i, A])
            a_dis = a_dis / (len(A) - 1)
            b_dis[ind] = np.sum(proximity_matrix[i, B])
            b_dis = b_dis / (len(B))
            a_dis = a_dis - b_dis

        max_ind = np.argmax(a_dis)
        if (A[max_ind] > 0):
            B.append(A[max_ind])
            del A[max_ind]
        else:
            break

    a_max = 0
    b_max = 0

    '''a_max = distance(A,A)
    b_max = distance(B,B)'''

    for i in A:
        for j in A:
            if a_max < proximity_matrix[i, j]:
                a_max = proximity_matrix[i, j]

    for i in B:
        for j in B:
            if b_max < proximity_matrix[i, j]:
                b_max = proximity_matrix[i, j]

    return ((-1 * a_max, A), (-1 * b_max, B))


ini_arr = list(range(len(x_data)))
heap = []
heapq.heappush(heap, (100, ini_arr))



for k in range(len(x_data)-1):
    C = heapq.heappop(heap)
    lenC = len(C[1])
    if (k > 0):
        x, y = hash[tuple(C[1])]
        Z[x, y] = cluster_no
        cluster_no -= 1
    A, B = split(C[1])
    heapq.heappush(heap, A)
    heapq.heappush(heap, B)
    dist = distance(A[1],B[1])
    N = len(x_data)

    print(k)
    #print(A[1])
    #print(B[1])


    if (len(A[1])==1):
        Z[N-2-k, 0] = A[1][0]
    else:
        hash[tuple(A[1])] = (N - 2 - k, 0)

    if (len(B[1])==1):
        Z[N-2-k, 1] = B[1][0]
    else:
        hash[tuple(B[1])] = (N - 2 - k, 1)
    Z[N - 2 - k, 2] = dist
    Z[N - 2 - k, 3] = lenC
    '''print(A[0], A[1])
    print(B[0], B[1])
    print()'''

sch.dendrogram(Z)


