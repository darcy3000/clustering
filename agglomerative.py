import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict as OD
import pickle
import networkx as nx
import random
#import seaborn as sns
import scipy.cluster.hierarchy as sch
#Reading the file
file = open("data.txt","r")
data = file.read()
data.replace('\n',"")
x_data = data.split('>')



#Special case
del x_data[0]

points = OD()#Storing all points where key is chr_*** and values is the entire string
name_to_id = {}#Hashing from name to id
id_to_name ={}#Hashing from id to name

for i in x_data:
    sequence = i.split('\n')
    val = ""
    key = sequence[0]
    for j in range(1,len(sequence)):
        val=val+sequence[j]
    points[key]=val

ids = np.arange(0,len(x_data))
dendro_tree = np.arange(0,len(x_data))#Will be used to populate the z matrix
cluster_no = len(x_data)#will be used to keep track while renaming of clusters
Z = np.zeros((len(x_data)-1,4))


for ind,i in enumerate(points):
    name_to_id[i] = ind
    id_to_name[ind] = i


def similarity(a,b):
    gap = 2
    substitution = 1
    match =0
    opt =np.zeros((len(a)+1,len(b)+1))
    for i in range(1,len(a)+1):
        opt[i,0] = opt[i-1,0] + gap

    for j in range(1,len(b)+1):
        opt[0,j] = opt[0,j-1] + gap

    for i in range(1,len(a)+1):
        for j in range(1,len(b)+1):
            scoreDiag = opt[i-1,j-1] + (match if a[i-1]==b[j-1] else substitution)
            scoreLeft = opt[i,j-1]+gap
            scoreUp = opt[i-1,j] + gap

            opt[i,j] = min(scoreDiag,scoreLeft,scoreUp)

    #print(opt)
    return opt[len(a),len(b)]

def find(a):
    global ids
    global name_to_id
    start = name_to_id[a]
    while(start != ids[start]):
        ids[start]=ids[ids[start]]
        start = ids[start]
    return start

def merge(a,b):
    global ids
    global name_to_id
    root_a = find(a)
    root_b = find(b)
    ids[root_a] = root_b


proximity_matrix = np.zeros((len(x_data),len(x_data)))
#Calculating the lower triangular proximity matrix
"""
for in1,i in enumerate(points):
    for in2,j in enumerate(points):
        if(in2>=in1):
            break
        print(in1)
        proximity_matrix[name_to_id[i],name_to_id[j]] = similarity(points[i],points[j])
        proximity_matrix[name_to_id[j],name_to_id[i]] = proximity_matrix[name_to_id[i],name_to_id[j]]

pickle.dump(proximity_matrix,open('proximity.pkl','wb'))
"""

proximity_matrix = pickle.load(open('proximity.pkl','rb'))
#Building the denrodgram
"""
dendro = proximity_matrix
sns.clustermap(dendro,cmap = 'mako',robust = True)
"""

#Special case there exist 2 distinct points such that there mismatch is zero so clubbing them into one cluster first

merge(id_to_name[119],id_to_name[120])
Z[0,0] = 119
Z[0,1] = 120
Z[0,2] = 0
Z[0,3] = 2
dendro_tree[119] = cluster_no
dendro_tree[120] = cluster_no
cluster_no+=1


for k in range(1,len(x_data)-1):
    cluster_1 = 0
    cluster_2 = 0
    min_dist = 999999999999
    print(k)
    for i in range(len(x_data)):
        for j in range(len(x_data)):
            if(proximity_matrix[i,j]<min_dist and proximity_matrix[i,j]!=0):
                cluster_1 = i
                cluster_2 = j
                min_dist = proximity_matrix[i,j]

    merge(id_to_name[cluster_1],id_to_name[cluster_2])


    cluster_1_root = find(id_to_name[cluster_1])
    cluster_2_root = find(id_to_name[cluster_2])

    #Populating the Z-matrix

    Z[k,0] = dendro_tree[cluster_1]
    Z[k,1] = dendro_tree[cluster_2]
    Z[k,2] = proximity_matrix[cluster_1,cluster_2]
    count_of_cluster_points = 0
    for ind,p in enumerate(ids):
        if (find(id_to_name[ind])==cluster_1_root or find(id_to_name[ind])==cluster_2_root):
            count_of_cluster_points+=1
            dendro_tree[ind]=cluster_no
    cluster_no+=1
    Z[k,3] = count_of_cluster_points

    #Updating the proximity matrix
    #Making entries in the same cluster zero
    cluster_two = [p for p in range(len(x_data)) if (cluster_2_root == find(id_to_name[p]))]
    for x in cluster_two:
        for y in cluster_two:
            proximity_matrix[x,y]=0

    for i in range(len(x_data)):
        root_i = find(id_to_name[i])
        cluster_one = [p for p in range(len(x_data))if(root_i==find(id_to_name[p]))]
        #If they are already in same cluster
        if(root_i == cluster_2_root):
            continue
        num = 0
        check_max = -9999
        for x in cluster_one:
            for y in cluster_two:
                check_max=max(check_max,proximity_matrix[x,y])

        val_to_insert = check_max

        for x in cluster_one:
            for y in cluster_two:
                proximity_matrix[x,y] = val_to_insert
                proximity_matrix[y,x] = val_to_insert




#pickle.dump(ids,open('id_5.pkl','wb'))
sch.dendrogram(Z)

ids = pickle.load(open("id_5.pkl",'rb'))

G = nx.DiGraph()
nodes = np.arange(len(x_data))
edges = []
for i in range(len(x_data)):
    edges.append([i,ids[i]])

G.add_nodes_from(nodes)
G.add_edges_from(edges)


color_list =['r','b','g','c','m']
s = {find(id_to_name[i]) for i in ids}
color_map = {}

for i in s:
    color_map[i]=random.choice(color_list)

color_nodes = []
for i in range(len(x_data)):
    color_nodes.append(color_map[find(id_to_name[i])])


pos = nx.spring_layout(G)
nx.draw_networkx(G,pos,node_size = 100,node_color=color_nodes,labels={})

