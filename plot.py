import numpy as np
import matplotlib.pyplot as plt
import json
from sklearn.decomposition import PCA

locks = {'A':[],'B':[],'C':[],'D':[],'E':[],'F':[],'G':[],'H':[],'I':[],'J':[],'.':[],'R':[]}
keys = []
rows = []
with open('test_data/VD1iter1time1001.dat') as f:
    for line in f:
        data_line = line.split()
        position = json.loads(data_line[0])
        chromosome = json.loads(data_line[1])
        phages = json.loads(data_line[2])

        for gene in chromosome:
            locks[gene['type']].append(gene['lock'])
        for phage in phages:
            keys.append(phage[0]['key'])

X = []
labels = []

for label, bitstrings in locks.items():
    X.extend(bitstrings)
    labels.extend([label] * len(bitstrings))

X = np.array(X)

pca = PCA(n_components=1)
X_pca = pca.fit_transform(X)

plt.scatter(X_pca[:, 0], X_pca[:, 1])

plt.show()
