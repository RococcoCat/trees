import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt

# calculate Jukes-Cantor distance
def jc_dist(s1, s2):
    length = min(len(s1), len(s2))
    pt = [i != j for i, j in zip(s1[:length], s2[:length])].count(True)
    dist = -3/4 * np.log(1 - 4/3 * pt / length)
    return dist

# Ultrametric tree
def is_ultrametric(df):
    cols = df.columns
    for i in range(len(cols)):
        for j in range(i+1, len(cols)):
            for k in range(j+1, len(cols)):
                a = (df.iloc[i][j], df.iloc[j][k], df.iloc[i][k])
                if len(set(a)) != 2:
                    print(a, i, j, k)
                    return False
    return True

def upgma(df, out_file, d):
    # recursive
    # base case: only 2 taxa left
    if len(df.columns) == 2:
        a,b = df.stack().idxmin()
        dist = 0.5 * df.loc[a][b]
        with open(out_file, "w") as out_file:
            if d == True:
                out_file.write(f"({a}:{dist}, {b}:{dist}) \n")
            else:
                out_file.write(f"({a}, {b}) \n")
        return (', '.join(df.columns))
    # recursive case: more than 2 taxa left
    else:
        # find the 3 taxa that are most closely related
        df = df.replace(0, np.NaN)
        a,b = df.stack().idxmin()
        dist = 0.5 * df.loc[a][b]
        if d == True: col_name = f"({a}:{dist}, {b}:{dist})" 
        else: col_name = f"({a}, {b})"       
        ab = []
        for i in df.columns:
            if i != a and i != b:
                dist = 0.5 * (df.loc[a][i] + df.loc[b][i])
                ab.append(dist)
            elif i != b:
                ab.append(np.NaN)
        df = df.drop(b, axis=1)
        df = df.drop(b, axis=0)
        df[a] = ab
        df.loc[[a]] = ab
        df = df.rename(columns={a: col_name})
        df = df.rename(index={a: col_name})
        return upgma(df, out_file, d)
    
# Neighbor joining
def neighbor_joining(df, out_file, d=True):
    if len(df.columns) == 2:
        a,b = df.stack().idxmin()
        dist = abs(0.5 * df.loc[a][b])
        with open(out_file, "w") as w:
            if d == True:
                w.write(f"({a}:{dist}, {b}:{dist}) \n")
            else:
                w.write(f"({a}, {b}) \n")
        return (', '.join(df.columns))
    else:
        u = {}
        for i in df.columns:
            u[i] = (1/(len(df.columns)-2)) * list(df[i].isnull()).count(False)
        df = df.replace(0, np.NaN)
        a,b = df.stack().idxmin()
        dist = abs(0.5 * df.loc[a][b])
        if d == True: col_name = f"({a}:{dist}, {b}:{dist})" 
        else: col_name = f"({a}, {b})"      
        for i in df.columns:
            for j in df.columns:
                if i != j:
                    df[i][j] = df[i][j] - u[i] - u[j]
        df = df.drop(b, axis=1)
        df = df.drop(b, axis=0)        
        df = df.rename(columns={a: col_name})
        df = df.rename(index={a: col_name})
        return neighbor_joining(df, out_file, d)

def display_tree(df):
    # Calculate linkage matrix
    Z = linkage(df)

    # Create dendrogram
    dendro = dendrogram(Z)
    dendrogram(Z, labels=df.columns, orientation='left')

    plt.title('Phylogenetic Tree')
    plt.xlabel('Species')
    plt.ylabel('Distance')
    plt.show()
