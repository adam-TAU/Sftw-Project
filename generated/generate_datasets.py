import numpy as np
import matplotlib.pyplot as plt
import sys

from sklearn import datasets

def save_dataset(dataset, filename):
    X, y = dataset
    lines = map(stringify_row, X)
    with open(filename + ".csv", "w") as f:
        f.writelines(lines)
    
    k = max(y) + 1
    for cluster in range(k):
        points = np.array([X[i,:] for i in range(X.shape[0]) if y[i] == cluster])
        plt.scatter(points[:,0], points[:,1])
    plt.savefig(filename + ".png")

def stringify_row(row):
    row = map(str, row)
    return ",".join(row) + "\n"

def circles(n_samples, other):
    other = 0.3 if other is None else other
    return datasets.make_circles(n_samples=n_samples, factor=other, noise=0.05)

def moons(n_samples, other):
    other = 0.05 if other is None else other
    return datasets.make_moons(n_samples=n_samples, noise=other)

def blobs(n_samples, other):
    other = 1.0 if other is None else other
    return datasets.make_blobs(n_samples=n_samples, cluster_std=other)

def aniso(n_samples, other):
    other = 1.0 if other is None else other
    X, y = datasets.make_blobs(n_samples=n_samples, cluster_std=other)
    transformation = [[0.6, -0.6], [-0.4, 0.8]]
    X_aniso = np.dot(X, transformation)
    return X_aniso, y

def varied(n_samples, _):
    return datasets.make_blobs(n_samples=n_samples, cluster_std=[1.0, 2.5, 0.5])
    

if __name__ == "__main__":
    args = sys.argv
    if len(args) < 3:
        print("Help: python3 generate_datasets.py <kind> <n_samples> [scale] [other]\n")
        print("- kind:      any of the following: circles/moons/blobs/aniso/varied")
        print("- n_samples: number of points. must be a positive integer")
        print("- scale:     scale factor to scale every output point by")
        print("- other:     misc (float) option. changes behavior differently depending on the kind\n")
        print("The output will be saved as `kind.csv`, along with an image `kind.png`.")

    else:
        argc = len(args)
        kind = args[1]
        try:
            n_samples = int(args[2])
            scale = float(args[3]) if argc > 3 else 1.0
            other = float(args[4]) if argc > 4 else None
        except Exception:
            print("Invalid input. Verify that `n_samples`, `scale` and `other` are valid numbers.")
            exit()

        if kind not in ["circles", "moons", "blobs", "aniso", "varied"]:
            print(f"Invalid kind: {kind}")
            exit()
        
        if n_samples < 1:
            print(f"Invalid number of samples: {n_samples}")
        
        func = {"circles":circles, "moons":moons, "blobs":blobs, "aniso":aniso, "varied":varied}[kind]
        data = func(n_samples, other)
        data = (data[0] * scale, data[1])

        save_dataset(data, kind)
