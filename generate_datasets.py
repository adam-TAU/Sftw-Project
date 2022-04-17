import numpy as np
import matplotlib.pyplot as plt

from sklearn import datasets

def save_dataset(dataset, filename):
    lines = map(stringify_row, dataset)
    with open(filename + ".csv", "w") as f:
        f.writelines(lines)

def stringify_row(row):
    row = map(str, row)
    return ",".join(row) + "\n"

np.random.seed(0)

# ============
# Generate datasets. We choose the size big enough to see the scalability
# of the algorithms, but not too big to avoid too long running times
# ============
n_samples = 100
noisy_circles = datasets.make_circles(n_samples=n_samples, factor=0.3, noise=0.05)
noisy_moons = datasets.make_moons(n_samples=n_samples, noise=0.05)
blobs = datasets.make_blobs(n_samples=n_samples, random_state=8)

# Anisotropicly distributed data
random_state = 170
X, y = datasets.make_blobs(n_samples=n_samples, random_state=random_state)
transformation = [[0.6, -0.6], [-0.4, 0.8]]
X_aniso = np.dot(X, transformation)
aniso = (X_aniso, y)

# blobs with varied variances
varied = datasets.make_blobs(
        n_samples=n_samples, cluster_std=[1.0, 2.5, 0.5], random_state=random_state
        )

plt.scatter(blobs[0][:,0], blobs[0][:,1])
plt.savefig("blobs.png")

# Save the datasets to disk
save_dataset(10 * blobs[0], "blobs")

# Save 1000 * 10 points to disk
save_dataset(datasets.make_blobs(n_samples=1000, centers=3, n_features=10, shuffle=True, random_state=31)[0], "1000_blobs_10.csv")

