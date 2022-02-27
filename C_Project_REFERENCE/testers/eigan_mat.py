import numpy as np
from numpy import linalg as LA

np.random.seed(0)


def generate_symmetric_mat_to_csv(csv_name, dim):
    # generates pseudo-random symmetric matrix sized dim x dim and
    # writes it to a csv file named csv_name_data.csv
    mat = np.random.rand(dim, dim)
    mat = np.tril(mat) + np.tril(mat, -1).T
    np.savetxt(csv_name + "_data.csv", mat, delimiter=",", fmt="%.4f")
    return mat


def eigan_vals_vecs_to_csv(csv_name, mat):
    # computes eigan-vectors and eigan-values, and wrties them to
    # csv_name_eigan.csv the eigan-value first row and the eigan-vectors as columns
    values, vectors = LA.eig(mat)
    np.savetxt(csv_name + "_eigan.csv", [values], delimiter=",", fmt="%.4f")
    with open(csv_name + "_eigan.csv", "a") as f:
        np.savetxt(f, vectors, delimiter=",", fmt="%.4f")


def generate_data_and_eigan_data(csv_name, dim):
    eigan_vals_vecs_to_csv(
        csv_name,
        generate_symmetric_mat_to_csv(csv_name, dim))


def main():
    generate_data_and_eigan_data("data1", 4)
    #generate_data_and_eigan_data("data2", 40)


if __name__ == "__main__":
    main()
