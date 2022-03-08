import math
import os
import numpy as np
from math import e
import os.path

C_file_name = "spkmeans"

def get_vectors(filename):
    with open(filename,"r") as f:
        vectors = [[float(value) for value in line.split(',')] for line in f.readlines() if line != []]
    return vectors


def build_weight_matrix(vectors):
    N_size = vectors.shape[0]
    result = np.fromfunction(lambda i, j: e ** (calc_norm(vectors[i], vectors[j])/(-2)), (N_size, N_size), dtype=int)
    return result - np.eye(N_size)


def build_diagonal_matrix(vectors):
    wMatrix = build_weight_matrix(vectors)
    N_size = wMatrix.shape[0]
    result = np.zeros((N_size, N_size))
    for i in range(len(wMatrix)):
        result[i][i] = np.sum(wMatrix[i])
    return result


def build_lnorm_matrix(vectors):
    wMatrix = build_weight_matrix(vectors)
    dMatrix = build_diagonal_matrix(vectors)
    for i in range(len(dMatrix)):
        dMatrix[i][i] = 1 / np.sqrt(dMatrix[i][i])
    return np.eye(dMatrix.shape[0]) - dMatrix @ wMatrix @ dMatrix


def build_jacobi_matrix(lmatrix):
    epsilon = 0.001
    jacobi_values = lmatrix
    n_size = len(jacobi_values)
    jacobi_vectors = np.eye(n_size)
    cur_f_norm = calc_off(jacobi_values)
    for iter in range(100):
        i, j = find_off_diag_max(jacobi_values)
        psi = (jacobi_values[j][j] - jacobi_values[i][i]) / (2 * jacobi_values[i][j])
        t = (-1 if psi < 0 else 1) / (abs(psi) + math.sqrt(psi ** 2 + 1))
        c = 1 / math.sqrt(1 + t ** 2)
        s = t * c
        P_matrix = build_P_matrix(n_size,i,j,c,s)
        jacobi_values = P_matrix.transpose() @ jacobi_values @ P_matrix
        jacobi_vectors = jacobi_vectors @ P_matrix
        prev_f_norm = cur_f_norm
        cur_f_norm = calc_off(jacobi_values)
        if abs(cur_f_norm - prev_f_norm) < epsilon:
            break
    return jacobi_values, jacobi_vectors


def build_P_matrix(N, i, j, c, s):
    result = np.eye(N)
    result[i][i] = c
    result[j][j] = c
    result[i][j] = s
    result[j][i] = -s
    return result


def calc_off(matrix):
    result = 0
    for i in range(len(matrix)):
        for j in range(i+1, len(matrix)):
            result += matrix[i][j]**2 + matrix[j][i]**2
    return result


def find_off_diag_max(matrix):
    matrix_max = -1
    max_i = -1
    max_j = -1
    for i in range(len(matrix)):
        for j in range(i+1, len(matrix)):
            if abs(matrix[i][j]) > matrix_max:
                matrix_max = abs(matrix[i][j])
                max_i = i
                max_j = j
    return max_i, max_j


def calc_norm(x1, x2):
    return np.sqrt(np.sum((x1 - x2)**2, axis=2))


def check_equality(matrix1, matrix2):
    epsilon = 1e-6
    if len(matrix1) != len(matrix2):
        return False
    for i in range(len(matrix1)):
        if len(matrix1[i]) != len(matrix2[i]):
            return False
        for j in range(len(matrix1[0])):
            if abs(matrix1[i][j] - matrix2[i][j]) > epsilon:
                return False
    return True


def print_mat(matrix):
    for i in range(len(matrix)):
        for j in range(0,len(matrix[i])-1):
            print(f'{matrix[i][j]:.4f},',end="")
        print(f'{matrix[i][len(matrix[0])-1]:.4f}')


if __name__ == "__main__":
    goals = ["wam", "ddg", "lnorm", "jacobi", "spk"]
    exec = ["python3 spkmeans.py", os.path.join(".", C_file_name)]
    result_file = "tmp.txt"
    for index in range(11):
        vectors_filename = os.path.join(".", "tests", f"test{index}.csv")
        lnorm_filename = os.path.join(".", "tests", f"test{index}_lnorm_output.txt")
        for goal in goals:
            for ex in exec:
                lng = "P" if ex == "python3 spkmeans.py" else "C"
                # k matter only for spk
                if goal == "spk":
                    k_values = [0, 3]
                else:
                    k_values = [0]
                for k_value in k_values:
                    if goal == "jacobi":
                        curr_file = lnorm_filename
                    else:
                        curr_file = vectors_filename
                    os.system(f"{ex} {k_value} {goal} {curr_file} > {result_file}")
                    vectors = np.array(get_vectors(vectors_filename))
                    if goal == "spk":
                        correct_matrix = get_vectors(os.path.join(".", "tests", f"test{index}_{goal}_{k_value}_output_{lng}.txt"))
                    else:
                        if goal == "wam":
                            correct_matrix = build_weight_matrix(vectors)
                        elif goal == "ddg":
                            correct_matrix = build_diagonal_matrix(vectors)
                        elif goal == "lnorm":
                            correct_matrix = build_lnorm_matrix(vectors)
                        elif goal == "jacobi":
                            lmatrix = np.array(get_vectors(lnorm_filename))
                            matrix1, matrix2 = build_jacobi_matrix(lmatrix)
                            correct_matrix = [[matrix1[i][i] for i in range(len(matrix1))]]
                            for i in range(len(matrix2)):
                                correct_matrix.append([matrix2[j][i] for j in range(len(matrix2))])
                        correct_matrix = np.round(correct_matrix, 4).tolist()
                    result_matrix = get_vectors(result_file)
                    res = check_equality(correct_matrix, result_matrix)
                    res_str = "Passed" if res else "Failed"
                    print(f"check file={vectors_filename} \tlanguage={lng} \tk_value={k_value} \tgoal={goal} \tResult={res_str}")

                    if not res:
                        print("Your output:")
                        print_mat(result_matrix)
                        print("correct output:")
                        print_mat(correct_matrix)

