import numpy as np
import matplotlib.pyplot as plt

# График для проверки расстояний в матрице расстояний
def check_matrix(dist_matrix, read_length):
    matrix_length = len(dist_matrix)
    dists = np.zeros(read_length)
    for line_number in range(matrix_length):
        for dist_number in range(line_number + 1, matrix_length):
            dist = dist_matrix[line_number, dist_number]
            dists[dist] += 1
    y = np.zeros(read_length)
    for i in range(read_length):
        y[i] = i
    plt.plot(y, dists)
    plt.show()

def find_peak(dist_matrix, read_length):
    matrix_length = len(dist_matrix)
    dists = np.zeros(read_length)
    peak_dists = np.zeros(2)
    for line_number in range(matrix_length):
        for dist_number in range(line_number + 1, matrix_length):
            dist = dist_matrix[line_number, dist_number]
            dists[dist] += 1
    middle_value = int(len(dists) / 2) + 1
    little_dists = np.zeros(middle_value)
    big_dists = np.zeros(middle_value)
    for i in range(middle_value):
        little_dists[i] = dists[i]
    for j in range(middle_value, len(dists)):
        big_dists[j - middle_value] = dists[j]
    peak_dists[0] = np.argmax(little_dists)
    peak_dists[1] = np.argmax(big_dists)
    return peak_dists
