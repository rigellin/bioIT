from __future__ import print_function
import numpy as np
import clustering_functions as clust
import warnings
warnings.filterwarnings('ignore')

# Объявление переменных
# Массив данных, считанных из файла. В нём содержатся и значения нуклеотидов, и вероятности
data = []
# Массив данных, полученных из data[]. В нём содержатся значения нуклеотидов
value = []
# Массив данных, полученных из data[]. В нём содержатся значения вероятностей
probability = []
# Массив оригинальных аллелей
alleles = []
# Матрица расстояний
dist_matrix_list = []
# Данные, конвертированные из строк в последовательность байтов
byte_value = []
# Количество кластеров
number_of_clusters_mas = [3, 4, 5, 6]

# Считывание данных из файла. На вход принимаются два типа данных
# 1. fastq-файл. Для него считываются риды и вероятности
# 2. Набор сгенерированных данных. В набор входят исходные аллели, сгенерированные данные и матрица расстояний
filename = "data/generate1"
if filename.find('fastq') != -1:
    filename_data = filename
    filename_matrix = filename + "_matrix.txt"
else:
    filename_data = filename + "_data.txt"
    filename_alleles = filename + "_alleles.txt"
    filename_matrix = filename + "_matrix.txt"

f = open(filename_data, 'r')
print('Reading file...')
line = f.readline()
i = 1

if filename_data.find('fastq') != -1:
    while line:
        i += 1
        line = f.readline()
        if (i % 2) == 0:
            data.append(line[0:-1])
    f.close()
    print('Extracting data from fastq file...')
    i = 0
    for j in data:
        if (i%2) == 0:
            value.append(j)
        else:
            probability.append(j)
        i += 1
else:
    print("Reading generated data...")
    while line:
        i += 1
        line = f.readline()
        data.append(line[0:-1])
    f.close()

    print('Reading alleles...')
    value.extend(data)
    allele_f = open(filename_alleles, "r")
    allele_line = allele_f.readline()
    while allele_line:
        alleles.append(allele_line)
        allele_line = allele_f.readline()
    allele_f.close()

    print("Reading distance matrix...")
    with open(filename_matrix, 'r') as file:
        dist_matrix_list = file.read().splitlines()
        file.close()
        buf = dist_matrix_list[0].split(',')
        dist_matrix = np.zeros((len(dist_matrix_list), len(buf)))
        for line in range(len(dist_matrix_list)):
            buf = dist_matrix_list[line].split(',')
            for position in range(len(buf)):
                dist_matrix[line, position] = int(buf[position])

print("Convert string to byte...")
byte_value = clust.convert_str_to_int(value)

#print("Start silhouette plots...")
#clust.silhouette(byte_value, number_of_clusters_mas)

print("Start spectral clustering...")
clust.spectral_clustering(dist_matrix, number_of_clusters_mas)

#print("Start DBSCAN...")
#clust.dbscan(dist_matrix, len(value[0]))