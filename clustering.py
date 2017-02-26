import numpy as np
import matplotlib.pyplot as plt
import clustering_functions as clust

# Объявление переменных
# Массив данных, считанных из файла. В нём содержатся и значения нуклеотидов, и вероятности
data = []
# Массив данных, полученных из data[]. В нём содержатся значения нуклеотидов
value = []
# Массив данных, полученных из data[]. В нём содержатся значения вероятностей
probability = []
# Количество кластеров
number_of_clusters = 10
# Центры кластеров
cluster_centers = []
# Список для кластеров
clusters = []
# Список вероятностей
probability_list = ['!','\"','#','$','%','&','\'','(',')','*','+',',','-','.','/','0','1','2','3','4','5','6','7','8','9',':',';','<','=','>','?','@','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','[','\\',']','^','_','`','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','{','|','}','~']

# Чтение из файла в массив. Считываем только чётные строки.
# В результате получаем полный массив данных - значение нуклеотида и вероятность
f = open('ash110.__class1_amp1__.R1.fastq', 'r')
# Считываем первую строку
print('Reading file...')
line = f.readline()
# Итератор для подсчёта строк
i = 1
# Цикл, считывающий каждую чётную строку в массив
while line:
    i += 1
    line = f.readline()
    if (i % 2) == 0:
        data.append(line[0:-1])
f.close()

# Значения и вероятности разносим в отдельные списки из data[]
print('Extracting data...')
i = 0
for j in data:
    if (i%2) == 0:
        value.append(j)
    else:
        probability.append(j)
    i += 1

# Получаем первые значения центров кластеров
cluster_centers.extend(clust.calculate_first_centers(value, len(value), number_of_clusters))

# Первые значения получены, начинаем кластерный анализ
clusters.extend(clust.clustering_data(value, cluster_centers))
