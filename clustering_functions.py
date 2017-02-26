import numpy as np

# Вычисление расстояния Хэмминга
def hemming(str1, str2):
    return (sum(1 for x,y in zip(str1, str2) if x!=y))

# Функция для составления статистики по массиву ридов - вычисляет, какой элемент
# чаще всего встречается на данной позиции во входном наборе данных
def statistic_read(input_data, read_length):
    print('Get statistic information...')
    # Выходное значение - строка, состоящая из средних значений
    result_read = ''
    # Существует 4 нуклеотида - A, T, G, C
    a_count = [0] * read_length
    t_count = [0] * read_length
    g_count = [0] * read_length
    c_count = [0] * read_length
    for i in input_data:
        # Позиция нуклеотида в риде
        nucl_index = 0
        # Подсчитываем количество каждого возможного значения нуклеотида для
        # каждой позиции в риде
        for j in i:
            if j == 'A':
                a_count[nucl_index] += 1
            if j == 'G':
                g_count[nucl_index] += 1
            if j == 'T':
                t_count[nucl_index] += 1
            if j == 'C':
                c_count[nucl_index] += 1
            nucl_index += 1
    # Для каждой позиции выбираем элемент, который встречается чаще остальных
    for i in range(0, read_length):
        # Количество вхождений нуклеотида для определённой позиции
        nucl_count = [a_count[i], g_count[i], t_count[i], c_count[i]]
        if nucl_count.index(max(nucl_count)) == 0:
            result_read += 'A'
        if nucl_count.index(max(nucl_count)) == 1:
            result_read += 'G'
        if nucl_count.index(max(nucl_count)) == 2:
            result_read += 'T'
        if nucl_count.index(max(nucl_count)) == 3:
            result_read += 'C'
    return result_read


# Вычисление первых значений центра кластеров
def calculate_first_centers(value, read_count, number_of_clusters):
    print('Calculating first centers...')
    centers = []
    # Вычисляем примерное количество ридов в каждом кластере
    number_of_read = read_count // number_of_clusters
    # Создаём массив "границ" кластера
    cluster_border = []
    # Создаём массив значений предполагаемого кластера
    cluster_value = []
    cluster_border.append(0)
    # Определяем границы кластеров
    while read_count > (cluster_border[-1] + number_of_read):
        cluster_border.append(cluster_border[-1] + number_of_read)
    # Добавляем последнюю границу
    cluster_border[-1] = read_count - 1
    # Вычисляем центры кластеров - для каждой позиции рида выбираем самое часто встречаемое значение в пределе области
    number_of_clusters_iter = 0
    while number_of_clusters_iter < number_of_clusters:
        for j in range(cluster_border[number_of_clusters_iter], cluster_border[number_of_clusters_iter + 1]):
            cluster_value.append(value[j])
        # Получаем среднее значение рида
        middle_read = statistic_read(cluster_value, len(cluster_value[0]))
        centers.append(middle_read)
        number_of_clusters_iter += 1
    return centers

# Для каждой строки из value вычисляем расстояние Хэмминга и относим к определённому кластеру
def clustering_data(value, cluster_centers):
    print('Calculating cluster numbers...')
    # Значение кластера для каждого рида
    clusters = [0] * len(value)
    # Счётчик номера рида
    read_num = 0
    for i in value:
        # Счётчик номера кластера
        cluster_num = 0

        for j in cluster_centers:
            if cluster_num == 0:
                # Счётчик для значения кластера для рида. По умолчанию относим рид к 1 кластеру
                cluster_value = 0
                hem = hemming(i, cluster_centers[0])
            else:
                new_hem = hemming(i, j)
                if new_hem < hem:
                    hem = new_hem
                    cluster_value = cluster_num
            cluster_num += 1
        clusters[read_num] = cluster_value
        read_num += 1
    return clusters