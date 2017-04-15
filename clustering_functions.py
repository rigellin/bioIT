from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import statistic_function as stat
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.cluster import DBSCAN
from sklearn.cluster import SpectralClustering

# Вычисление расстояния Хэмминга
def hemming(str1, str2):
    return (sum(1 for x,y in zip(str1, str2) if x!=y))

def calculate_dist_matrix(value):
    dist_matrix = np.zeros((len(value), len(value)))
    for i in range(len(value)):
        for j in range(len(value)):
            dist_matrix[i, j] = hemming(value[i], value[j])
            if dist_matrix[i, j] == 0:
                dist_matrix[i, j] = 1
    return dist_matrix

# Функция перевода строки из буквенного в числовой формат
def convert_str_to_int(value):
    current_value = np.zeros((len(value[0]) * 2))
    new_value = np.zeros((len(value), len(current_value)))
    array_index = 0
    for read in value:
        nucl_index = 0
        for nucl in read:
            if nucl == 'A':
                current_value[nucl_index] = 0
                current_value[nucl_index + 1] = 0
            if nucl == 'G':
                current_value[nucl_index] = 0
                current_value[nucl_index + 1] = 1
            if nucl == 'T':
                current_value[nucl_index] = 1
                current_value[nucl_index + 1] = 0
            if nucl == 'C':
                current_value[nucl_index] = 1
                current_value[nucl_index + 1] = 1
            nucl_index += 2
        new_value[array_index] = current_value
        array_index += 1
    return new_value

def silhouette(new_value, range_n_clusters):
    coeff_list = np.zeros(len(range_n_clusters))
    iteration = 0
    for n_clusters in range_n_clusters:
        # Create a subplot with 1 row and 2 columns
        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.set_size_inches(18, 7)

        # The 1st subplot is the silhouette plot
        # The silhouette coefficient can range from -1, 1 but in this example all
        # lie within [-0.1, 1]
        ax1.set_xlim([-0.1, 1])
        # The (n_clusters+1)*10 is for inserting blank space between silhouette
        # plots of individual clusters, to demarcate them clearly.
        ax1.set_ylim([0, len(new_value) + (n_clusters + 1) * 10])

        # Initialize the clusterer with n_clusters value and a random generator
        # seed of 10 for reproducibility.
        clusterer = KMeans(n_clusters=n_clusters, random_state=10)
        cluster_labels = clusterer.fit_predict(new_value)

        # The silhouette_score gives the average value for all the samples.
        # This gives a perspective into the density and separation of the formed
        # clusters
        silhouette_avg = silhouette_score(new_value, cluster_labels)
        coeff_list[iteration] = silhouette_avg
        print("For n_clusters =", n_clusters,
              "The average silhouette_score is :", silhouette_avg)

        # Compute the silhouette scores for each sample
        sample_silhouette_values = silhouette_samples(new_value, cluster_labels)

        y_lower = 10
        for i in range(n_clusters):
            # Aggregate the silhouette scores for samples belonging to
            # cluster i, and sort them
            ith_cluster_silhouette_values = \
                sample_silhouette_values[cluster_labels == i]

            ith_cluster_silhouette_values.sort()

            size_cluster_i = ith_cluster_silhouette_values.shape[0]
            y_upper = y_lower + size_cluster_i

            color = cm.spectral(float(i) / n_clusters)
            ax1.fill_betweenx(np.arange(y_lower, y_upper),
                              0, ith_cluster_silhouette_values,
                              facecolor=color, edgecolor=color, alpha=0.7)

            # Label the silhouette plots with their cluster numbers at the middle
            ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

            # Compute the new y_lower for next plot
            y_lower = y_upper + 10  # 10 for the 0 samples

        ax1.set_title("The silhouette plot for the various clusters.")
        ax1.set_xlabel("The silhouette coefficient values")
        ax1.set_ylabel("Cluster label")

        # The vertical line for average silhouette score of all the values
        ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

        ax1.set_yticks([])  # Clear the yaxis labels / ticks
        ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

        # 2nd Plot showing the actual clusters formed
        colors = cm.spectral(cluster_labels.astype(float) / n_clusters)
        ax2.scatter(new_value[:, 0], new_value[:, 1], marker='.', s=30, lw=0, alpha=0.7,
                    c=colors)

        # Labeling the clusters
        centers = clusterer.cluster_centers_
        # Draw white circles at cluster centers
        ax2.scatter(centers[:, 0], centers[:, 1],
                    marker='o', c="white", alpha=1, s=200)

        for i, c in enumerate(centers):
            ax2.scatter(c[0], c[1], marker='$%d$' % i, alpha=1, s=50)

        ax2.set_title("The visualization of the clustered data.")
        ax2.set_xlabel("Feature space for the 1st feature")
        ax2.set_ylabel("Feature space for the 2nd feature")

        plt.suptitle(("Silhouette analysis for KMeans clustering on sample data "
                      "with n_clusters = %d" % n_clusters),
                     fontsize=14, fontweight='bold')

        plt.show()
        iteration += 1
    max_coeff = np.argmax(coeff_list)
    print("Estimated number of clusters: ", range_n_clusters[max_coeff])

def dbscan(X, read_length):
    peak_dists = stat.find_peak(X, read_length)
    Eps = peak_dists[0] + 5
    stop = 'false'
    last_black_count = 0
    while stop == 'false':
        db = DBSCAN(eps=Eps, min_samples=2, metric='precomputed').fit(X)
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        labels = db.labels_
        # Первый критерий остановки - количество "чёрных" элементов равно 0
        black_count = 0
        for label in labels:
            if label == -1:
                black_count += 1
        if black_count == 0:
            stop = 'true'
        else:
            # Второй критерий остановки - количество "чёрных" элементов не изменилось
            if black_count == last_black_count:
                stop = 'true'
            else:
                # Третий критерий остановки - слишком близко подошли ко второму пику
                if Eps > (peak_dists[1] - 10):
                    stop = 'true'
                else:
                    # Продолжаем искать
                    Eps += 5
                    last_black_count = black_count

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    print('Estimated number of clusters: %d' % n_clusters_)
    # print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
    # print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
    # print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
    # print("Adjusted Rand Index: %0.3f"
    #      % metrics.adjusted_rand_score(labels_true, labels))
    # print("Adjusted Mutual Information: %0.3f"
    #      % metrics.adjusted_mutual_info_score(labels_true, labels))
    #print("Silhouette Coefficient: %0.3f"
          #% metrics.silhouette_score(X, labels))

    ##############################################################################
    # Plot result
    import matplotlib.pyplot as plt

    # Black removed and is used for noise instead.
    unique_labels = set(labels)
    colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
    for k, col in zip(unique_labels, colors):
        if k == -1:
            # Black used for noise.
            col = 'k'

        class_member_mask = (labels == k)

        xy = X[class_member_mask & core_samples_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
                 markeredgecolor='k', markersize=14)

        xy = X[class_member_mask & ~core_samples_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
                 markeredgecolor='k', markersize=6)

    plt.title('Estimated number of clusters: %d' % n_clusters_)
    plt.show()

def spectral_clustering(dist_matrix, number_of_clusters_mas):
    for cluster_count in number_of_clusters_mas:
        spectral = SpectralClustering(n_clusters=cluster_count, eigen_solver='arpack',
                                              affinity="nearest_neighbors").fit(dist_matrix)

        y_pred = spectral.labels_.astype(np.int)
        stat.dan_index(dist_matrix, y_pred)
        colors = np.array([x for x in 'bgrcmykbgrcmykbgrcmykbgrcmyk'])
        colors = np.hstack([colors] * 20)
        plt.scatter(dist_matrix[:, 0], dist_matrix[:, 1], color=colors[y_pred].tolist(), s=10)
        plt.show()