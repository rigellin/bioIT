import random
import clustering_functions as clust
import statistic_function as stat

# Генерируем 6 случайных ридов
print("Generating alleles...")
source_reads = []
value = []
random.seed()
allele_numbers = random.randint(3, 6)
filename = "data/generate1"
while len(source_reads) < allele_numbers:
    read = []
    while len(read) < 285:
        read_number = random.randint(0, 3)
        if read_number == 0:
            read.append('A')
        if read_number == 1:
            read.append('G')
        if read_number == 2:
            read.append('T')
        if read_number == 3:
            read.append('C')
    source_reads.append(read)
read_count = random.randint(100, 200)
print(read_count)
read_iter = 0
source = []
export_file = open(filename + "_data.txt", "a")
print("Generating reads...")
while read_iter < read_count:
    number_of_errors = random.randint(0, int(len(source_reads[0])/3))
    source.clear()
    random_read_number = random.randint(0, len(source_reads) - 1)
    source.extend(source_reads[random_read_number])
    for i in range(number_of_errors):
        error_nucl = random.randint(0, len(source) - 1)
        new_nucl = random.randint(0, 3)
        if new_nucl == 0:
            source[error_nucl] = 'A'
        if new_nucl == 1:
            source[error_nucl] = 'G'
        if new_nucl == 2:
            source[error_nucl] = 'T'
        if new_nucl == 3:
            source[error_nucl] = 'C'
    buf_read = ''
    for nucl1 in source:
        buf_read += nucl1
    buf_read += '\n'
    export_file.write(buf_read)
    value.append(buf_read)
    read_iter += 1
export_file.close()
print("Write alleles...")
export_alleles = open(filename + "_alleles.txt", "a")
for allele in source_reads:
    buf_allele = ''
    for nucl in allele:
        buf_allele += nucl
    buf_allele += '\n'
    export_alleles.write(buf_allele)
export_alleles.close()
print("Write reads...")

print("Calculating distance matrix...")
dist_matrix = clust.calculate_dist_matrix(value)
file = open(filename + "_matrix.txt", 'a')
for i in dist_matrix:
    string = ''
    for j in range(len(i)):
        if j == 0:
            string += str(int(i[j]))
        else:
            string += "," + str(int(i[j]))
    string += "\n"
    file.write(string)
file.close()

print("Checking matrix...")
stat.check_matrix(dist_matrix, 285)