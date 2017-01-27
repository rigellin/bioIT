import numpy as np
import matplotlib.pyplot as plt

#Чтение из файла в массив
f = open('ash110.__class1_amp1__.R1.fastq', 'r')
line = f.readline()
i = 1
data = []
while line:
    i+=1
    line = f.readline()
    if (i%2)==0:
        data.append(line[0:-1])
f.close()

#Значения и вероятности в отдельные списки
value = []
probability = []
i=0
for j in data:
    if (i%2)==0:
        value.append(j)
    else:
        probability.append(j)
    i+=1

#Список вероятностей
probability_list = ['!','\"','#','$','%','&','\'','(',')','*','+',',','-','.','/','0','1','2','3','4','5','6','7','8','9',':',';','<','=','>','?','@','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','[','\\',']','^','_','`','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','{','|','}','~']

#Формируем массивы значений для графиков
plot_a = [0]*100
plot_c = [0]*100
plot_g = [0]*100
plot_t = [0]*100
plot_median = [0]*len(value)

for current_value,current_probability in zip(value, probability):
    median = 0
    for i,j in zip(current_value,current_probability):
        if i=='A':
            plot_a[probability_list.index(j)]+=1
        if i=='C':
            plot_c[probability_list.index(j)]+=1
        if i=='G':
            plot_g[probability_list.index(j)]+=1
        if i=='T':
            plot_t[probability_list.index(j)]+=1
        median += probability_list.index(j)
    plot_median[value.index(current_value)] = round(median/len(current_value))

#Массив значений аргумента
x=np.arange(0,100,1)

#Строим количественный график
plt.subplot(211)
plt.plot(x,plot_a, label=r"А") 
plt.plot(x,plot_c, label=r"C") 
plt.plot(x,plot_g, label=r"G") 
plt.plot(x,plot_t, label=r"T")
plt.title(r'$График\ качества\ нуклеотидов$')
plt.xlabel(r'$Вероятности$')
plt.ylabel(r'$Количество\ нуклеотидов$')
plt.grid(True) 

#Расстояние Хэмминга
hem = []
for k in range(len(value)-1):
    hem.append(sum(1 for x,y in zip(value[k], value[k+1]) if x!=y))
hem.append(0)
plt.subplot(212)
plt.hist(hem, 50)
plt.title(r'$Расстояние\ Хэмминга$')
plt.xlabel(r'$Значение$')
plt.ylabel(r'$Количество$')
plt.show()
        
