import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('performance.csv', delimiter = ',')

unique = df.nr_cores.unique()

mean_time = df.groupby(['nr_cores', 'str_length']).mean()
std_v = df.groupby(['nr_cores', 'str_length']).std()


for i in unique:
    values = df['nr_cores'] == i
    subset = df[values]
    plt.plot(subset['str_length'], subset['time'], 'go--', color = 'blue', linewidth=2, markersize=12)
    plt.title(('Benchmarking with {} Cores').format(i))
    plt.xlabel('Stringe size')
    plt.ylabel('time [s]')

    #plt.show()                                                        #print(df[values])

#print(list(mean_time))
#print(mean_time)

values = mean_time['nr_cores'] == 4
#print(values)

#for i in unique:
#    values = mean_time['nr_cores'] == i
#    subset = mean_time[values]

#print()
