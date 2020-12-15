import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

df = pd.read_csv('performance.csv', delimiter = ',')

unique = df.nr_cores.unique()

mean_time = df.groupby(['nr_cores', 'str_length']).mean()
std_v = df.groupby(['nr_cores', 'str_length']).std()

values = mean_time['nr_cores'] == 4

mean_time =mean_time.drop('Unnamed: 0', axis = 1).reset_index()

sns.pointplot(data=mean_time, x="str_length", y="time", hue="nr_cores", legend='Nr. Cores')
df_raw = df.drop('Unnamed: 0', axis =1)
sns.pointplot(data=mean_time, x="str_length", y="time", hue="nr_cores", err_style="bars", ci=90)
