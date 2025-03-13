import pandas as pd
import matplotlib.pyplot as plt

# Configurações do gráfico
plt.rcParams['font.size'] = 12
plt.rcParams['figure.figsize'] = (12, 6)

# Ler o CSV
df = pd.read_csv('results.csv')

# Ordenar por número de threads
df = df.sort_values('ENV_VALUE')

# Calcular speedup (baseado no tempo com 1 thread)
baseline_global = df[df['ENV_VALUE'] == 1]['avg_global_time'].values[0]
baseline_update = df[df['ENV_VALUE'] == 1]['avg_time_update'].values[0]

df['speedup_global'] = baseline_global / df['avg_global_time']
df['speedup_update'] = baseline_update / df['avg_time_update']

# Criar gráficos
fig, (ax1, ax2) = plt.subplots(1, 2)

# Gráfico para Global Time
ax1.plot(df['ENV_VALUE'], df['speedup_global'], 
        marker='o', linestyle='--', color='blue', label='Measured')
ax1.plot(df['ENV_VALUE'], df['ENV_VALUE'], 
        linestyle=':', color='red', label='Ideal')
ax1.set_title('Speedup - Global Time')
ax1.set_xlabel('Number of Threads')
ax1.set_ylabel('Speedup')
ax1.set_ylim(1,1.5)
ax1.legend()
ax1.grid(True)

# Gráfico para Update Time
ax2.plot(df['ENV_VALUE'], df['speedup_update'], 
        marker='s', linestyle='--', color='green', label='Measured')
ax2.plot(df['ENV_VALUE'], df['ENV_VALUE'], 
        linestyle=':', color='red', label='Ideal')
ax2.set_title('Speedup - Update Time')
ax2.set_xlabel('Number of Threads')
ax2.set_ylabel('Speedup')
ax2.set_ylim(1,2)
ax2.legend()
ax2.grid(True)

# Ajustar layout e mostrar
plt.tight_layout()
plt.savefig('speedup_analysis.png', dpi=300, bbox_inches='tight')
plt.show()

# Mostrar dados tabulares
print("\nDados utilizados para os gráficos:")
print(df[['ENV_VALUE', 'speedup_global', 'speedup_update']])
