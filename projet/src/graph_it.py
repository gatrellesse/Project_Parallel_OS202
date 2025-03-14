import matplotlib.pyplot as plt


# Função para ler os valores de um arquivo
def ler_valores(arquivo):
    with open(arquivo, 'r') as f:
        valores = [float(valor) for valor in f.read().split()]
    return valores


# Ler os valores dos dois arquivos
valores_arquivo1 = ler_valores('update-time-1.txt')
valores_arquivo2 = ler_valores('display-time-1.txt')  # Substitua pelo nome do segundo arquivo

# Gerar as iterações (índices)
iteracoes = range(len(valores_arquivo1))
iteracoes2 = range(len(valores_arquivo2))

# Plotar os valores
plt.plot(iteracoes, valores_arquivo1, label='update time')
plt.plot(iteracoes2, valores_arquivo2, label='display time')

# Adicionar título e labels
plt.title('Loop time per iteration')
plt.xlabel('Iteration')
plt.ylabel('Time (s)')
# Adicionar legenda
plt.legend()

# Mostrar o gráfico
plt.savefig('time_analysis.png', dpi=300)
plt.show()
