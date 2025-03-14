import matplotlib.pyplot as plt


# Função para ler os valores de um arquivo
def ler_valores(arquivo):
    with open(arquivo, 'r') as f:
        valores = [float(valor) for valor in f.read().split()]
    return valores


# Ler os valores dos dois arquivos
valores_arquivo1 = ler_valores('update-time-troisieme.txt')
valores_arquivo2 = ler_valores('display-time-troisieme.txt')  # Substitua pelo nome do segundo arquivo

# Gerar as iterações (índices)
iteracoes = range(len(valores_arquivo1))
iteracoes2 = range(len(valores_arquivo2))

# Plotar os valores
plt.plot(iteracoes, valores_arquivo1, label='update time')
plt.plot(iteracoes2, valores_arquivo2, label='display time')

# Adicionar título e labels
plt.title('Valores em função da iteração')
plt.xlabel('Iteração')
plt.ylabel('Valores')

# Adicionar legenda
plt.legend()

# Mostrar o gráfico
plt.show()
