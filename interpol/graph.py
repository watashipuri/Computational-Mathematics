import pandas as pd
import matplotlib.pyplot as plt

# Загрузка данных из CSV
data = pd.read_csv("results.csv")

# Создаем словарь для перевода длины отрезка в формат дробей
length_labels = {
    2.0: "2",
    1.0: "1",
    0.5: "1/2",
    0.25: "1/4",
    0.125: "1/8",
    0.0625: "1/16",
    0.03125: "1/32",
    0.015625: "1/64"
}

# Сортируем данные по убыванию длины отрезка
data['L'] = pd.Categorical(data['L'], categories=sorted(length_labels.keys(), reverse=True), ordered=True)
data = data.sort_values('L', ascending=False)

# Построение графиков для каждого значения N
for N in sorted(data['N'].unique()):
    subset = data[data['N'] == N]
    plt.plot(subset['L'].map(length_labels), subset['Max Error'], marker='o', label=f'N={N}')

plt.xlabel('Длина отрезка L')
plt.ylabel('Максимальная ошибка интерполяции')
plt.yscale('log')  # Устанавливаем логарифмическую шкалу для оси Y
plt.title('Зависимость ошибки интерполяции от длины отрезка')
plt.gca().invert_xaxis()  # Инвертируем ось X, чтобы длина шла на уменьшение
plt.legend()
plt.grid(True, which="both", linestyle='--')  # Добавляем сетку для логарифмической шкалы
plt.show()
