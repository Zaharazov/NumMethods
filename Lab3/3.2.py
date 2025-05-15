# Метод прогонки
def tridiagonal_solve(a, b, c, d):
    n = len(b)

    is_condition_met = False

    for i in range(n):
        if i == 0:
            if abs(b[i]) <= abs(c[i]):
                return
            if abs(b[i]) > abs(c[i]):
                is_condition_met = True

        elif i == n - 1:
            if abs(b[i]) <= abs(a[i - 1]):
                return
            if abs(b[i]) > abs(a[i - 1]): # добавит в отчет преобразование формул
                is_condition_met = True

        else:
            if abs(b[i]) < abs(a[i - 1]) + abs(c[i]):
                return
            if abs(b[i]) > abs(a[i - 1]) + abs(c[i]):
                is_condition_met = True

    if not is_condition_met:
        return

    # прямой ход
    alpha = [0] * n
    beta = [0] * n

    alpha[0] = c[0] / b[0]
    beta[0] = d[0] / b[0]

    for i in range(1, n - 1):
        denom = b[i] - a[i - 1] * alpha[i - 1]
        alpha[i] = c[i] / denom
        beta[i] = (d[i] - a[i - 1] * beta[i - 1]) / denom

    beta[n - 1] = (d[n - 1] - a[n - 2] * beta[n - 2]) / (b[n - 1] - a[n - 2] * alpha[n - 2])

    # обратный ход
    x = [0] * n
    x[n - 1] = beta[n - 1]

    for i in range(n - 2, -1, -1):
        x[i] = beta[i] - alpha[i] * x[i + 1]

    x = [round(val, 6) for val in x]

    return x

# Входные данные
x_star = 0.1

xi = [-0.4, -0.1, 0.2, 0.5, 0.8]
fi = [1.9823, 1.6710, 1.3694, 1.0472, 0.64350]
n = len(xi) - 1

# h[i]
h = [xi[i] - xi[i-1] for i in range(1, n+1)]

# Матрица системы
A = [0 for _ in range(n-1)]
B = [0 for _ in range(n-1)]
C = [0 for _ in range(n-1)]
D = [0 for _ in range(n-1)]

for i in range(1, n):
    hi_1 = h[i-1]
    hi = h[i]
    A[i-1] = hi_1
    B[i-1] = 2 * (hi_1 + hi)
    C[i-1] = hi
    D[i-1] = 3 * ((fi[i+1] - fi[i]) / hi - (fi[i] - fi[i-1]) / hi_1)

# Решение системы
c_internal = tridiagonal_solve(A[1:], B, C[:-1], D)

c = [0] + c_internal + [0]

# Коэффициенты a, b, d
a = [0] * n
b = [0] * n
d = [0] * n

for i in range(n):
    a[i] = fi[i]
    b[i] = ((fi[i+1] - fi[i]) / h[i]) - (h[i] / 3) * (2 * c[i] + c[i+1])
    d[i] = (c[i+1] - c[i]) / (3 * h[i])

# Вычисление значения сплайна в x*
def evaluate_spline(x_star):
    for i in range(n):
        if xi[i] <= x_star <= xi[i+1]:
            dx = x_star - xi[i]
            return a[i] + b[i]*dx + c[i]*dx**2 + d[i]*dx**3
    raise ValueError("x* вне диапазона интерполяции")

# Пример вызова
y_star = evaluate_spline(x_star)
print(f"S({x_star}) = {y_star:.6f}")


import matplotlib.pyplot as plt
x_dense = []
step = 0.005
x_val = xi[0]
while x_val <= xi[-1]:
    x_dense.append(round(x_val, 5))
    x_val += step

y_dense = [evaluate_spline(x) for x in x_dense]

# График
plt.plot(xi, fi, 'ro', label='Узлы интерполяции')
plt.plot(x_dense, y_dense, 'b-', label='Кубический сплайн')
plt.title("Интерполяция кубическим сплайном")
plt.xlabel("x")
plt.ylabel("S(x)")
plt.legend()
plt.grid(True)
plt.show()
