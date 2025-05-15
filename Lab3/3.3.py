import matplotlib.pyplot as plt

# Входные данные
x_vals = [-0.7, -0.4, -0.1, 0.2, 0.5, 0.8]
y_vals = [2.3462, 1.9823, 1.671, 1.3694, 1.0472, 0.6435]

n = len(x_vals)

# LU-азложение из первой лабы
def solve_system(A, b):
    n = len(A)
    A = [row[:] for row in A]

    # LU-разложение
    for k in range(n):
        for i in range(k + 1, n):
            A[i][k] /= A[k][k]
            for j in range(k + 1, n):
                A[i][j] -= A[i][k] * A[k][j]

    # Прямая подстановка (Ly = b)
    y = [0] * n
    for i in range(n):
        y[i] = b[i] - sum(A[i][j] * y[j] for j in range(i))

    # Обратная подстановка (Ux = y)
    x = [0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(A[i][j] * x[j] for j in range(i + 1, n))) / A[i][i]

    return x

def normal_system1(x, y):
    n = len(x)
    sum_x = sum(x)
    sum_x2 = sum(xi**2 for xi in x)
    sum_y = sum(y)
    sum_xy = sum(x[i]*y[i] for i in range(n))

    A = [
        [n,      sum_x],
        [sum_x,  sum_x2]
    ]
    b = [sum_y, sum_xy]
    return A, b

def normal_system2(x, y):
    n = len(x)
    sum_x = sum(x)
    sum_x2 = sum(xi**2 for xi in x)
    sum_x3 = sum(xi**3 for xi in x)
    sum_x4 = sum(xi**4 for xi in x)
    sum_y = sum(y)
    sum_xy = sum(x[i]*y[i] for i in range(n))
    sum_x2y = sum((x[i]**2)*y[i] for i in range(n))

    A = [
        [n,        sum_x,   sum_x2],
        [sum_x,    sum_x2,  sum_x3],
        [sum_x2,   sum_x3,  sum_x4]
    ]
    b = [sum_y, sum_xy, sum_x2y]
    return A, b

def compute_error(x, y, coeffs):
    error = 0.0
    for i in range(len(x)):
        fx = sum(coeffs[j] * x[i]**j for j in range(len(coeffs)))
        error += (fx - y[i]) ** 2
    return error

A1, b1 = normal_system1(x_vals, y_vals)
A2, b2 = normal_system2(x_vals, y_vals)

coeffs1 = solve_system(A1, b1)
coeffs2 = solve_system(A2, b2)

print("Многочлен 1-й степени:")
print(f"F(x) ≈ {coeffs1[0]:.4f} + ({coeffs1[1]:.4f})*x")

print("Многочлен 2-й степени:")
print(f"F(x) ≈ {coeffs2[0]:.4f} + ({coeffs2[1]:.4f})*x + ({coeffs2[2]:.4f})*x^2")

# Вычисление ошибок
error1 = compute_error(x_vals, y_vals, coeffs1)
error2 = compute_error(x_vals, y_vals, coeffs2)

print(f"Сумма квадратов ошибок (1-я степень): {error1:.6f}")
print(f"Сумма квадратов ошибок (2-я степень): {error2:.6f}")

step = 0.01
x_plot = [min(x_vals)-0.1 + i*step for i in range(int((max(x_vals)-min(x_vals)+0.2)/step))]

y_deg1 = [sum(coeffs1[j] * x**j for j in range(len(coeffs1))) for x in x_plot]
y_deg2 = [sum(coeffs2[j] * x**j for j in range(len(coeffs2))) for x in x_plot]

plt.scatter(x_vals, y_vals, color='black', label='Табличные точки')
plt.plot(x_plot, y_deg1, label='Приближение 1-й степени', color='blue')
plt.plot(x_plot, y_deg2, label='Приближение 2-й степени', color='red', linestyle='dashed')
plt.title("МНК: Аппроксимация табличной функции")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True)
plt.legend()
plt.show()
