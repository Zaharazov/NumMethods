def tridiagonal_solve(a, b, c, d):
    n = len(b)

    is_condition_met = False # проверка сходимости

    for i in range(n):
        if i == 0:
            if abs(b[i]) <= abs(c[i]):
                return
            if abs(b[i]) > abs(c[i]):
                is_condition_met = True

        elif i == n - 1:
            if abs(b[i]) <= abs(a[i - 1]):
                return
            if abs(b[i]) > abs(a[i - 1]):
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
        alpha[i] = c[i] / (b[i] - a[i - 1] * alpha[i - 1])
        beta[i] = (d[i] - a[i - 1] * beta[i - 1]) / (b[i] - a[i - 1] * alpha[i - 1])

    beta[n - 1] = (d[n - 1] - a[n - 2] * beta[n - 2]) / (b[n - 1] - a[n - 2] * alpha[n - 2])

    # обратный ход
    x = [0] * n
    x[n - 1] = beta[n - 1]

    for i in range(n - 2, -1, -1):
        x[i] = beta[i] - alpha[i] * x[i + 1]

    x = [round(val, 6) for val in x]

    return x

# Размер системы
n = 5

# Коэффициенты матрицы
a = [-2, 2, -8, -7]
b = [8, 12, -9, 17, 13]
c = [-4, -7, 1, -4]
d = [32, 15, -10, 133, -76]

x = tridiagonal_solve(a, b, c, d)

print("Решение системы:", x)
