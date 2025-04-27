import math

X_star = 0.1
true_value = math.acos(X_star)

# Набор (а)
x_a = [-0.4, -0.1, 0.2, 0.5]
y_a = [math.acos(x) for x in x_a]

# Набор (б)
x_b = [-0.4, 0.0, 0.2, 0.5]
y_b = [math.acos(x) for x in x_b]

# Лагранж
def lagrange(x, y, X_star):
    n = len(x)
    result = 0.0
    for i in range(n):
        term = y[i]
        for j in range(n):
            if (x[i] - x[j] != 0):
                term *= (X_star - x[j]) / (x[i] - x[j])
        result += term
    return result

# Ньютон + разделенная разность
def newton(x, y, X_star):
    n = len(x)

    div_diff = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        div_diff[i][0] = y[i]

    for j in range(1, n):
        for i in range(n - j):
            xi, xij = x[i], x[i + j]
            f1, f0 = div_diff[i + 1][j - 1], div_diff[i][j - 1]
            div_diff[i][j] = (f0 - f1) / (xi - xij)

    result = div_diff[0][0]
    term = 1.0
    for j in range(1, n):
        term *= (X_star - x[j - 1])
        result += div_diff[0][j] * term

    return result


# Вычисления для набора (а)
L3_a = lagrange(x_a, y_a, X_star)
N3_a = newton(x_a, y_a, X_star)
error_L_a = abs(true_value - L3_a)
error_N_a = abs(true_value - N3_a)

# Вычисления для набора (б)
L3_b = lagrange(x_b, y_b, X_star)
N3_b = newton(x_b, y_b, X_star)
error_L_b = abs(true_value - L3_b)
error_N_b = abs(true_value - N3_b)

# Вывод
print(f"\nТочное значение arccos({X_star}) = {true_value:.5f}\n")

print("Набор (а): x = [-0.4, -0.1, 0.2, 0.5]")
print(f"L = {L3_a:.5f}, погрешность = {error_L_a:.5f}")
print(f"N = {N3_a:.5f}, погрешность = {error_N_a:.5f}\n")

print("Набор (б): x = [-0.4, 0.0, 0.2, 0.5]")
print(f"L = {L3_b:.5f}, погрешность = {error_L_b:.5f}")
print(f"N = {N3_b:.5f}, погрешность = {error_N_b:.5f}")
