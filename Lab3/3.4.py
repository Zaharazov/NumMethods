# Входные данные

X = [-1.0, 0.0, 1.0, 2.0, 3.0]
f = [2.3562, 1.5708, 0.7854, 0.46365, 0.32175]

x_star = 1.0

i = X.index(x_star)

x1, x2, x3 = X[i - 1], X[i], X[i + 1]
f1, f2, f3 = f[i - 1], f[i], f[i + 1]

# Первая производная
first_derivative_left = (f2 - f1) / (x2 - x1)
first_derivative_right = (f3 - f2) / (x2 - x1)

first_derivative = (f3 - f1) / 2*(x2 - x1)

# Вторая производная

second_derivative_full = 2 * (((f3 - f2)/(x3 - x2)) - ((f2 - f1)/(x2 - x1))) / (x3 - x1)

second_derivative = (f3 - 2 * f2 + f1) / ((x3 - x2) ** 2)

# Вывод результатов
print("Первая производная слева в x =", x2, "равна", first_derivative_left)
print("Первая производная справа в x =", x2, "равна", first_derivative_right)
print("Первая производная в x =", x2, "равна", first_derivative)
print("Вторая производная в x =", x2, "равна", second_derivative)
print("Вторая производная (полная формула) в x =", x2, "равна", second_derivative_full)
