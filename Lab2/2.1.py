import matplotlib.pyplot as plt
import math


def phi(x):
    expr = -x ** 2 + 2 * x + 1
    if expr < 0:
        return float('nan')
    return expr ** (1 / 3)

def simple_iteration_method(a, b, epsilon, max_iter=100):
    x = b

    for val in [a + i * 0.01 for i in range(int(abs(b - a) / 0.01) + 1)]:
        phi_val = phi(val)
        if math.isnan(phi_val) or not (a <= phi_val <= b):
            print(f"Функция φ(x) выходит за пределы.")
            return

    for k in range(max_iter):
        x_new = phi(x)

        if abs(x_new - x) < epsilon:
            print(f"Корень найден: {round(x_new, 7)} на итерации {k + 1}")
            return
        x = x_new

    print("Не удалось найти корень за максимальное количество итераций")
    return


def f(x):
    return x**3 + x**2 - 2*x - 1

def df(x):
    return 3*x**2 + 2*x - 2

def ddf(x):
    return 6*x + 2

def newton_method(a, b, epsilon, max_iter=100):

    x = (a+b)/2

    if f(a) * f(b) >= 0:
        print("Начальное приближение не подходит")
        quit()

    check_a = False
    check_b = False
    if abs(f(a) * ddf(a)) < df(a) * df(a):
        check_a = True
    if abs(f(b) * ddf(b)) < df(b) * df(b):
        check_b = True

    if check_b and check_a:
        if f(b) * ddf(b) <= 0:
            print("Начальное приближение не подходит")
            quit()

    if check_a:
        x = a
    if check_b:
        x = b

    for k in range(max_iter):
        x_new = x - f(x) / df(x)

        if abs(x_new - x) < epsilon:
            print(f"Корень найден: {round(x_new, 7)} на итерации {k + 1}")
            return
        x = x_new

    print("Не удалось найти корень за максимальное количество итераций")
    return


def f1(x):
    return x**3 + x**2

def f2(x):
    return 2*x + 1

x = [i * 0.01 for i in range(-300, 300)]
y1 = [f1(i) for i in x]
y2 = [f2(i) for i in x]

plt.figure(figsize=(8, 6))
plt.plot(x, y1, label='$f_1(x) = x^3 + x^2$')
plt.plot(x, y2, label='$f_2(x) = 2x + 1$')
plt.axhline(0, color='black', linewidth=0.8)
plt.axvline(0, color='black', linewidth=0.8)
plt.grid()
plt.legend()
plt.show()

x0_a = 0
x0_b = 2
epsilon = 1e-6

#x0 = float(input())

newton_method(x0_a, x0_b, epsilon)
simple_iteration_method(x0_a, x0_b, epsilon)
