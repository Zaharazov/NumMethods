import math
import matplotlib.pyplot as plt

# Задаем начальные условия
x0 = 1.0
xn = 2.0
y0 = 2 * math.exp(1)  # y(1)
z0 = 2 * math.exp(1)  # y'(1)
h = 0.1  # Шаг

# Правая часть уравнения y'' = f(x, y, y')
def f(x, y, z):
    return (1 / math.sqrt(x)) * z - (1 / (4 * x ** 2)) * (x + math.sqrt(x) - 8) * y

# Точное решение
def exact_solution(x):
    return (x**2 + 1 / x) * math.exp(math.sqrt(x))

def euler_method(f, x0, xn, y0, z0, h):
    xs = [x0] # x
    ys = [y0] # y
    zs = [z0] # y'

    x = x0
    y = y0
    z = z0

    while x < xn - 1e-10:
        y_new = y + h * z
        z_new = z + h * f(x, y, z)
        x += h

        xs.append(x)
        ys.append(y_new)
        zs.append(z_new)

        y, z = y_new, z_new

    return xs, ys, zs

def runge_kutta_4(f, x0, xn, y0, z0, h):
    xs = [x0]
    ys = [y0]
    zs = [z0]

    x = x0
    y = y0
    z = z0

    while x < xn - 1e-10:
        k1y = z
        k1z = f(x, y, z)

        k2y = z + h * k1z / 2
        k2z = f(x + h/2, y + h * k1y / 2, z + h * k1z / 2)

        k3y = z + h * k2z / 2
        k3z = f(x + h/2, y + h * k2y / 2, z + h * k2z / 2)

        k4y = z + h * k3z
        k4z = f(x + h, y + h * k3y, z + h * k3z)

        y += h / 6 * (k1y + 2*k2y + 2*k3y + k4y)
        z += h / 6 * (k1z + 2*k2z + 2*k3z + k4z)
        x += h

        xs.append(x)
        ys.append(y)
        zs.append(z)

    return xs, ys, zs

def adams_4(f, x0, xn, y0, z0, h):
    # Используем метод Рунге-Кутты 4-го порядка для первых 4 шагов
    xs, ys, zs = runge_kutta_4(f, x0, x0 + 3*h, y0, z0, h)

    # Вычисляем значения правой части (т.е. f_k) для первых 4 узлов
    fs = [f(xs[i], ys[i], zs[i]) for i in range(4)]

    # Начинаем с индекса 3 (уже есть 4 значения)
    k = 3
    while xs[-1] + h <= xn + 1e-12:
        x_k   = xs[k]
        y_k   = ys[k]
        z_k   = zs[k]
        f_k   = fs[k]

        z_k1  = zs[k-1]
        z_k2  = zs[k-2]
        z_k3  = zs[k-3]

        f_k1  = fs[k-1]
        f_k2  = fs[k-2]
        f_k3  = fs[k-3]

        y_next = y_k + h / 24 * (55 * z_k - 59 * z_k1 + 37 * z_k2 - 9 * z_k3)
        z_next = z_k + h / 24 * (55 * f_k - 59 * f_k1 + 37 * f_k2 - 9 * f_k3)
        x_next = x_k + h

        xs.append(x_next)
        ys.append(y_next)
        zs.append(z_next)
        fs.append(f(x_next, y_next, z_next))

        k += 1

    return xs, ys, zs


def runge_romberg(yh, yh2, p):
    return [
        yh2[2 * i] + (yh2[2 * i] - yh[i]) / (2 ** p - 1)
        for i in range(len(yh))
    ]

def mae(y_exact, y_approx):
    errors = [abs(ye - ya) for ye, ya in zip(y_exact, y_approx)]
    return sum(errors) / len(errors)

def plot_results(name, xs, ys, xs_half, ys_half, y_exact, ys_rr):
    plt.figure(figsize=(8, 5))
    plt.plot(xs, y_exact, 'k-', label='Точное решение', linewidth=2)
    plt.plot(xs, ys, 'bo-', label=f'Численное h={h}')
    plt.plot(xs_half, ys_half, 'r.-', label=f'Численное h={h/2}')
    plt.plot(xs, ys_rr, 'g--', label='Рунге–Ромберг')
    plt.title(name)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.show()

def analyze_method(name, method_func, f, x0, xn, y0, z0, h, p):
    # Решение с шагом h
    xs, ys, _ = method_func(f, x0, xn, y0, z0, h)

    # Решение с шагом h/2
    xs_half, ys_half, _ = method_func(f, x0, xn, y0, z0, h / 2)

    # Уточнённое решение по Рунге–Ромбергу
    ys_rr = runge_romberg(ys, ys_half, p)

    # Точное решение
    y_exact = [exact_solution(x) for x in xs]

    # MAE
    mae_orig = mae(y_exact, ys)
    mae_rr = mae(y_exact, ys_rr)

    # Таблица
    print(f"\n{name}")
    print(f"{'x':>6} | {'y числ.':>14} | {'y точное':>14} | {'Погрешность':>12} | {'y Р-Р':>14} | {'Погр. Р-Р':>12}")
    print("-" * 85)
    for x, y, y_rr, y_ex in zip(xs, ys, ys_rr, y_exact):
        err = abs(y - y_ex)
        err_rr = abs(y_rr - y_ex)
        print(f"{x:6.2f} | {y:14.8f} | {y_ex:14.8f} | {err:12.2e} | {y_rr:14.8f} | {err_rr:12.2e}")

    print(f"\nСредняя ошибка (MAE):")
    print(f"  Для метода: {mae_orig:.2e}")
    print(f"  Для Рунге–Ромберга : {mae_rr:.2e}")

    plot_results(name, xs, ys, xs_half, ys_half, y_exact, ys_rr)

analyze_method("Метод Эйлера", euler_method, f, x0, xn, y0, z0, h, p=1)
analyze_method("Метод Рунге–Кутты", runge_kutta_4, f, x0, xn, y0, z0, h, p=4)
analyze_method("Метод Адамса", adams_4, f, x0, xn, y0, z0, h, p=4)
