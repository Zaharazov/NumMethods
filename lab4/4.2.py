import math
import matplotlib.pyplot as plt

# Правая часть уравнения y'' = f(x, y, y')
def f(x, y, dy):
    return ((2 * x + 1) / x) * dy - ((x + 1) / x) * y

# Точное решение
def exact_solution(x):
    return math.exp(x) * (x ** 2 + 1)

# Интервал и шаг
x0 = 1e-5
x1 = 1.0
h = 0.05
y0_prime = 1  # y'(x0)
alpha0 = 0.0
alpha1 = 20.0

# Краевые условия
# y'(x0) = 1
# y'(x1) - 2*y(x1) = 0

def tridiagonal_solve(a, b, c, d):
    n = len(b)

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

def shooting_method(f, x0, x1, y0_prime, h, alpha0, alpha1, eps=1e-6, max_iter=50):
    def phi(alpha):
        # Решаем задачу Коши при данной догадке alpha = y(x0)
        _, ys, zs = runge_kutta_4(f, x0, x1, alpha, y0_prime, h)
        y1 = ys[-1]
        z1 = zs[-1]
        return z1 - 2 * y1  # это граница: y'(1) - 2y(1)

    a0 = alpha0
    a1 = alpha1
    phi0 = phi(a0)
    phi1 = phi(a1)

    for _ in range(max_iter):
        if abs(phi1) < eps:
            break
        # Метод секущих
        a2 = a1 - phi1 * (a1 - a0) / (phi1 - phi0)
        a0, phi0 = a1, phi1
        a1, phi1 = a2, phi(a2)

    # Финальное решение с найденным y(x0)
    xs, ys, zs = runge_kutta_4(f, x0, x1, a1, y0_prime, h)
    return xs, ys, zs

def finite_difference_method(x0, x1, h, y0_prime):
    # Ручная генерация xs — гарантированно доходит до x1
    xs = []
    x = x0
    while x < x1:
        xs.append(round(x, 10))  # защита от накопления ошибки
        x += h
    xs.append(x1)
    N = len(xs) - 1

    a_full = [0.0] * (N + 1)  # поддиагональ (a[1]..a[N])
    b = [0.0] * (N + 1)       # главная диагональ
    c_full = [0.0] * (N + 1)  # наддиагональ (c[0]..c[N-1])
    d = [0.0] * (N + 1)       # правая часть

    # Левая граница: y'(x0) = y0_prime ≈ (y1 - y0) / h
    b[0] = -1 / h
    c_full[0] = 1 / h
    d[0] = y0_prime

    # Внутренние узлы
    for i in range(1, N):
        xi = xs[i]
        pi = -(2 * xi + 1) / xi
        qi = (xi + 1) / xi

        a_full[i] = 1 / h**2 - pi / (2 * h)
        b[i] = -2 / h**2 + qi
        c_full[i] = 1 / h**2 + pi / (2 * h)
        d[i] = 0.0

    # Правая граница: y'(x1) - 2*y(x1) = 0
    a_full[N] = -1 / h
    b[N] = 1 / h - 2
    d[N] = 0.0

    # Подготовим срезы поддиагонали и наддиагонали
    a = a_full[1:]     # от a[1] до a[N] → длина N
    c = c_full[:-1]    # от c[0] до c[N-1] → длина N

    y = tridiagonal_solve(a, b, c, d)
    return xs, y

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

def analyze_method(name, method_func, f, x0, x1, y0, z0, h, p):
    # Решение с шагом h
    xs, ys = method_func(f, x0, x1, y0, z0, h)

    # Решение с шагом h/2
    xs_half, ys_half = method_func(f, x0, x1, y0, z0, h / 2)

    # Уточнение по Рунге–Ромбергу
    ys_rr = runge_romberg(ys, ys_half, p)

    y_exact = [exact_solution(x) for x in xs]

    mae_orig = mae(y_exact, ys)
    mae_rr = mae(y_exact, ys_rr)

    print(f"\n{name}")
    print(f"{'x':>6} | {'y числ.':>14} | {'y точное':>14} | {'Погрешность':>12} | {'y Р-Р':>14} | {'Погр. Р-Р':>12}")
    print("-" * 85)
    for x, y, y_rr, y_ex in zip(xs, ys, ys_rr, y_exact):
        err = abs(y - y_ex)
        err_rr = abs(y_rr - y_ex)
        print(f"{x:6.2f} | {y:14.8f} | {y_ex:14.8f} | {err:12.2e} | {y_rr:14.8f} | {err_rr:12.2e}")

    print(f"\nСредняя ошибка (MAE):")
    print(f"  Без Рунге–Ромберга: {mae_orig:.2e}")
    print(f"  С Рунге–Ромбергом : {mae_rr:.2e}")

    plot_results(name, xs, ys, xs_half, ys_half, y_exact, ys_rr)

# Обёртка под метод стрельбы
def shooting_wrapper(f, x0, x1, y0, y0_prime, h):
    xs, ys, _ = shooting_method(f, x0, x1, y0_prime, h, alpha0, alpha1)
    return xs, ys

# Обёртка под КРМ
def finite_diff_wrapper(f, x0, x1, y0, y0_prime, h):
    xs, ys = finite_difference_method(x0, x1, h, y0_prime)
    return xs, ys

analyze_method("Метод стрельбы", shooting_wrapper, f, x0, x1, None, y0_prime, h, p=4)
analyze_method("Конечно-разностный метод", finite_diff_wrapper, f, x0, x1, None, y0_prime, h, p=2)
