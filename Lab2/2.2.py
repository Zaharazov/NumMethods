import math
import matplotlib.pyplot as plt

# Простые итерации
def phi1(x1, x2):
    return math.sqrt(16 - x2 ** 2)

def phi2(x1, x2):
    return math.log(x1 + 4)

def dphi1_dx1(x1, x2):
    return 0

def dphi1_dx2(x1, x2):
    return -x2 / math.sqrt(16 - x2**2)

def dphi2_dx1(x1, x2):
    return 1 / (x1 + 4)

def dphi2_dx2(x1, x2):
    return 0

def simple_iteration_method(x1, x2, epsilon, max_iter=100):

    G_x1_min, G_x1_max = x1 - 0.5, x1 + 0.5
    G_x2_min, G_x2_max = x2 - 0.5, x2 + 0.5

    q1 = abs(dphi1_dx1(x1, x2)) + abs(dphi1_dx2(x1, x2))
    q2 = abs(dphi2_dx1(x1, x2)) + abs(dphi2_dx2(x1, x2))

    q = max(q1, q2)

    if q >= 1:
        print("Корень найти нельзя")
        return

    for k in range(max_iter):
        x1_new = phi1(x1, x2)
        x2_new = phi2(x1, x2)

        if abs(x1_new - x1) < epsilon and abs(x2_new - x2) < epsilon:
            print(f"Корень найден: x1 = {round(x1_new, 7)}, x2 = {round(x2_new, 7)} на итерации {k + 1}")
            return

        if not (G_x1_min <= x1_new <= G_x1_max and G_x2_min <= x2_new <= G_x2_max):
            print("Последовательные приближения вышли за пределы области G")
            return

        x1, x2 = x1_new, x2_new

    print("Не удалось найти корень за максимальное количество итераций")
    return

# Ньютон
def f1(x1, x2):
    return x1 ** 2 + x2 ** 2 - 16

def f2(x1, x2):
    return x1 - math.exp(x2) + 4

def df1_dx1(x1, x2):
    return 2 * x1

def df1_dx2(x1, x2):
    return 2 * x2

def df2_dx1(x1, x2):
    return 1

def df2_dx2(x1, x2):
    return -math.exp(x2)

def newton_method(x1, x2, epsilon, max_iter=100):
    for k in range(max_iter):
        F1 = f1(x1, x2)
        F2 = f2(x1, x2)

        # Матрица Якоби
        J11 = df1_dx1(x1, x2)
        J12 = df1_dx2(x1, x2)
        J21 = df2_dx1(x1, x2)
        J22 = df2_dx2(x1, x2)

        det_J = J11 * J22 - J12 * J21

        if abs(det_J) < 1e-10:
            print("Якобиан вырожден, метод не применим")
            return

        # Обратная матрица
        inv_J11 = J22 / det_J
        inv_J12 = -J12 / det_J
        inv_J21 = -J21 / det_J
        inv_J22 = J11 / det_J

        # f_k + J_k * dx = 0
        # J_k * dx = - f_k
        # dx = -(J^-1_k * f_k)
        dx1 = - (inv_J11 * F1 + inv_J12 * F2)
        dx2 = - (inv_J21 * F1 + inv_J22 * F2)

        # x_k+1 = x_k + dx_k
        x1_new = x1 + dx1
        x2_new = x2 + dx2

        if abs(x1_new - x1) < epsilon and abs(x2_new - x2) < epsilon:
            print(f"Решение найдено: x1 = {round(x1_new, 6)}, x2 = {round(x2_new, 6)} на итерации {k + 1}")
            return x1_new, x2_new

        x1, x2 = x1_new, x2_new

    print("Не удалось достичь требуемой точности за максимальное количество итераций")
    return

# График
def f1_graph(x1):
    return math.sqrt(16 - x1**2)

def f2_graph(x1):
    return -math.sqrt(16 - x1**2)

def f3_graph(x2):
    return math.exp(x2) - 4

x1_values = [i * 0.01 for i in range(-400, 401)]

y1_values = [f1_graph(x1) for x1 in x1_values]
y2_values = [f2_graph(x1) for x1 in x1_values]

x2_values = [i * 0.01 for i in range(-240, 240)]
x1_values_from_g = [f3_graph(x2) for x2 in x2_values]

plt.figure(figsize=(8, 6))

plt.plot(x1_values, y1_values, label=r'$x_1^2 + x_2^2 = 16$', color='b')
plt.plot(x1_values, y2_values, color='b')

plt.plot(x1_values_from_g, x2_values, label=r'$x_1 = e^{x_2} - 4$', color='r')

plt.axhline(0, color='black', linewidth=0.8)
plt.axvline(0, color='black', linewidth=0.8)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xlabel(r'$x_1$', fontsize=12)
plt.ylabel(r'$x_2$', fontsize=12)
plt.legend(loc='best', fontsize=10)
plt.title('Графики системы нелинейных уравнений', fontsize=14)

plt.show()

x0_1 = 3.2
x0_2 = 2
epsilon = 1e-6

newton_method(x0_1, x0_2, epsilon)
simple_iteration_method(x0_1, x0_2, epsilon)
