import numpy as np
import matplotlib.pyplot as plt

# Предобработка

def analytic_solution(x, t):
    return (1 / np.pi**2) * (1 - np.exp(-np.pi**2 * t)) * np.sin(np.pi * x)

def initial_condition(x):
    return np.zeros_like(x)

def source_function(x):
    return np.sin(np.pi * x)

left_boundary = 0.0
right_boundary = 0.0

N = 50
L = 1.0
a = 1.0
h = L / N
tau = 0.0001
T = -np.log(0.05) / (np.pi**2)
K = int(np.ceil(T/tau))
sigma = tau / h**2
plot_times = [0, T/10, T/8, T/6, T/4, T/3, T/2, T]

# Солвер

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


def explicit_scheme():
    x = np.linspace(0, L, N + 1)
    u = initial_condition(x)
    f = source_function(x)

    results = {}
    time_points = [int(round(t / tau)) for t in plot_times]

    for k in range(K + 1):
        t = k * tau
        if k in time_points:
            results[round(t, 6)] = u.copy()

        if k == K:
            break

        u_new = u.copy()
        for i in range(1, N):
            u_new[i] = (sigma * u[i + 1] +
                        (1 - 2 * sigma) * u[i] +
                        sigma * u[i - 1] +
                        tau * f[i])

        u_new[0] = left_boundary
        u_new[N] = right_boundary

        u = u_new

    return x, results


def implicit_scheme():
    x = np.linspace(0, L, N + 1)
    u = initial_condition(x)
    f = source_function(x)

    results = {}
    time_points = [int(round(t / tau)) for t in plot_times]

    m = N - 1

    for k in range(K + 1):
        t = k * tau
        if k in time_points:
            results[round(t, 6)] = u.copy()

        if k == K:
            break

        a_sub = sigma * np.ones(m - 1)
        b_diag = -(1 + 2 * sigma) * np.ones(m)
        c_sup = sigma * np.ones(m - 1)

        d = np.zeros(m)
        for i in range(1, N):
            j = i - 1
            d[j] = -(u[i] + tau * f[i])

        d[0] -= sigma * left_boundary
        d[m - 1] -= sigma * right_boundary

        u_inner = tridiagonal_solve(a_sub, b_diag, c_sup, d)

        u_new = np.zeros(N + 1)
        u_new[0] = left_boundary
        u_new[N] = right_boundary
        u_new[1:N] = u_inner

        u = u_new

    return x, results


def crank_nicolson_scheme():
    x = np.linspace(0, L, N + 1)
    u = initial_condition(x)
    f = source_function(x)

    results = {}
    time_points = [int(round(t / tau)) for t in plot_times]

    m = N - 1
    theta = 0.5

    for k in range(K + 1):
        t = k * tau
        if k in time_points:
            results[round(t, 6)] = u.copy()

        if k == K:
            break

        a_sub = -theta * sigma * np.ones(m - 1)
        b_diag = (1 + 2 * theta * sigma) * np.ones(m)
        c_sup = -theta * sigma * np.ones(m - 1)

        d = np.zeros(m)
        for i in range(1, N):
            j = i - 1
            explicit_part = ((1 - theta) * sigma * (u[i + 1] - 2 * u[i] + u[i - 1]) + u[i])
            d[j] = explicit_part + tau * f[i]

        d[0] += theta * sigma * left_boundary
        d[m - 1] += theta * sigma * right_boundary

        u_inner = tridiagonal_solve(a_sub, b_diag, c_sup, d)

        u_new = np.zeros(N + 1)
        u_new[0] = left_boundary
        u_new[N] = right_boundary
        u_new[1:N] = u_inner

        u = u_new

    return x, results


# Постобработка

def plot_results(x, results, scheme_name):
    plt.figure(figsize=(10, 6))

    for t, u_num in sorted(results.items()):
        plt.plot(x, u_num, 'o-', markersize=3, label=f'{scheme_name} t={t:.3f}')
        u_exact = analytic_solution(x, t)
        plt.plot(x, u_exact, '--', alpha=0.7, label=f'Analytic t={t:.3f}')

    plt.title(f'{scheme_name}')
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()

def compute_error_vs_time(scheme_func, choice):
    N = 50
    L = 1.0
    h = L / N
    tau = 0.0001
    T_max = -np.log(0.05) / (np.pi**2)
    K = int(np.ceil(T_max / tau))

    times = np.arange(0, K + 1) * tau
    errors = []

    x = np.linspace(0, L, N + 1)
    if choice == "1":
        x, results = explicit_scheme()
    elif choice == "2":
        x, results = implicit_scheme()
    elif choice == "3":
        x, results = crank_nicolson_scheme()

    for t in sorted(results.keys()):
        u_num = results[t]
        u_exact = analytic_solution(x, t)
        error = np.max(np.abs(u_num - u_exact))
        errors.append(error)

    # график
    plt.figure(figsize=(8, 5))
    plt.plot(sorted(results.keys()), errors, 'o-', label=scheme_name)
    plt.xlabel('Время T')
    plt.ylabel('Максимальная ошибка')
    plt.title(f'Зависимость ошибки от времени ({scheme_name})')
    plt.grid(True, alpha=0.6)
    plt.legend()
    plt.show()


# Main

print(f"Параметры: N={N}, h={h:.4f}, tau={tau:.6f}, T={T:.4f}")
print(f"Число Куранта: sigma={sigma:.4f}")

print("1 - Явная, 2 - Неявная, 3 - Кранк-Николсон")
choice = input()

if choice == "1":
    scheme_name = "Явная схема"
    x, results = explicit_scheme()
    compute_error_vs_time(explicit_scheme, choice)
elif choice == "2":
    scheme_name = "Неявная схема"
    x, results = implicit_scheme()
    compute_error_vs_time(implicit_scheme, choice)
elif choice == "3":
    scheme_name = "Схема Кранка-Николсона"
    x, results = crank_nicolson_scheme()
    compute_error_vs_time(crank_nicolson_scheme, choice)

print("\nРезультаты:")

for t, u_num in sorted(results.items()):
    u_exact = analytic_solution(x, t)
    error = np.max(np.abs(u_num - u_exact))
    print(f"t={t:.4f}: max_error={error:.2e}")


indices = [2, 5, 6]
results_subset = {t: results[t] for i, t in enumerate(sorted(results.keys())) if i in indices}

plot_results(x, results_subset, scheme_name)
