import numpy as np
import matplotlib.pyplot as plt

# Предобработка

def left_boundary(t):
    return -np.sin(a * t)

def right_boundary(t):
    return np.sin(a * t)

def u0(x):
    return np.sin(x)

def ut0(x):
    return -a * np.cos(x)

def analytic_solution(x, t):
    return np.sin(x - a * t)

a = 6.0
L = np.pi
T = 10.0
N = 50
h = L / N
tau = 0.0005
K = int(np.ceil(T/tau))
sigma = (a * tau / h)**2

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

    return x


def explicit_scheme(order):
    x = np.linspace(0, L, N + 1)
    u = np.zeros((K + 1, N + 1))

    u[0, :] = u0(x)

    if order == 1:
        u[1, :] = u[0, :] + tau * ut0(x)
    else:
        u_tt = -a**2 * np.sin(x)
        u[1, :] = u[0, :] + tau * ut0(x) + 0.5 * tau**2 * u_tt

    u[0, 0] = left_boundary(0)
    u[0, -1] = right_boundary(0)
    u[1, 0] = left_boundary(tau)
    u[1, -1] = right_boundary(tau)

    for k in range(1, K):
        for i in range(1, N):
            u[k+1, i] = 2 * u[k, i] - u[k-1, i] + sigma * (u[k, i+1] - 2*u[k, i] + u[k, i-1])

        u[k+1, 0] = left_boundary((k+1)*tau)
        u[k+1, -1] = right_boundary((k+1)*tau)

    return x, u


def implicit_scheme(order):
    x = np.linspace(0, L, N + 1)
    u = np.zeros((K + 1, N + 1))
    u[0, :] = u0(x)
    if order == 1:
        u[1, :] = u[0, :] + tau * ut0(x)
    else:
        u_tt = -a**2 * np.sin(x)
        u[1, :] = u[0, :] + tau * ut0(x) + 0.5 * tau**2 * u_tt

    for k in range(1, K):
        n = N - 1
        a_sub = -sigma * np.ones(n-1)
        b_diag = (1 + 2 * sigma) * np.ones(n)
        c_sub = -sigma * np.ones(n-1)

        d = 2*u[k,1:N] - u[k-1,1:N]
        d[0] += sigma * left_boundary((k+1)*tau)
        d[-1] += sigma * right_boundary((k+1)*tau)

        u_new_inner = tridiagonal_solve(a_sub, b_diag, c_sub, d)
        u[k+1, 1:N] = u_new_inner
        u[k+1, 0] = left_boundary((k+1)*tau)
        u[k+1, -1] = right_boundary((k+1)*tau)
    return x, u

# Постобработка

def compute_error(u, x):
    errors = []
    times = np.linspace(0, T, K+1)
    for k, t in enumerate(times):
        exact = analytic_solution(x, t)
        errors.append(np.max(np.abs(u[k,:] - exact)))
    return times, errors

# Main

print(f"Параметры: N={N}, h={h:.4f}, tau={tau:.6f}, T={T:.4f}")
print(f"sigma={sigma:.6f}")

print("1 - явная, 2 - неявная")
choice = input()
print("Порядок аппроксимации второго начального условия: ")
order = int(input())

if choice == "1":
    scheme_name = "Явная схема"
    x, u = explicit_scheme(order)
else:
    scheme_name = "Неявная схема"
    x, u = implicit_scheme(order)

times, errors = compute_error(u, x)


plt.figure(figsize=(12, 5))
plt.plot(times, errors, label=f'{scheme_name}, порядок {order}')
plt.xlabel("t")
plt.ylabel("Максимальная ошибка")
plt.title("Зависимость ошибки от времени")
plt.grid(True)
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))

for t_plot in [T/4, T/2, T]:
    k = int(t_plot / tau)
    plt.plot(x, u[k], 'o-', markersize=3, label=f'u(x,t), t={t_plot:.2f}')
    plt.plot(x, analytic_solution(x, t_plot), '--', alpha=0.8, label=f'U(x,t), t={t_plot:.2f}')

plt.xlabel("x")
plt.ylabel("u(x,t)")
plt.title(f"{scheme_name}, порядок {order}")
plt.grid(True)
plt.legend()
plt.show()
