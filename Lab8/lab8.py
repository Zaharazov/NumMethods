import numpy as np
import matplotlib.pyplot as plt

# Предобработка

def left_boundary(y, t, a):
    return np.cosh(y) * np.exp(-3 * a * t)

def right_boundary(y, t, a):
    return np.zeros_like(y)

def bottom_boundary(x, t, a):
    return np.cos(2 * x) * np.exp(-3 * a * t)

def top_boundary(x, t, a):
    return (5/4) * np.cos(2 * x) * np.exp(-3 * a * t)

def ut0(x, y):
    return np.cos(2 * x) * np.cosh(y)

def analytic_solution(x, y, t, a):
    return np.cos(2 * x) * np.cosh(y) * np.exp(-3 * a * t)

a = 1.0
Lx = np.pi / 4
Ly = np.log(2)
T = 1.0
Nx = 15
Ny = 15
Nt = 100
hx = Lx / Nx
hy = Ly / Ny
tau = T / Nt

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


def MPN(Nx, Ny, Nt, hx, hy, tau, a):
    x = np.linspace(0, Lx, Nx + 1)
    y = np.linspace(0, Ly, Ny + 1)
    X, Y = np.meshgrid(x, y)
    u = ut0(X, Y)

    snapshots = {}
    errors = [np.max(np.abs(u - analytic_solution(X, Y, 0, a)))]

    r_x = a * tau / (2 * hx ** 2)
    r_y = a * tau / (2 * hy ** 2)

    snapshot_times = [0.01, T/2, T]

    for n in range(1, Nt + 1):
        t_half = (n - 0.5) * tau

        u_star = u.copy()
        for j in range(1, Ny):
            a_diag = -r_x * np.ones(Nx - 1)
            b_diag = (1 + 2 * r_x) * np.ones(Nx - 1)
            c_diag = -r_x * np.ones(Nx - 1)
            d = r_y * u[j + 1, 1:-1] + (1 - 2 * r_y) * u[j, 1:-1] + r_y * u[j - 1, 1:-1]
            d[0] += r_x * left_boundary(y[j], t_half, a)
            d[-1] += r_x * right_boundary(y[j], t_half, a)
            u_star[j, 1:-1] = tridiagonal_solve(a_diag, b_diag, c_diag, d)

        u_new = u_star.copy()
        for i in range(1, Nx):
            a_diag = -r_y * np.ones(Ny - 1)
            b_diag = (1 + 2 * r_y) * np.ones(Ny - 1)
            c_diag = -r_y * np.ones(Ny - 1)
            d = r_x * u_star[1:-1, i + 1] + (1 - 2 * r_x) * u_star[1:-1, i] + r_x * u_star[1:-1, i - 1]
            d[0] += r_y * bottom_boundary(x[i], t_half, a)
            d[-1] += r_y * top_boundary(x[i], t_half, a)
            u_new[1:-1, i] = tridiagonal_solve(a_diag, b_diag, c_diag, d)

        u = u_new.copy()
        u[:, 0] = left_boundary(y, n * tau, a)
        u[:, -1] = right_boundary(y, n * tau, a)
        u[0, :] = bottom_boundary(x, n * tau, a)
        u[-1, :] = top_boundary(x, n * tau, a)

        U_exact = analytic_solution(X, Y, n * tau, a)
        errors.append(np.max(np.abs(u - U_exact)))

        if any(np.isclose(n*tau, ts, atol=tau/2) for ts in snapshot_times):
            snapshots[n*tau] = u.copy()

    return x, y, u, snapshots, errors


def MDSH(Nx, Ny, Nt, hx, hy, tau, a):
    x = np.linspace(0, Lx, Nx + 1)
    y = np.linspace(0, Ly, Ny + 1)
    X, Y = np.meshgrid(x, y)
    u = ut0(X, Y)

    snapshots = {}
    errors = [np.max(np.abs(u - analytic_solution(X, Y, 0, a)))]

    r_x = a * tau / hx ** 2
    r_y = a * tau / hy ** 2

    snapshot_times = [0.01, T/2, T]

    for n in range(1, Nt + 1):
        t_half = (n - 0.5) * tau

        u_star = u.copy()
        for j in range(1, Ny):
            a_diag = -r_x * np.ones(Nx - 1)
            b_diag = (1 + 2 * r_x) * np.ones(Nx - 1)
            c_diag = -r_x * np.ones(Nx - 1)
            d = u[j, 1:-1]
            d[0] += r_x * left_boundary(y[j], t_half, a)
            d[-1] += r_x * right_boundary(y[j], t_half, a)
            u_star[j, 1:-1] = tridiagonal_solve(a_diag, b_diag, c_diag, d)

        u_new = u_star.copy()
        for i in range(1, Nx):
            a_diag = -r_y * np.ones(Ny - 1)
            b_diag = (1 + 2 * r_y) * np.ones(Ny - 1)
            c_diag = -r_y * np.ones(Ny - 1)
            d = u_star[1:-1, i]
            d[0] += r_y * bottom_boundary(x[i], t_half, a)
            d[-1] += r_y * top_boundary(x[i], t_half, a)
            u_new[1:-1, i] = tridiagonal_solve(a_diag, b_diag, c_diag, d)

        u = u_new.copy()
        u[:, 0] = left_boundary(y, n * tau, a)
        u[:, -1] = right_boundary(y, n * tau, a)
        u[0, :] = bottom_boundary(x, n * tau, a)
        u[-1, :] = top_boundary(x, n * tau, a)

        U_exact = analytic_solution(X, Y, n * tau, a)
        errors.append(np.max(np.abs(u - U_exact)))

        if any(np.isclose(n*tau, ts, atol=tau/2) for ts in snapshot_times):
            snapshots[n*tau] = u.copy()

    return x, y, u, snapshots, errors

# Постобработка

def plot_3d_snapshots(x, y, snapshots, method_name):
    X, Y = np.meshgrid(x, y)
    keys = sorted(snapshots.keys())
    n_snap = len(keys)

    # находим глобальные min и max для Z
    all_values = np.array([snapshots[t] for t in keys])
    z_min = np.min(all_values)
    z_max = np.max(all_values)

    fig = plt.figure(figsize=(6 * n_snap, 5))
    for idx, t in enumerate(keys):
        ax = fig.add_subplot(1, n_snap, idx + 1, projection='3d')
        ax.plot_surface(X, Y, snapshots[t], cmap='viridis', alpha=0.8)
        ax.plot_surface(X, Y, analytic_solution(X, Y, t, a), cmap='coolwarm', alpha=0.3)
        ax.plot_wireframe(X, Y, analytic_solution(X, Y, t, a), color='green', linewidth=0.5)
        ax.set_title(f'{method_name}, t={t:.2f}')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('u')
        ax.set_zlim(z_min, z_max)
    plt.show()

# Main

print("1 - схема переменных направлений, 2 - схема дробных шагов")
choice = input()

if choice == "1":
    method_name = "Схема переменных направлений"
    x, y, u_final, snapshots, errors = MPN(Nx, Ny, Nt, hx, hy, tau, a)
elif choice == "2":
    method_name = "Схема дробных шагов"
    x, y, u_final, snapshots, errors = MDSH(Nx, Ny, Nt, hx, hy, tau, a)

plot_3d_snapshots(x, y, snapshots, method_name)

time_steps = np.linspace(0, T, len(errors))

plt.figure(figsize=(6,4))
plt.plot(time_steps, errors)
plt.xlabel("Время (t)")
plt.ylabel("Максимальная ошибка")
plt.xlim(0, T)
plt.grid(True, which="both", ls="-", alpha=0.5)
plt.title(f"{method_name}")
plt.show()
