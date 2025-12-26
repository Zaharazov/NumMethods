import numpy as np
import matplotlib.pyplot as plt

# Предобработка

def left_boundary(y):
    return np.cos(y)

def right_boundary(y):
    return np.zeros_like(y)

def bottom_boundary(x):
    return np.cos(x)

def top_boundary(x):
    return np.zeros_like(x)

def analytic_solution(x, y):
    return np.cos(x) * np.cos(y)

Lx = np.pi / 2
Ly = np.pi / 2
Nx = 12
Ny = 12
hx = 0.13
hy = 0.13

max_iter = 5000
eps = 1e-6
omega = 1.6

# Солвер

def libman_method(Nx, Ny, hx, hy, max_iter, eps):
    x = np.linspace(0, Lx, Nx+1)
    y = np.linspace(0, Ly, Ny+1)
    u = np.zeros((Ny+1, Nx+1))
    u[0, :]  = bottom_boundary(x)
    u[-1, :] = top_boundary(x)
    u[:, 0]  = left_boundary(y)
    u[:, -1] = right_boundary(y)

    Ax = 1.0 / hx**2
    Ay = 1.0 / hy**2

    snapshots = {}
    errors = []

    X, Y = np.meshgrid(x, y)
    U_exact = analytic_solution(X, Y)

    for it in range(max_iter):
        u_new = u.copy()
        max_err = 0.0

        for j in range(1, Ny):
            for i in range(1, Nx):
                u_new[j, i] = (Ax*(u[j, i+1] + u[j, i-1]) + Ay*(u[j+1, i] + u[j-1, i]) + 2*u[j, i]) / (2 * (Ax + Ay))
                max_err = max(max_err, abs(u_new[j, i] - u[j, i]))

        u = u_new

        errors.append(np.max(np.abs(u - U_exact)))

        if it == 1 or it == 100:
            snapshots[it] = u.copy()

        if max_err < eps:
            snapshots[it] = u.copy()
            break

    return x, y, u, snapshots, errors


def zeidel_method(Nx, Ny, hx, hy, max_iter, eps):
    x = np.linspace(0, Lx, Nx+1)
    y = np.linspace(0, Ly, Ny+1)
    u = np.zeros((Ny+1, Nx+1))
    u[0, :] = bottom_boundary(x)
    u[-1, :] = top_boundary(x)
    u[:, 0] = left_boundary(y)
    u[:, -1] = right_boundary(y)

    Ax = 1.0 / hx**2
    Ay = 1.0 / hy**2

    snapshots = {}
    errors = []

    X, Y = np.meshgrid(x, y)
    U_exact = analytic_solution(X, Y)

    for it in range(max_iter):
        max_err = 0.0
        for j in range(1, Ny):
            for i in range(1, Nx):
                old = u[j, i]
                u[j, i] = (Ax*(u[j, i+1] + u[j, i-1]) + Ay*(u[j+1, i] + u[j-1, i]) + 2*old) / (2 * (Ax + Ay))
                max_err = max(max_err, abs(u[j, i] - old))

        errors.append(np.max(np.abs(u - U_exact)))

        if it == 1 or it == 100:
            snapshots[it] = u.copy()
        if max_err < eps:
            snapshots[it] = u.copy()
            break

    return x, y, u, snapshots, errors

def sor_method(Nx, Ny, hx, hy, omega, max_iter, eps):
    x = np.linspace(0, Lx, Nx+1)
    y = np.linspace(0, Ly, Ny+1)
    u = np.zeros((Ny+1, Nx+1))
    u[0, :] = bottom_boundary(x)
    u[-1, :] = top_boundary(x)
    u[:, 0] = left_boundary(y)
    u[:, -1] = right_boundary(y)

    Ax = 1.0 / hx**2
    Ay = 1.0 / hy**2

    snapshots = {}
    errors = []

    X, Y = np.meshgrid(x, y)
    U_exact = analytic_solution(X, Y)

    for it in range(max_iter):
        max_err = 0.0
        for j in range(1, Ny):
            for i in range(1, Nx):
                old = u[j, i]
                u_temp = (Ax*(u[j, i+1] + u[j, i-1]) + Ay*(u[j+1, i] + u[j-1, i]) + 2*old ) / (2 * (Ax + Ay))
                u[j, i] = (1 - omega) * old + omega * u_temp
                max_err = max(max_err, abs(u[j, i] - old))

        errors.append(np.max(np.abs(u - U_exact)))

        if it == 1 or it == 25:
            snapshots[it] = u.copy()
        if max_err < eps:
            snapshots[it] = u.copy()
            break

    return x, y, u, snapshots, errors

# Постобработка

def plot_3d_snapshots(x, y, snapshots, method_name):
    X, Y = np.meshgrid(x, y)
    keys = sorted(snapshots.keys())
    n_snap = len(keys)

    U_analytic = analytic_solution(X, Y)

    fig = plt.figure(figsize=(6 * n_snap, 5))

    for idx, k in enumerate(keys):
        ax = fig.add_subplot(1, n_snap, idx + 1, projection='3d')

        ax.plot_surface(X, Y, snapshots[k], cmap='viridis', alpha=0.8)

        ax.plot_surface(X, Y, U_analytic, cmap='coolwarm', alpha=0.3)
        ax.plot_wireframe(X, Y, U_analytic, color='green', linewidth=0.5)

        ax.set_title(f'{method_name}, шаг {k}')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('u')

    plt.show()


# Main

print("1 - простые итерации (либмана), 2 - зейдель, 3 - простые итерации с верхней релаксацией (SOR)")
choice = input()

if choice == "1":
    method_name = "Метод Либмана"
    x, y, u, snaps, errors = libman_method(Nx, Ny, hx, hy, max_iter, eps)

elif choice == "2":
    method_name = "Метод Зейделя"
    x, y, u, snaps, errors = zeidel_method(Nx, Ny, hx, hy, max_iter, eps)

elif choice == "3":
    method_name = "Метод SOR"
    x, y, u, snaps, errors = sor_method(Nx, Ny, hx, hy, omega, max_iter, eps)

plot_3d_snapshots(x, y, snaps, method_name)

plt.figure(figsize=(6,4))
plt.semilogy(errors)
plt.xlabel("Итерация")
plt.ylabel("Максимальная ошибка")
plt.grid(True)
plt.title(method_name)
plt.show()
