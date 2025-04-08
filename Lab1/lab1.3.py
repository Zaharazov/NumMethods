import math

def norm(vector):
    """Евклидова норма вектора."""
    return math.sqrt(sum(x ** 2 for x in vector))


def check_solution(A, b, x):
    n = len(A)

    # Вычисляем b_calc = A * x
    b_calc = [sum(A[i][j] * x[j] for j in range(n)) for i in range(n)]

    # Выводим результаты
    print("\nПроверка решения:")
    for i in range(n):
        print(f"Ожидаемое: {b[i]:.6f}, Полученное: {b_calc[i]:.6f}")


def check1(A):
    """Проверка на диагональное преобладание по строкам или столбцам."""
    n = len(A)

    # Проверка по строкам
    def check_rows():
        for i in range(n):
            if abs(A[i][i]) <= sum(abs(A[i][j]) for j in range(n) if i != j):
                return False
        return True

    # Проверка по столбцам
    def check_columns():
        for j in range(n):
            if abs(A[j][j]) <= sum(abs(A[i][j]) for i in range(n) if i != j):
                return False
        return True

    return check_rows() or check_columns()

def method_of_simple_iterations(A, b, epsilon=1e-6, max_iter=2000):
    """Метод простых итераций."""
    assert check1(A)
    n = len(A)

    # Начальное приближение
    x = b[:]
    prev = x[:]

    na = max(sum(abs(A[i][j] / A[i][i]) for j in range(n) if i != j) for i in range(n))

    for _ in range(max_iter):
        for i in range(n):
            x[i] = b[i] / A[i][i] - sum((A[i][j] / A[i][i]) * prev[j] for j in range(n) if j != i)

        if (na / (1 - na)) * norm([x[i] - prev[i] for i in range(n)]) < epsilon:
            return x

        prev = x[:]

    print("Метод не сошелся за максимальное количество итераций.")
    return x


def method_of_zeidel(A, b, epsilon=1e-6, max_iter=2000):
    """Метод Зейделя."""
    assert check1(A)
    n = len(A)

    # Начальное приближение
    x = b[:]
    prev = x[:]

    na = max(sum(abs(A[i][j] / A[i][i]) for j in range(n) if i != j) for i in range(n))
    nc = max(sum(abs(A[i][j]) for j in range(n) if i != j) / abs(A[i][i]) for i in range(n))

    for _ in range(max_iter):
        for i in range(n):
            sum1 = sum((A[i][j] / A[i][i]) * x[j] for j in range(i))
            sum2 = sum((A[i][j] / A[i][i]) * prev[j] for j in range(i + 1, n))
            x[i] = b[i] / A[i][i] - sum1 - sum2

        if (na < 1 and nc / (1 - na)) * norm([x[i] - prev[i] for i in range(n)]) < epsilon or (
                na >= 1 and norm([x[i] - prev[i] for i in range(n)]) < epsilon):
            return x

        prev = x[:]

    print("Метод не сошелся за максимальное количество итераций.")
    return x


A = [
    [12, -3, -1, 3],
    [5, 20, 9, 1],
    [6, -3, -21, -7],
    [8, -7, 3, -27]
]
b = [-31, 90, 119, 71]

epsilon = 1e-6

x_simple = method_of_simple_iterations(A, b, epsilon)
print("Решение методом простых итераций:", [round(x, 4) for x in x_simple])
check_solution(A, b, x_simple)
print("------------------------------------------------")
x_zeidel = method_of_zeidel(A, b, epsilon)
print("Решение методом Зейделя:", [round(x, 4) for x in x_zeidel])
check_solution(A, b, x_zeidel)