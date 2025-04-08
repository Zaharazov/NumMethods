import math

# Поработать во время отчета

# Функция для выполнения QR-алгоритма
def qr_algorithm(A, epsilon=1e-6, max_iter=1000):
    A = [row[:] for row in A]
    n = len(A)

    iterations = 0
    # QR-итерации
    for _ in range(max_iter):
        # Разлагаем A
        # Q - ортогональная матрица
        # R - верхняя треугольная матрица
        Q, R = qr_decomposition(A)

        # Обновляем матрицу A
        A = multiply_matrices(R, Q)

        iterations = iterations + 1

        # Проверка на сходимость (элементы ниже диагонали малы)
        off_diagonal = sum(sum(abs(A[i][j]) for j in range(i)) for i in range(1, n))
        if off_diagonal < epsilon:
            break

    # Собственные значения на диагонали
    print(iterations)
    return [A[i][i] for i in range(n)]

# QR-разложение
def qr_decomposition(A):
    n = len(A)
    m = len(A[0])

    Q = [[0] * m for _ in range(n)]  # Ортогональная матрица
    R = [[0] * m for _ in range(m)]  # Верхняя треугольная матрица

    for j in range(m):
        # Берем j-й столбец матрицы A
        v = [A[i][j] for i in range(n)]

        for i in range(j):
            # Вычисляем скалярное произведение векторов v и q_i
            R[i][j] = sum(x1 * x2 for x1, x2 in zip(Q[i], v))

            # Вычитаем проекцию на q_i
            v = [x1 - R[i][j] * x2 for x1, x2 in zip(v, Q[i])]

        # Нормализуем вектор v и помещаем в столбец Q
        R[j][j] = math.sqrt(sum(x ** 2 for x in v))
        Q[j] = [x / R[j][j] for x in v]

    return Q, R

# Функция для умножения двух матриц
def multiply_matrices(A, B):
    m, n = len(A), len(A[0])
    n2, p = len(B), len(B[0])

    if n != n2:
        raise ValueError("Матрицы не могут быть перемножены")

    result = [[0] * p for _ in range(m)]
    for i in range(m):
        for j in range(p):
            result[i][j] = sum(A[i][k] * B[k][j] for k in range(n))

    return result

# Пример использования
A = [
    [-5, -8, 4],
    [4, 2, 6],
    [-2, 5, -6]
]

eigenvalues = qr_algorithm(A)
print("Собственные значения:", eigenvalues)
