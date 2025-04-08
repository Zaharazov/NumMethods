import math

def is_symmetric(matrix):
    n = len(matrix)
    for i in range(n):
        for j in range(n):
            if matrix[i][j] != matrix[j][i]:
                return False
    return True


def max_off_diagonal(A):
    max_val = 0
    p, q = 0, 1
    for i in range(n):
        for j in range(i + 1, n): # только элементы выше диагонали
            if abs(A[i][j]) > max_val:
                max_val = abs(A[i][j])
                p, q = i, j
    return p, q, max_val

def rotate(A, V, p, q):
    if A[p][p] == A[q][q]:
        theta = math.pi / 4
    else:
        theta = 0.5 * (A[q][q] - A[p][p]) / A[p][q]
        t = 1 / (abs(theta) + math.sqrt(theta ** 2 + 1))
        if theta < 0:
            t = -t
        c = 1 / math.sqrt(1 + t ** 2)
        s = t * c

    for i in range(n):
        if i != p and i != q: # не диагональные элементы
            A_ip = A[i][p]
            A_iq = A[i][q]
            A[i][p] = A_ip * c - A_iq * s
            A[p][i] = A[i][p]
            A[i][q] = A_ip * s + A_iq * c
            A[q][i] = A[i][q]

    A_pp = A[p][p]
    A_qq = A[q][q]
    A_pq = A[p][q]

    A[p][p] = c ** 2 * A_pp - 2 * s * c * A_pq + s ** 2 * A_qq # диагональные элементы
    A[q][q] = s ** 2 * A_pp + 2 * s * c * A_pq + c ** 2 * A_qq
    A[p][q] = A[q][p] = 0

    for i in range(n):
        V_ip = V[i][p]
        V_iq = V[i][q]
        V[i][p] = V_ip * c - V_iq * s
        V[i][q] = V_ip * s + V_iq * c


def create_matrix(eigen_values, eigen_vectors):
    """Check: create matrix A from eigenvalues and eigenvectors"""
    n = len(eigen_values)
    A = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                A[i][j] += eigen_vectors[i][k] * eigen_values[k] * eigen_vectors[j][k]
    return A

def print_matrix(matrix):
    """Выводит матрицу в красивом виде."""
    for row in matrix:
        print("  ".join(f"{val:.5f}" for val in row))




A = [
    [4, 7, -1],
    [7, -9, -6],
    [-1, -6, -4]
]

epsilon = 1e-6
n = len(A)
V = [[1 if i == j else 0 for j in range(n)] for i in range(n)] # здесь собственные векторы

if not is_symmetric(A):
    print("Ошибка: матрица не симметрична!")
else:
    iterations = 0
    while True:
        p, q, max_val = max_off_diagonal(A)
        if max_val < epsilon:
            break
        rotate(A, V, p, q)
        iterations += 1

    eigenvalues = [A[i][i] for i in range(n)]
    eigenvectors = V

    print("Собственные значения:", eigenvalues)
    print("Собственные векторы:")
    for vec in eigenvectors:
        print(vec)

    A_reconstructed = create_matrix(eigenvalues, eigenvectors)

    print("Восстановленная матрица A:")
    print_matrix(A_reconstructed)


    # if A[p][p] == A[q][q]:
    #     theta = math.pi / 4
    # else:
    #     theta = 0.5 * (A[q][q] - A[p][p]) / A[p][q]
    #     t = 1 / (abs(theta) + math.sqrt(theta ** 2 + 1))
    #     if theta < 0:
    #         t = -t
    #     c = 1 / math.sqrt(1 + t ** 2)
    #     s = t * c