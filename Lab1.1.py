def lu_decomposition(A): # преобразуем матрицу А в матрицу L + U
    n = len(A)
    for k in range(n):
        for i in range(k + 1, n):
            A[i][k] /= A[k][k]
            for j in range(k + 1, n):
                A[i][j] -= A[i][k] * A[k][j]


def forward_substitution(A, b): # прямая подстановка
    n = len(A)
    y = [0] * n
    for i in range(n): # считаем первый y, потом при подсчете второго вычитаем предыдущие и т.д.
        y[i] = b[i] - sum(A[i][j] * y[j] for j in range(i))
    return y


def backward_substitution(A, y): # обратная подстановка
    n = len(A)
    x = [0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(A[i][j] * x[j] for j in range(i + 1, n))) / A[i][i] # аналогично + делим на коэф перед x
    return x


def solve_lu(A, b): # решаем задачу
    lu_decomposition(A)
    y = forward_substitution(A, b) # Ly=b
    x = backward_substitution(A, y) # Ux=y
    return x


def determinant(A): # определитель = произведение всех диагональных элементов матрицы U
    lu_decomposition(A)
    det = 1
    for i in range(len(A)):
        det *= A[i][i]
    return det


def inverse_matrix(A):
    n = len(A)
    lu_decomposition(A)
    inv_A = [[0] * n for _ in range(n)]
    for i in range(n):
        e = [1 if j == i else 0 for j in range(n)]
        y = forward_substitution(A, e) # L⋅y=e
        x = backward_substitution(A, y) # U⋅x=y
        for j in range(n):
            inv_A[j][i] = x[j] # фомируем матрицу по столбцам
    return inv_A


def verify_lu(A_original):
    n = len(A_original)
    A = [row[:] for row in A_original]
    lu_decomposition(A)

    result = [[sum((1 if i == k else A[i][k] if i > k else 0) * (A[k][j] if k <= j else 0) for k in range(n)) for j in
               range(n)] for i in range(n)] # произведение L на U, составленных из матрицы A
    return result

A = [[-7, -9, 1, -9],
     [-6, -8, -5, 2],
     [-3, 6, 5, -9],
     [-2, 0, -5, -9]]
b = [29, 42, 11, 75]

x = solve_lu([row[:] for row in A], b)
print("Решение СЛАУ:")
print([round(value, 6) for value in x])

det_A = determinant([row[:] for row in A])
print("Определитель:", round(det_A, 6))

inv_A = inverse_matrix([row[:] for row in A])
print("Обратная матрица:")
for row in inv_A:
    print(row)

print("Проверка LU-разложения (L * U):")
LU_product = verify_lu(A)
for row in LU_product:
    print([round(value, 6) for value in row])
