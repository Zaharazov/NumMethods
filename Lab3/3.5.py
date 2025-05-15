import math

def f(x):
    return x / (x**2 + 9)

def rectangle_method(a, b, h):
    n = int((b - a) / h)
    result = 0
    for i in range(0, n):
        xi = a + (i+0.5) * h
        result += f(xi)
    return result * h

def trapez_method(a, b, h):
    n = int((b - a) / h)
    result = 0.5 * (f(a)+f(b))
    for i in range(1, n):
        xi = a + i * h
        result += f(xi)
    return result * h

def simpson_method(a, b, h):
    n = int((b - a) / h)

    result = f(a) + f(b)
    for i in range(1, n):
        xi = a + i * h
        coeff = 4 if i % 2 != 0 else 2
        result += coeff * f(xi)
    return result * h / 3

def runge_romberg(Fh, Fkh, k, p):
    return Fh + (Fh - Fkh) / (k**p - 1)

def absolute_error(F_exact, F_approx):
    return abs(F_exact - F_approx)

# Входные данные
a = 0
b = 2
h1 = 0.5
h2 = 0.25
k = int(h1 / h2)
F_true = 0.5 * math.log(13 / 9)

# Метод прямоугольников
F_rect_h1 = rectangle_method(a, b, h1)
F_rect_h2 = rectangle_method(a, b, h2)
F_rect_rr = runge_romberg(F_rect_h2, F_rect_h1, k, p=2)

# Метод трапеций
F_trap_h1 = trapez_method(a, b, h1)
F_trap_h2 = trapez_method(a, b, h2)
F_trap_rr = runge_romberg(F_trap_h2, F_trap_h1, k, p=2)

# Метод Симпсона
F_simp_h1 = simpson_method(a, b, h1)
F_simp_h2 = simpson_method(a, b, h2)
F_simp_rr = runge_romberg(F_simp_h2, F_simp_h1, k, p=4)

# Вывод
print("Точное значение интеграла: {:.8f}\n".format(F_true))

print("Метод прямоугольников:")
print(f"  h = {h1}: {F_rect_h1:.8f}")
print(f"  h = {h2}: {F_rect_h2:.8f}")
print(f"  Рунге-Ромберг: {F_rect_rr:.8f}")
print(f"  Погрешность: {absolute_error(F_true, F_rect_rr):.8f}\n")

print("Метод трапеций:")
print(f"  h = {h1}: {F_trap_h1:.8f}")
print(f"  h = {h2}: {F_trap_h2:.8f}")
print(f"  Рунге-Ромберг: {F_trap_rr:.8f}")
print(f"  Погрешность: {absolute_error(F_true, F_trap_rr):.8f}\n")

print("Метод Симпсона:")
print(f"  h = {h1}: {F_simp_h1:.8f}")
print(f"  h = {h2}: {F_simp_h2:.68f}")
print(f"  Рунге-Ромберг: {F_simp_rr:.8f}")
print(f"  Погрешность: {absolute_error(F_true, F_simp_rr):.8f}")
