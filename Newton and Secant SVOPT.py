import sympy as sp
x = sp.symbols('x')
# -------------------------
# Newton Method
# -------------------------
def newton_method(expr, x0, tol=1e-5, max_iter=10):
    f = sp.sympify(expr)
    f_prime = sp.diff(f, x)
    f_double_prime = sp.diff(f_prime, x)

    print("\nNewton's Method (Table Form)")
    headers = ["Iter", "xk", "f'(xk)", "f''(xk)", "x(k+1)", "Check |f'(xk)|<ε"]
    print("{:<5}{:<15}{:<15}{:<15}{:<15}{:<20}".format(*headers))
    print("-"*90)

    xk = x0
    for k in range(max_iter):
        f1 = float(f_prime.subs(x, xk))
        f2 = float(f_double_prime.subs(x, xk))
        if f2 == 0:
            print("Division by zero. Stopping.")
            return
        x_new = xk - f1 / f2
        check = "Yes" if abs(f1) < tol else "No"

        print(f"{k:<5}{xk:<15.6f}{f1:<15.6f}{f2:<15.6f}{x_new:<15.6f}{check:<20}")

        if abs(f1) < tol:
            print(f"\n Converged at iteration {k}, minimizer ≈ {x_new:.6f}")
            return
        xk = x_new

    print(f"\n Stopped after {max_iter} iterations, last x ≈ {xk:.6f}")


# -------------------------
# Secant Method
# -------------------------
def secant_method(expr, x0, x1, tol=1e-5, max_iter=10):
    f = sp.sympify(expr)
    f_prime = sp.diff(f, x)

    print("\nSecant Method (Table Form)")
    headers = ["Iter", "x1", "x2", "f'(x1)", "f'(x2)", "z", "f'(z)", "Check |f'(z)|<ε", "New Region"]
    print("{:<5}{:<12}{:<12}{:<15}{:<15}{:<12}{:<15}{:<18}{:<15}".format(*headers))
    print("-"*120)

    for k in range(max_iter):
        f1 = float(f_prime.subs(x, x0))
        f2 = float(f_prime.subs(x, x1))
        if f2 - f1 == 0:
            print("Division by zero. Stopping.")
            return
        z = x1 - f2 * (x1 - x0) / (f2 - f1)
        fz = float(f_prime.subs(x, z))

        check = "Yes" if abs(fz) < tol else "No"

        if f1 * fz < 0:
            new_region = f"({x0:.6f},{z:.6f})"
            x1 = z
        else:
            new_region = f"({z:.6f},{x1:.6f})"
            x0 = z

        print(f"{k:<5}{x0:<12.6f}{x1:<12.6f}{f1:<15.6f}{f2:<15.6f}{z:<12.6f}{fz:<15.6f}{check:<18}{new_region:<15}")

        if abs(fz) < tol:
            print(f"\n Converged at iteration {k}, minimizer ≈ {z:.6f}")
            return

    print(f"\nStopped after {max_iter} iterations, last z ≈ {z:.6f}")


# -------------------------
# Main Menu Loop
# -------------------------
if __name__ == "__main__":
    print("📌 Single Variable Optimization")
    print("Symbols you can use in function input:")
    print("  Power: x**2   (NOT x^2)")
    print("  Exponential: exp(x)")
    print("  Logarithm: log(x)")
    print("  Trigonometry: sin(x), cos(x), tan(x)")
    print("Examples:  x**2 + 54/x   |   exp(x) - 2*x   |   x**2/2 - sin(x)")
    print()

    while True:
        print("\n==== MENU ====")
        print("1. Newton's Method")
        print("2. Secant Method")
        print("3. Exit")
        choice = input("Enter your choice (1/2/3): ")

        if choice == "1":
            expr = input("\nEnter function f(x): ")
            x0 = float(input("Enter initial guess x0: "))
            tol = float(input("Enter tolerance ε: "))
            newton_method(expr, x0, tol)

        elif choice == "2":
            expr = input("\nEnter function f(x): ")
            x0 = float(input("Enter initial guess x0: "))
            x1 = float(input("Enter initial guess x1: "))
            tol = float(input("Enter tolerance ε: "))
            secant_method(expr, x0, x1, tol)

        elif choice == "3":
            print("\nExiting... ")
            break

        else:
            print(" Invalid choice. Please try again.")
