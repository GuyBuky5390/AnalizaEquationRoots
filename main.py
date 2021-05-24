from sympy import symbols, diff
from sympy.utilities.lambdify import lambdify

EXP = 0.0001  # Epsilon


# returns true if there is an intersection point for function f
# between start and end
def has_intersection_point(f, start, end):
    return True if f(start) * f(end) < 0 else False


# returns a list of tuples when each tuple is (x, f(x), f'(x))
def find_solutions(start_point, end_point, f, diff_f):
    lst = []
    while start_point <= end_point + EXP:
        y = f(start_point)
        y_tag = diff_f(start_point)
        lst.append((start_point, y, y_tag))
        start_point += 0.1
    return lst


def bisection_method(a, b, f, diff_f):
    i = 0
    if diff_f(a) * diff_f(b) < 0:
        f = diff_f

    while abs(b) - abs(a) > EXP:
        c = (a + b) / 2
        f_a, f_b, f_c = f(a), f(b), f(c)
        if f_a * f_c > 0:
            a = c
        else:
            b = c

        c = (a + b) / 2
        i += 1
    return c, i


def newton_raphson_method(a, b, f, diff_f):
    i = 0
    x = (a + b) / 2
    f_x = f(x)
    f_tag_x = diff_f(x)

    while f_x > EXP:
        x -= f_x / f_tag_x
        f_x = f(x)
        f_tag_x = diff_f(x)
        i += 1
    return x, i


def secant_method(a, b, f, diff_f):
    i = 0
    x = a
    next_x = b
    while abs(next_x - x) > EXP:
        temp = next_x
        next_x = (x * f(next_x) - next_x * f(x)) / (f(next_x) - f(x))
        x = temp
        i += 1

    return x, i


# the main function
def main():
    # here we create the function and labmdifying it so we can invoke
    x = symbols('x')
    f = (x ** 3 - x ** 2 - 2 * x - 4)
    print(f"The function: {f}")
    polinom = lambdify(x, f)
    derivative = lambdify(x, diff(f, x))

    # start and end point
    start_point = int(input("Enter start point: "))
    end_point = int(input("Enter end point: "))

    # counter for amount of roots
    count = 0

    # if there is no intersection point between start and end point
    # we exit the program
    if not has_intersection_point(polinom, start_point, end_point):
        print("Intersection point does not exist! exiting ..")
        exit(1)

    # if there is, we choose one of the methods
    print("\nIntersection point exists!")
    option = int(input("Choose method:\n1) Bisection method\n2) Newton Raphson\n3) Secant method\nChoose option: "))
    if option == 1:
        print("BiSection Method\n")
        selected_method = bisection_method
    elif option == 2:
        print("Newton Raphson Method\n")
        selected_method = newton_raphson_method
    else:
        print("Secant Method\n")
        selected_method = secant_method

    # gets the list of solutions for given x from start to end
    # intervals of 0.1
    solutions_list = find_solutions(start_point, end_point, polinom, derivative)  # list of tuples -> (x, f(x), f'(x))
    for i in range(1, len(solutions_list)):

        # if there is a sign change (+-) from f(x) and f(x-1)
        # or f'(x) and f'(x-1)
        # we find the root and add 1 to the counter
        if solutions_list[i][1] * solutions_list[i - 1][1] < 0 or solutions_list[i][2] * solutions_list[i - 1][2] < 0:
            root, iterations = selected_method(solutions_list[i - 1][0], solutions_list[i][0], polinom, derivative)
            if polinom(root) <= EXP:
                count += 1
                print(f"Root #{count} : {root}\nNum of iterations : {iterations}\n")


if __name__ == '__main__':
    main()
