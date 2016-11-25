import numpy as np
from functools import partial


def wrapper_func(func, args, x):
    """
    Define a wrapper function so that we can freeze some of the arguments
    using partial.
    :param func: The function that is being wrapped.
    :param args: arguments of function @func
    :param x: The input to the function
    :return:
    """
    return func(x, *args)


def check_gradient(obj_func, grad_func, x, args=[],
                   epsilon=0.0001, tolerance=0.001):
    """
    Checks if the gradient function matches the approximate gradient computed
    using the objective function @obj_func
    :param obj_func:
    :param grad_func:
    :param x:
    :param args:
    :param epsilon:
    :return: True if the gradient matches, False otherwise
    """

    d = x.size
    E = np.eye(d) * epsilon
    X = np.tile(x, (d, 1))
    Xplus = X + E
    Xminus = X - E

    f_objective = partial(wrapper_func, obj_func, args)
    f_gradient = partial(wrapper_func, grad_func, args)

    f2 = np.apply_along_axis(f_objective, 1, Xplus)
    f1 = np.apply_along_axis(f_objective, 1, Xminus)
    approx_grad = (f2 - f1)/(2 * epsilon)

    grad = f_gradient(x)

    if grad.shape != approx_grad.shape:
        print "Gradient:", grad.round(4)
        print "Approx Gradient:", approx_grad.round(4)
        raise Exception('Shapes of Gradient and Approximate Gradient dont match')

    max_diff = np.max(np.abs(grad - approx_grad))
    if max_diff > tolerance:
        print "Gradient:", grad.round(4)
        print "Approx Gradient:", approx_grad.round(4)
        print "Maximum absolute difference between grad and approx grad:", round(max_diff, 6)
        return False

    return True

