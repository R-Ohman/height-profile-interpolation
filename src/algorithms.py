from math import cos, pi
from typing import Callable

from matrix import Matrix
import pandas as pd


def get_x_values(data: list[tuple[float, float]]) -> list[float]:
    """
    Returns the x values from data.
    :param data: list of tuples of x and y values
    :return: list: list of x values
    """
    return [x for x, _ in data]


def get_y_values(data: list[tuple[float, float]]) -> list[float]:
    """
    Returns the y values from data.
    :param data: list of tuples of x and y values
    :return: list: list of y values
    """
    return [y for _, y in data]


def get_selection_points(nodes: list[tuple[float, float]], strategy: Callable[[float, float, int], list[float]], n: int) -> list[tuple[float, float]]:
    """
    Returns the selection points based on the strategy.
    :param nodes: list of tuples of x and y values
    :param strategy: strategy to select the points (e.g., Chebyshev nodes)
    :param n: number of selection points
    :return: list: list of selection points
    """
    correct_nodes = strategy(nodes[0][0], nodes[-1][0], n)

    # Find the closest nodes to correct_nodes from nodes.
    result = []
    for node in correct_nodes:
        closest_node = min(nodes, key=lambda x: abs(x[0] - node))
        if closest_node not in result:
            result.append(closest_node)
    return result


def get_linspace_points(nodes: list[tuple[float, float]], n: int) -> list[float]:
    """
    Returns the interpolation points.
    :param nodes: list of tuples of x and y values
    :param n: number of interpolation points
    :return: list: list of interpolation points
    """
    return linspace(nodes[0][0], nodes[-1][0], n)


def read_nodes(path: str) -> list[tuple[float, float]]:
    """
    Read nodes from a CSV file.
    :param path: path to the CSV file
    :return: list: list of tuples of x and y values
    """
    if path.endswith('.csv'):
        return [(row[0], row[1]) for row in pd.read_csv(path).values]

    with open(path, 'r') as file:
        return [(float(x), float(y)) for x, y in (line.strip().split() for line in file)]


def linspace(start: float, stop: float, num: int) -> list:
    """
    Create an array of equally spaced values.
    :param start: start value
    :param stop: end value
    :param num: number of values
    :return: list: array of equally spaced values
    """
    return [start + (stop - start) * i / (num - 1) for i in range(num)]


def get_chebyshev_nodes(start: float, stop: float, num: int) -> list[float]:
    """
    Create an array of Chebyshev nodes.
    :param start: start value
    :param stop: end value
    :param num: number of nodes
    :return: list: array of Chebyshev nodes
    """
    return [(start + stop) / 2 + (stop - start) / 2 * cos((2 * i + 1) / (2 * num) * pi) for i in range(num)]


def lagrange_interpolation(data: list[tuple[float, float]], x_values: list[float]) -> list[float]:
    """
    Interpolate the data using Lagrange interpolation.
    :param data: list of tuples of x and y values
    :param x_values: list of x values to interpolate
    :return: list: list of interpolated y values
    """
    n = len(data)
    y_values = [0] * len(x_values)

    for vi in range(len(x_values)):
        for i in range(n):
            xi, yi = data[i]
            for j in range(n):
                if i != j:
                    xj, _ = data[j]
                    yi *= (x_values[vi] - xj) / (xi - xj)
            y_values[vi] += yi
    return y_values


def cubic_spline_interpolation(data: list[tuple[float, float]], x_values: list[float]) -> list[float]:
    """
    Interpolate the data using cubic spline interpolation.
    :param data: list of tuples of x and y values
    :param x_values: list of x values to interpolate
    :return: list: list of interpolated y values
    """
    subranges_number = len(data) - 1
    # the size of the matrix
    n = 4 * subranges_number

    A = Matrix.create_matrix(n, n)
    b = Matrix.create_matrix(n, 1)

    # the row index
    row_id = 0

    # loop over each subrange
    for i in range(subranges_number):
        # Set the first coefficient of the i-th cubic polynomial to 1
        A[row_id][i * 4] = 1
        # Set the corresponding value in vector b to the y value of the i-th data point
        b[row_id] = data[i][1]
        row_id += 1

        # Calculate the difference in x values between the i-th and (i+1)-th data points
        h = data[i + 1][0] - data[i][0]
        # Set the coefficients of the i-th cubic polynomial for the (i+1)-th data point
        for j in range(4):
            A[row_id][i * 4 + j] = h ** j
        # Set the corresponding value in vector b to the y value of the (i+1)-th data point
        b[row_id] = data[i + 1][1]
        row_id += 1

        # If not the last subrange, set the conditions for the derivatives to be equal at the inner points
        if i < subranges_number - 1:
            # Set the condition for the first derivatives to be equal
            for j in range(1, 4):
                A[row_id][i * 4 + j] = j * h ** (j - 1)
            A[row_id][(i + 1) * 4 + 1] = -1
            row_id += 1

            # Set the condition for the second derivatives to be equal
            A[row_id][i * 4 + 2] = 2
            A[row_id][i * 4 + 3] = 6 * h
            A[row_id][(i + 1) * 4 + 2] = -2
            row_id += 1

    # Set the condition for the second derivative to be zero at the first point
    A[row_id][2] = 1
    row_id += 1

    # Set the condition for the second derivative to be zero at the last point
    A[row_id][n - 2] = 2
    A[row_id][n - 1] = 6 * (data[-1][0] - data[-2][0])
    row_id += 1

    # Solve the system of equations to find the coefficients of the cubic polynomials
    x = A.solve(b).flatten()

    # the list of results
    results = []
    # Loop over each x value to interpolate
    for x_value in x_values:
        found = False
        # Find the subrange that the x value belongs to
        for i in range(subranges_number):
            if data[i][0] <= x_value <= data[i + 1][0]:
                # Calculate the difference in x values
                h = x_value - data[i][0]
                # Calculate the interpolated y value using the coefficients of the i-th cubic polynomial
                y_interpolated = sum(x[i * 4 + j] * h ** j for j in range(4))
                results.append(y_interpolated)
                found = True
                break
        # If the x value is outside the range of the data, use the last calculated y value
        if not found:
            results.append(results[-1])
    return results
