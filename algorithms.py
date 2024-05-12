from math import cos, pi
from typing import Callable

from matrix import Matrix
import pandas as pd


def get_selection_points(nodes: list[tuple[float, float]], strategy: Callable[[float, float, int], list[float]], n: int) -> list[tuple[float, float]]:
    """
    Returns the selection points based on the strategy.
    """
    correct_nodes = strategy(nodes[0][0], nodes[-1][0], n)

    # Find the closest nodes to correct_nodes from nodes.
    result = []
    for node in correct_nodes:
        closest_node = min(nodes, key=lambda x: abs(x[0] - node))
        if closest_node not in result:
            result.append(closest_node)
    return result


def get_interpolation_points(nodes: list[tuple[float, float]], n: int) -> list[float]:
    """
    Returns the interpolation points.
    """
    return linspace(nodes[0][0], nodes[-1][0], n)


def read_nodes(path: str) -> list[tuple[float, float]]:
    """
    Read nodes from a CSV file.
    """
    return [(row.iloc[0], row.iloc[1]) for _, row in pd.read_csv(path).iterrows()]


def linspace(start: float, stop: float, num: int) -> list:
    """
    Create an array of equally spaced values.
    :param start: start value
    :param stop: end value
    :param num: number of values
    :return: list: array of equally spaced values
    """
    return [start + (stop - start) * i / (num - 1) for i in range(num)]


def chebyshev_nodes(start: float, stop: float, num: int) -> list[float]:
    """
    Create an array of Chebyshev nodes.
    :param start: start value
    :param stop: end value
    :param num: number of nodes
    :return: list: array of Chebyshev nodes
    """
    return [(start + stop) / 2 + (stop - start) / 2 * cos((2 * i + 1) / (2 * num) * pi) for i in range(num)]


def lagrange_interpolation(data: list[tuple[float, float]], xs: list[float]) -> list[float]:
    """
    Interpolate the data using Lagrange interpolation.
    :param data: list of tuples of x and y values
    :param xs: list of x values to interpolate
    :return: list: list of interpolated y values
    """
    n = len(data)
    ys = [0] * len(xs)

    for m in range(len(xs)):
        for i in range(n):
            xi, yi = data[i]
            term = yi
            for j in range(n):
                if j != i:
                    xj, _ = data[j]
                    term *= (xs[m] - xj) / (xi - xj)
            ys[m] += term
    return ys

