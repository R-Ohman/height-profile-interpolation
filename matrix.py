from typing import Union


class Matrix:
    def __init__(self, matrix: list[list[int]]):
        """
        Initialize the matrix.
        :param matrix: 2D list of integers
        """
        self._matrix = matrix
        self._rows = len(matrix)
        self._columns = len(matrix[0])

    def __getitem__(self, item: int) -> list[int]:
        """
        Get a row of the matrix.
        :param item: index of the row
        :return: list[int]: row of the matrix
        """
        return self._matrix[item]

    def _arithmetical_operation(self, other: 'Matrix', operation: callable) -> 'Matrix':
        """
            Apply an operation (addition or subtraction) between two matrices.
            :param other: Matrix object
            :param operation: A callable representing the operation to apply
            :return: Matrix: result of the operation
            """
        if self._rows != other.rows or self._columns != other.columns:
            raise ValueError("Matrices must have the same dimensions")

        return Matrix([[operation(self._matrix[i][j], other.matrix[i][j])
                        for j in range(self._columns)] for i in range(self._rows)])

    def __add__(self, other: 'Matrix') -> 'Matrix':
        """
        Add another matrix to this matrix.
        :param other: Matrix object
        :return: Matrix: sum of two matrices
        """
        return self._arithmetical_operation(other, lambda x, y: x + y)

    def __sub__(self, other: 'Matrix') -> 'Matrix':
        """
        Subtract another matrix from this matrix.
        :param other: Matrix object
        :return: Matrix: difference of two matrices
        """
        return self._arithmetical_operation(other, lambda x, y: x - y)

    def __mul__(self, other: Union['Matrix', int]) -> 'Matrix':
        """
        Multiply this matrix by another matrix or a scalar.
        :param other: Matrix object or integer scalar
        :return: Matrix: product of two matrices or scalar multiplication
        """
        if isinstance(other, int):
            return Matrix([[element * other for element in row] for row in self._matrix])

        if isinstance(other, Matrix):
            if self._columns != other.rows:
                raise ValueError("Number of columns of the first matrix must be equal to the number of rows of the second matrix")

            result = Matrix.create_matrix(self._rows, other.columns)
            for i in range(self._rows):
                for j in range(other.columns):
                    for k in range(self._columns):
                        result.matrix[i][j] += self._matrix[i][k] * other.matrix[k][j]
            return result

        raise ValueError("Invalid type for multiplication")

    def __copy__(self):
        """
        Create a copy of the matrix.
        :return: Matrix: copy of the matrix
        """
        return Matrix([row.copy() for row in self._matrix])

    def lu_decomposition(self) -> tuple['Matrix', 'Matrix', 'Matrix']:
        """
        Perform LU decomposition of the matrix.
        :return: tuple[Matrix, Matrix, Matrix]: L, U, P matrices
        """
        n = self._rows
        l = Matrix.create_diagonal(n)
        u = self.__copy__()
        p = Matrix.create_diagonal(n)

        for i in range(n):
            Matrix.pivot(l, u, p, i)
            for j in range(i + 1, n):
                l[j][i] = u[j][i] / u[i][i]
                for k in range(i, n):
                    u[j][k] -= l[j][i] * u[i][k]
        return l, u, p

    def solve(self, b: 'Matrix') -> 'Matrix':
        """
        Solve the system of linear equations using LU decomposition.
        :param b: Matrix object
        :return: Matrix: solution of the system of linear equations
        """
        l, u, p = self.lu_decomposition()
        n = self._rows
        y = Matrix.create_matrix(n, 1)
        x = Matrix.create_matrix(n, 1)

        # Solve Ly = Pb
        for i in range(n):
            y[i][0] = b[p[i][0]][0] - sum(l[i][j] * y[j][0] for j in range(i))

        # Solve Ux = y
        for i in range(n - 1, -1, -1):
            x[i][0] = y[i][0] - sum(u[i][j] * x[j][0] for j in range(i + 1, n))
            x[i][0] /= u[i][i]

        return x

    @property
    def rows(self) -> int:
        return self._rows

    @property
    def columns(self) -> int:
        return self._columns

    @property
    def matrix(self) -> list[list[int]]:
        return self._matrix

    @staticmethod
    def create_matrix(rows: int, columns: int, value: int = 0) -> 'Matrix':
        """
        Create a matrix with the specified dimensions and fill it with the specified value.
        :param rows: number of rows
        :param columns: number of columns
        :param value: value to fill the matrix with
        :return: Matrix: matrix with the specified dimensions and filled with the specified value
        """
        return Matrix([[value for _ in range(columns)] for _ in range(rows)])

    @staticmethod
    def create_diagonal(n: int, value: float = 1) -> 'Matrix':
        """
        Create a diagonal matrix with the specified value on the diagonal.
        :param n: size of the matrix
        :param value: value to fill the diagonal with
        :return: Matrix: diagonal matrix with the specified value on the diagonal
        """
        matrix = Matrix.create_matrix(n, n)
        for i in range(n):
            matrix[i][i] = value
        return matrix

    @staticmethod
    def from_list(vector: list) -> 'Matrix':
        """
        Create a matrix from a list of lists.
        :param vector: 2D list of integers
        :return: Matrix: matrix created from the list
        """
        return Matrix([[element] for element in vector])

    def flatten(self) -> list:
        """
        Flatten the matrix into a list.
        :return: list: flattened matrix
        """
        return [element for row in self._matrix for element in row]

    @staticmethod
    def pivot(l: 'Matrix', u: 'Matrix', p: 'Matrix', i: int) -> None:
        n = u.columns
        max_val = 0
        max_ind = i
        for k in range(i, n):
            if abs(u[k][i]) > max_val:
                max_val = abs(u[k][i])
                max_ind = k

        if max_ind != i:
            u[i], u[max_ind] = u[max_ind], u[i]
            l[i], l[max_ind] = l[max_ind], l[i]
            p[i], p[max_ind] = p[max_ind], p[i]