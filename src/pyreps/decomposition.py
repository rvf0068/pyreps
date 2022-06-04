import sympy as sp

from pyreps.matrix import MatrixRepresentation


def block(M):
    """A function that returns the indices of the columns where the blocks of a
    matrix end.

    Args:
        M (Matrix): A sympy (square) matrix

    Returns:
        v (list): A list that shows where the blocks
            of a matrix end.

    Examples:
        >>> from sympy.matrices import Matrix
        >>> from pyreps.decomposition import block
        >>> M = Matrix([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 1, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1]])
        >>> block(M)
        [0, 4, 5, 6, 7, 8, 9]

    """
    v = []
    c1 = 0
    i = 0
    n = M.shape[0]
    while (c1 < n):
        c = 0
        for j in range(c1, n):
            if (M[i, j] != 0 or M[j, i] != 0):
                if (sp.Abs(i-j) > c):
                    c = sp.Abs(i-j)
                    # c is (size of the first block) - 1
        if (c == 0):
            v.append(c1)
            c1 = c1+1
            i = c1
        else:
            bloques = False
            while not bloques:
                bloques = True
                for j in range(c1, c1+c+1):
                    for k in range(c1+c+1, n):
                        if (M[j, k] != 0 or M[k, j] != 0):
                            if (sp.Abs(i-k) > c):
                                c = sp.Abs(i-k)
            v.append(c1+c)
            c1 = c1+c+1
            i = c1
    return v


def blockI(M, n, i):
    """A function that inserts the matrix `M`, starting at the entry `(i, i)` of a
    identity matrix of degree `n`.

    Args:
        M (Matrix): The matrix which will be put in the entry `(i, i)` of
            a identity matrix of degree n.

        n (int): Size of the resulting matrix.

        i (int): A integer that will indicated the entry `(i, i)` of the
            identity matrix.
    Returns:
        N (Matrix): The identity matrix which contains the matrix `M`
            in the entry `(i, i)`.

    Raises:
        IndexError: If the number of the columns or the rows plus i are bigger
            to n.

    Examples:
        To use this function use ``blockI(M, n, i)``.

        >>> from sympy.matrices import Matrix
        >>> from pyreps.decomposition import blockI
        >>> M = Matrix([[1, 1, 1], [1, 1, 1]])
        >>> blockI(M, 4, 0)
        Matrix([
        [1, 1, 1, 0],
        [1, 1, 1, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]])
    
    """
#    a=M.shape[0]
    N = sp.eye(n)
    for j in range(0, M.shape[0]):
        for k in range(0, M.shape[1]):
            N[j+i, k+i] = M[j, k]
    return N


def MTS(A):
    """Create a non singular upper triangular matrix V such that V*AV=I.

    Args:
        A (Matrix): A positive definite Hermitian matrix.

    Returns:
        V (Matrix): A non singular upper triangular matrix V that
            V*AV=I.

    Examples:
        >>> from sympy.matrices import Matrix
        >>> M = Matrix([[1, 0, 1], [2, -1, 3], [4, 3, 2]])
        >>> A = M.H*M

        .. Note:: A is positive definite Hermitian matrix.

        >>> from pyreps.decomposition import MTS
        >>> V = MTS(A)
        >>> V
        Matrix([
        [sqrt(21)/21, -sqrt(2310)/231, -12*sqrt(110)/11],
        [          0,  sqrt(2310)/110, 87*sqrt(110)/110],
        [          0,               0,        sqrt(110)]])
        >>> V.H*A*V
        Matrix([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]])

    """
    A1 = A
    n = A.shape[0]
    V = sp.eye(n)
    for i in range(n):
        C = sp.eye(n)
        C[i, i] = 1/sp.sqrt(A1[i, i])
        for j in range(i+1, n):
            C[i,  j] = -(1/A1[i, i])*A1[i, j]
        V = V*C
        A1 = (C.H)*A1*C
    return V


def unitary_representation(d):
    """Becomes a matrix representation in a unitary matrix representation.

    Args:
        d (pyreps.matrix.MatrixRepresentation): A pyreps representation

    Returns:
        __main__.MatrixRepresentation: A unitary matrix representation.

    """
    G = d.group
    n = d.degree
    A = sp.zeros(n, n)
    for g in d._elements:
        J = (d.map[g].H)*d.map[g]
        J = sp.expand(J)
        A = J+A
    V = MTS(A)
    M = {g: sp.ImmutableMatrix((V.inv())*d.map[g]*V) for g in d._elements}
    return MatrixRepresentation(M, G, n)


def is_irreducible(d):
    """Determines if a representation is irreducible.

    Args:
        d (pyreps.matrix.MatrixRepresentation): A pyreps representation

    Returns:
        True if the representation is irreducible,  otherwise, a non scalar
    matrix that reduces the matrix representation.

    Examples:
        To see if a representation is irreducible use
        ``is_irreducible(G, d)``.

        # >>> from sympy.combinatorics.named_groups import SymmetricGroup
        # >>> from pyreps.matrix import regular_representation
        # >>> from pyreps.decomposition import is_irreducible
        # >>> G = SymmetricGroup(3)
        # >>> rr = regular_representation(G)
        # >>> M = is_irreducible(G, rr)
        # >>> M
        # Matrix([
        # [  0, 1/3,   0,   0,   0,   0],
        # [1/3,   0,   0,   0,   0,   0],
        # [  0,   0,   0, 1/3,   0,   0],
        # [  0,   0, 1/3,   0,   0,   0],
        # [  0,   0,   0,   0,   0, 1/3],
        # [  0,   0,   0,   0, 1/3,   0]])
    """
    n = d.degree
    N = sp.eye(n)
    R = unitary_representation(d)
    for r in range(n):
        for s in range(n):
            H = sp.zeros(n)
            if (n-1-r == n-1-s):
                H[n-1-r, n-1-r] = 1
            else:
                if (n-1-r > n-1-s):
                    H[n-1-r, n-1-s] = 1
                    H[n-1-s, n-1-r] = 1
                else:
                    H[n-1-r, n-1-s] = 1*sp.I
                    H[n-1-s, n-1-r] = -1*sp.I
            M = sp.zeros(n)
            for g in R._elements:
                M = M+(R.map[g].H*H*R.map[g])
            M = (sp.sympify(1)/n)*M
            M = sp.expand(M)
            if (M != M[0, 0]*N):
                return M
    else:
        return True


def reduce(d):
    """Decompose a representation into irreducibles

    Args:
        d (dict): The representation.

    Returns:
        U (Matrix): A matrix which decomposes the representation.
    """
    G = d.group
    b = d.degree
    M = is_irreducible(d)
    if M is True:
        return(sp.eye(b))
    else:
        (P,  J) = M.jordan_form()
        P = sp.expand(P)
        w = [block(P.inv()*d.map[g]*P) for g in d._elements]
        length = len(w[0])
        au = w[0]
        for g in w:
            if (len(g) < length):
                length = len(g)
                au = g
        e = 0
        U = P
        for a in au:
            d1 = {g: sp.ImmutableMatrix((P.inv()*d.map[g]*P)[e:a+1, e:a+1])
                  for g in d._elements}
            U = U*blockI(reduce(MatrixRepresentation(d1, G, (a+1-e))), b, e)
            e = a+1
        return sp.ImmutableMatrix(U)
