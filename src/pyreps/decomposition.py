import sympy as sp

from sympy import eye, sqrt, I, sympify, expand, Abs
from sympy.matrices import zeros

from pyreps.matrix import MatrixRepresentation


def block(M):
    """A function that return where end the blocks of a matrix.

    Args:
        M (Matrix): The matrix which will be find their blocks.

    Returns:
        v (list): A list that indicated where end the blocks
            of a matrix.

    Examples:
        To find the blocks of a matrix use ``block(M)``.

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
                if (Abs(i-j) > c):
                    c = Abs(i-j)
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
                            if (Abs(i-k) > c):
                                c = Abs(i-k)
            v.append(c1+c)
            c1 = c1+c+1
            i = c1
    return v


def blockI(M, n, i):
    """A function that given a matrix, put it since the entry (i,i) of a
    identity matrix of degree n.

    Args:
        M (Matrix): The matrix which will be put in the entry (i,i) of
            a identity matrix of degree n.

        n (int): Determine the size of the identity matrix.

        i (int): A integer that will indicated the entry (i,i) of the
            identity matrix.
    Returns:
        N (Matrix): The identity matrix with contains the matrix M
            in the entry (i,i).

    Raises:
        IndexError: If the number of the columns or the raws plus i are bigger
            to n.

    Examples:
        To use this function use ``blockI(M, n, i)``.

        >>> from sympy.matrices import Matrix
        >>> from pyreps.decomposition import blockI
        >>> M = Matrix([[1, 1, 1], [1, 1, 1]])
        >>> blockI(M,4,0)
        Matrix([
        [1, 1, 1, 0],
        [1, 1, 1, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]])
    
    """
#    a=M.shape[0]
    N = eye(n)
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
        To create this matrix use ``MTS(A)``.

        >>> from sympy.matrices import Matrix
        >>> M = Matrix([[1, 0, 1], [2, -1, 3], [4, 3, 2]])
        >>> N = M.H
        >>> A = N*M

        .. Note:: A is positive definite Hermitian matrix.

        >>> from pyreps.decomposition import MTS
        >>> V = MTS(N*M)
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
    V = eye(n)
    for i in range(0, n):
        C = eye(n)
        C[i, i] = 1/sqrt(A1[i, i])
        for j in range(i+1, n):
            C[i,  j] = -(1/A1[i, i])*A1[i, j]
        V = V*C
        V.simplify()
        A1 = (C.H)*A1*C
        A1.simplify()
    return V


def unitary_representation(G, d):
    """Becomes a matrix representation in a unitary matrix representation.

    Args:
        d (dict): A dict that must contains the map of a group into the group
            of nonsingular linear transformations of some finite dimensional
            vector space.
        G (sympy.combinatorics.perm_groups.PermutationGroup): The group.

    Returns:
        __main__.MatrixRepresentation: A unitary matrix representation.

    """
    n = d.degree
    A = zeros(n, n)
    for g in d.map:
        J = (d.map[g].H)*d.map[g]
        J = sp.expand(J)
        A = J+A
    A1 = A
    V = eye(n)
    for i in range(0, n):
        C = eye(n)
        C[i, i] = 1/sqrt(A1[i, i])
        for j in range(i+1, n):
            C[i, j] = -(1/A1[i, i])*A1[i, j]
        V = V*C
        V = sp.expand(V)
        A1 = (C.H)*A1*C
        A1 = sp.expand(A1)
    V = MTS(A)
    M = {}
    for g in list(G.generate()):
        M[g] = sp.ImmutableMatrix((V.inv())*d.map[g]*V)
    return MatrixRepresentation(M, G, n)


def is_irreducible(G, d):
    """Determines if a representation is irreducible.

    Args:
        d (dict): A dict that must contains the map of a group into the group
            of nonsingular linear transformations of some finite dimensional
            vector space.

        G (sympy.combinatorics.perm_groups.PermutationGroup): The group.

    Returns:
        True if the representation is irreducible,  a matrix non
        escalar that reduce the matrix representation in otherwise.

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
    N = eye(n)
    for r in range(0, n):
        for s in range(0, n):
            H = zeros(n)
            if (n-1-r == n-1-s):
                H[n-1-r, n-1-r] = 1
            else:
                if (n-1-r > n-1-s):
                    H[n-1-r, n-1-s] = 1
                    H[n-1-s, n-1-r] = 1
                else:
                    H[n-1-r, n-1-s] = 1*I
                    H[n-1-s, n-1-r] = -1*I
            M = zeros(n, n)
            R = unitary_representation(G, d)
            for g in R.map:
                M = M+(R.map[g].H*H*R.map[g])
            M = (sympify(1)/n)*M
            M = expand(M)
            if (M != M[0, 0]*N):
                return M
    else:
        return True


def reduce(G, d):
    """Decompose a representation into irreducibles

    Args:
        G (Group): The group.

        d (dict): The representation.

    Returns:
        U (Matrix): A matrix wich descompose the representation.
    """
    M = is_irreducible(G, d)
    b = d.degree
    if M:
        return(eye(b))
    else:
        (P,  J) = M.jordan_form()
        P = expand(P)
        w = []
        for g in d.map:
            w.append(block(P.inv()*d.map[g]*P))
        length = len(w[0])
        au = w[0]
        for g in w:
            if (len(g) < length):
                length = len(g)
                au = g
        e = 0
        U = P
        for a in au:
            d1 = {}
            for g in list(G.generate()):
                d1[g] = sp.ImmutableMatrix((P.inv()*d.map[g]*P)[e:a+1, e:a+1])
            U = U*blockI(reduce(G, MatrixRepresentation(d1, G, (a+1-e))), b, e)
            e = a+1
        return U
