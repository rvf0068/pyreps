import math
import numpy as np

from sympy import eye
from sympy.matrices import zeros
from sympy.combinatorics import Permutation
from sympy.combinatorics.partitions import IntegerPartition
from sympy.core import S

from pyreps.pchains import P_chains
from pyreps.young import YoungTableaux


def tuple_sorted(a):
    """Sorted tuples of tuples.
    Args:
        a (tuple): A tuple the which will be sorted.

    Returns:
        (tuple): The tuple sorted.

    Examples:
        The function ``sorted`` does not sort tuples of tuples, but
        this function can do it.

        >>> from pyreps.utilities import tuple_sorted
        >>> a1 = ((6, 5), (1, 0), (3, 2))
        >>> a2 = (((4, 2), (1, 0), (5, 3)), ((2, 3), (1, 0), (6, 4)))
        >>> sorted(a1)
        [(1, 0), (3, 2), (6, 5)]
        >>> tuple_sorted(a1)
        ((0, 1), (2, 3), (5, 6))
        >>> sorted(a2)
        [((2, 3), (1, 0), (6, 4)), ((4, 2), (1, 0), (5, 3))]
        >>> tuple_sorted(a2)
        (((0, 1), (2, 3), (4, 6)), ((0, 1), (2, 4), (3, 5)))

    """
    if isinstance(a, int):
        return a
    if isinstance(a[0], int):
        return sorted(a)
    else:
        w = []
        for b in a:
            w.append(tuple(tuple_sorted(b)))
        return tuple(sorted(tuple(w)))


def tuple_permutation(v, P):
    """Determines the orientation of ``b`` taken the orientation of ``a`` positive.

    Args:
        a (tuple): The tuple which will under the Permutation ``p``.
        p (<class 'sympy.combinatorics.permutations.Permutation'>): The
        Permutation.
    Returns:
        (tuple): The tuple with their elements permutated under
        the permutation ``p``.

    Examples:
        To do act the Permutation on the tuple use
        ``tuple_permutation(tuple)``.

        >>> from pyreps.utilities import tuple_permutation
        >>> a1 = (0, 1, 2, 3, 4)
        >>> a2 = ((2, 4), (1, 5), (3, 0))
        >>> a3 = (((0, 1), (2, 4), (3, 5)), ((0, 5), (1, 3), (2, 4)))
        >>> tuple_permutation(a1, Permutation(0, 1, 2))
        (1, 2, 0, 3, 4)
        >>> tuple_permutation(a2, Permutation(1, 3))
        ((2, 4), (3, 5), (1, 0))
        >>> tuple_permutation(a3, Permutation(0, 1)(2, 3))
        (((1, 0), (3, 4), (2, 5)), ((1, 5), (0, 2), (3, 4)))

        .. Note:: The function return other tuple that represent
                  how the Permutation is acting in a natural way in the origin
                  tuple.

    """
    u = []
    w = list(v).copy()
    test = True
    for i in range(len(v)):
        if isinstance(v[i], int):
            if (v[i] in P):
                w[i] = P(v[i])
        else:
            u.append(tuple_permutation(tuple(v[i]), P))
            test = False
    if test:
        return tuple(w)
    else:
        return tuple(u)


def eq_elements(a, b):
    """A function that identify when tuples are equal except by orientation.

    Args:
        a (tuple): The first tuple.
        b (tuple): The second tuple.

    Returns:
        bool: True if the tuples are equal except by orientation,
        False otherwise.

    Raises:
        TypeError: If the tuples don't have the same structure, for
        example:
        a = ((0, 1), (2, 3), (5, 6))
        b = ((1, 0))

    Examples:
        To see if two tuples are equal use ``eq_elements``.

        >>> from pyreps.utilities import eq_elements
        >>> a1 = ((0, 1), (2, 3), (5, 6))
        >>> b1 = ((0, 3), (2, 1), (5, 6))
        >>> a2 = ((0, 1), (2, 3), (5, 6))
        >>> b2 = ((6, 5), (1, 0), (3, 2))
        >>> eq_elements(a1, b1)
        False
        >>> eq_elements(a2, b2)
        True

    """
    if isinstance(a, int):
        return a == b
    if isinstance(a[0], int):
        return (set() == set(a).difference(set(b)))
    else:
        for i in range(len(a)):
            test = False
            for j in range(len(b)):
                if eq_elements(a[i], b[j]):
                    test = True
            if not test:
                return False
        else:
            return True


def orientation_function(a, b, p):
    """Determines the orientation of ``b`` taken the orientation of ``a`` positive.

    Args:
        a (tuple): The first tuple.
        b (tuple): The second tuple.
        p (tuple): The dimension of the simplex.

    Returns:
        ValueError: If the tuples are not equal under the
        function ``eq_elements``.

    Examples:
        To see if two tuples are equal use ``eq_elements``.

        >>> from pyreps.utilities import orientation_function
        >>> a1 = (((0, 1), (2, 3), (4, 6)), ((0, 1), (2, 4), (3, 5)))
        >>> b1 = (((4, 2), (1, 0), (5, 3)), ((2, 3), (1, 0), (6, 4)))
        >>> a2 = ((0, 1), (2, 3), (5, 6))
        >>> b2 = ((6, 5), (1, 0), (3, 2))
        >>> orientation_function(a1, b1, 1)
        False
        >>> orientation_function(a2, b2, 2)
        True

        .. Note:: For ``a1`` and ``b1`` the function receive
                  the integer ``1``, and in the case of ``a2`` and ``b2``
                  receive ``2`` to indentify the dimension of the
                  p-simplex, in the practice the class determine this
                  number.

    """
    if (p == 0):
        return True
    else:
        v = np.zeros((len(a),), dtype=int)
        for i in range(len(a)):
            for j in range(len(b)):
                if eq_elements(a[i], b[j]):
                    v[j] = i
        P = Permutation(v)
        return P.is_even


def boundary_op_n(v):
    """Returns the action of the boundary operator on p-chains.

    Args:
        v ( __main__.P_chains): The p-chain of

    Returns:
        v ( __main__.P_chains): The (p-1)-chain lower the boundary operator.

    Examples:
        To use of the boundary operator, write ``boundary_op_n``.

        >>> from pyreps.pchains import P_chains
        >>> from pyreps.utilities import boundary_op_n
        >>> w = P_chains([(0,), (1,), (2,), (3,)], [1, 1, 1, 1])
        >>> v = P_chains([(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)], [1, 1, 1, 1])
        >>> u = boundary_op_n(v)
        >>> boundary_op_n(w).dic
        {}
        >>> u.dic
        {(1, 2): 2, (0, 2): 0, (0, 1): 2, (1, 3): 0, (0, 3): -2, (2, 3): 2}
        >>> boundary_op_n(u).dic
        {(2,): 0, (1,): 0, (0,): 0, (3,): 0}

        .. Note:: Above w, v are the 0-simplex, 2-simplex of the tetrahedron
                  respectively, if v is a p-simplex we denoted the boundary_op_n(v) =
                  partial_{p}(v), the theory
                  said that partial_{p-1}(partial_{p}(v)) = 0, that was checked in
                  the ``boundary_op_n(u)``. In the case when you use 0-simplex, the
                  result is a empty dictionary like in ``boundary_op_n(w)``.

   """
    p = len(list(v.dic.keys())[0]) - 1
    s = P_chains([], [])
    if (p != 0):
        for u in v.dic.keys():
            c = 0
            for i in u:
                w = list(u).copy()
                w.remove(i)
                if orientation_function(tuple(tuple_sorted(tuple(w))),
                                        tuple(w), p):
                    s1 = P_chains([tuple(tuple_sorted(tuple(w)))],
                                  [abs(v.dic[u])])
                    if (np.sign((v.dic[u])*(-1)**c) < 0):
                        s = s - s1
                    else:
                        s = s + s1
                    c = c+1
                else:
                    s1 = P_chains([tuple(tuple_sorted(tuple(w)))],
                                  [abs(v.dic[u])])
                    if (np.sign((v.dic[u])*(-1)**(c+1)) < 0):
                        s = s - s1
                    else:
                        s = s + s1
                    c = c+1
        return s
    else:
        return s


def nullspace(A):
    """Returns a  ``list`` of column vectors that span the nullspace of the matrix.
    Args:
        A (Matrix): The matrix which we will find the nullspace.
        p (<class 'sympy.combinatorics.permutations.Permutation'>): The
        Permutation.

    Returns:
        (list): A list of list with the generators of the kernel.

    Examples:
        To find the nullspace of a matrix, use ``nullspace(A)``.

        >>> from sympy.matrices import Matrix
        >>> from pyreps.utilities import nullspace
        >>> M1 = Matrix([[2, 4, 6, 6], [8, 20, 0, 1], [5, 0, 3, 2]])
        >>> M2 = Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, -1], [0, -1, 0], [-1, 0, 0]])
        >>> nullspace(M1)
        [[3/16, -1/8, -47/48, 1]]
        >>> nullspace(M2)
        [array([0, 0, 0])]

        .. Note:: Essentially the function only obtain the nullspace
                  with the function ``A.nullspace()`` and returns the trivial kernel
                  if ``A.nullspace()`` is a emtpy list.

    """
    u = A.nullspace()
    w = []
    for g in u:
        v = []
        for i in g:
            v.append(i)
        w.append(v)
    if (w == []):
        return [np.zeros((A.shape[1],), dtype=int)]
    else:
        return w


def columnspace(A):
    """Returns a ``list`` of column vectors that span the columnspace of the
       matrix.

    Args:
        A (Matrix): The matrix which we will find the columnspace.
        p (<class 'sympy.combinatorics.permutations.Permutation'>): The
        Permutation.

    Returns:
        (list): A list of list with the generators of the columnspace (image).

    Examples:
        To find the columnspace of a matrix, use ``columnspace(A)``.

        >>> from sympy.matrices import Matrix
        >>> from pyreps.utilities import nullspace, columnspace
        >>> M1 = Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
        >>> M2 = Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, -1], [0, -1, 0], [-1, 0, 0]])
        >>> columnspace(M1)
        [array([0, 0, 0])]
        >>> nullspace(M2)
        [array([0, 0, 0])]
        >>> columnspace(M2)
        [[1, 0, 0, 0, 0, -1], [0, 1, 0, 0, -1, 0], [0, 0, 1, -1, 0, 0]]

        .. Note:: Essentially the function only obtain the columnspace
                  with the function ``A.columnspace()`` and returns the trivial image
                  if ``A.columnspace()`` is a emtpy list. In the example ``M2`` is noted
                  that is right that:
                  ``dimension(nullspace(A))``+``dimension(columnspace(A))`` = ``number of
                  columns``.

    """
    u = A.columnspace()
    w = []
    for g in u:
        v = []
        for i in g:
            v.append(i)
        w.append(v)
    if (w == []):
        return [np.zeros((A.shape[0],), dtype=int)]
    else:
        return w


def Reduce(N):
    """Returns a row reduced form of a matrix and a matrix that save the
       operations.

    Args:
        N (Matrix): The matrix which will be operated.

    Returns:
        tuple: The first element is the matrix which to be multiplied
        by the right to the original matrix, return the row reduced form
        and the other object is the row reduced form of the origin matrix.

    Examples:
        To use this functio use ``Reduce(Matrix)``. We will use the help of the
        function ``rref`` to verify that the result is right.

        >>> from sympy.matrices import Matrix
        >>> from pyreps.utilities import Reduce
        >>> M = Matrix([[-1, -1, -1, -1, 0, 0, 0, 0], [ 1, 0, 0, 0, -1, -1, 0, 0],[ 0, 1, 0, 0, 1, 0, -1, -1], [ 0, 0, 1, 0, 0, 1, 1, 0], [ 0, 0, 0, 1, 0, 0, 0, 1]])
        >>> M.rref()
        (Matrix([
        [1, 0, 0, 0, -1, -1,  0,  0],
        [0, 1, 0, 0,  1,  0, -1, -1],
        [0, 0, 1, 0,  0,  1,  1,  0],
        [0, 0, 0, 1,  0,  0,  0,  1],
        [0, 0, 0, 0,  0,  0,  0,  0]]), (0, 1, 2, 3))
        >>> P = Reduce(M)
        >>> P
        (Matrix([
        [ 0,  1,  0,  0, 0],
        [ 0,  0,  1,  0, 0],
        [ 0,  0,  0,  1, 0],
        [-1, -1, -1, -1, 0],
        [ 1,  1,  1,  1, 1]]), Matrix([
        [1, 0, 0, 0, -1, -1,  0,  0],
        [0, 1, 0, 0,  1,  0, -1, -1],
        [0, 0, 1, 0,  0,  1,  1,  0],
        [0, 0, 0, 1,  0,  0,  0,  1],
        [0, 0, 0, 0,  0,  0,  0,  0]]))

        ..Note:: The first matrix is the row reduced form, and the second
                 is a matrix which if is multiplied the left size to the
                 origin matrix, then we obtain the row reduced form, like below.
                 print(P[0]*M == (M.rref())[0])
                 True

    """
    M = N.copy()
    lead = 0
    rowCount = M.shape[0]
    columnCount = M.shape[1]
    A = eye(rowCount)
    # v=[]
    # v.append(A)
    for r in range(rowCount):
        # display(r)
        B1 = eye(rowCount)
        # B2=eye(rowCount)
        # B3=eye(rowCount)
        if (columnCount <= lead):
            return A, M
        i = r
        while (M[i, lead] == 0):
            i = i + 1
            if (rowCount == i):
                i = r
                lead = lead + 1
                if (columnCount == lead):
                    return A, M
        B1.row_swap(i, r)
        M.row_swap(i, r)
        a = M[r, lead]
        for k in range(columnCount):
            M[r, k] = S(M[r, k])/a
            if (k < rowCount):
                B1[r, k] = S(B1[r, k])/a
        for i in range(0, rowCount):
            if (i != r):
                a = M[i, lead]
                for k in range(0, columnCount):
                    M[i, k] = M[i, k] - M[r, k]*a
                    if (k < rowCount):
                        B1[i, k] = B1[i, k] - B1[r, k]*a
        lead = lead + 1
        A = B1*A
    return A, M


def permutation_in_simplex_test(vec, P):
    r"""Returns a simplex under a permutation.

    Args:
        vec ( __main__.P_chains): A p-chain which the permutation will act.
        P ( sympy.combinatorics.permutations.Permutation): The permutation.

    Returns:
        (__main__.P_chains): A new p-chain that is the result of the
        permutation acting on the original p-chain ``vec``.

    Examples:
        To see how a permutation act on a p-simplex, use
        ``permutation_in_simplex_test(SimplicialComplex, Permutation)``.
        Also we must check that the boundary operator on a p-simplex
        (partial_{p}) is well-defined and that (if ``p-simplex`` := sigma)
        partial_{p}(-sigma) = - partial_{p}(sigma). For this purpose, it
        suffices to show that the right-hand side of:

        .. math::

                  \partial_{p}(\sigma) = \partial_{p}([v_{0},...,v_{p}]) =
                  = \sum_{i=0}^{p}(-1)^{i}[v_{0},...,v_{i},...
                  v_{p}].

        (where v_{i} means that the vertex v_{i} is to be deleted
        from the array)

        changes sign if we exchange two adjacent vertices in the array
        [v_{0},...,v_{p}] (important step will be explain according with the
        theory):

        >>> from pyreps.pchains import P_chains
        >>> from pyreps.utilities import boundary_op_n
        >>> u1 = P_chains([(0, 1, 2, 3)], [1])
        >>> u2 = P_chains([(0, 2, 1, 3)], [1])
        >>> bu1 = boundary_op_n(u1).dic
        >>> bu2 = boundary_op_n(u2).dic
        >>> bu1
        {(1, 2, 3): 1, (0, 2, 3): -1, (0, 1, 3): 1, (0, 1, 2): -1}
        >>> bu2
        {(1, 2, 3): -1, (0, 1, 3): -1, (0, 2, 3): 1, (0, 1, 2): 1}

        .. Note:: The p-simplex in u1 and u2 differ by a sign. And we could
                  see that the result changes sign, like is wanted.

        Now se must check that partial_{p}(rho(sigma)) =
        rho(partial_{p}(sigma))
        (where rho = Permutation ``P``). For this we will use some p-simplices
        associated with a graph.

        >>> from pyreps.simplicial import SimplicialComplex
        >>> from pyreps.utilities import boundary_op_n, permutation_in_simplex_test
        >>> G = nx.complete_graph(5)
        >>> sc = SimplicialComplex(G)
        >>> sigma = sc.basis_group_oriented_p_chains(1)
        >>> sigma.dic
        {(0, 1): 1, (0, 2): 1, (0, 3): 1, (0, 4): 1, (1, 2): 1, (1, 3): 1, (1, 4): 1, (2, 3): 1, (2, 4): 1, (3, 4): 1}
        >>> bo_sigma = boundary_op_n(sigma)
        >>> rho_bo_sigma = permutation_in_simplex_test(bo_sigma, Permutation(0, 1))
        >>> rho_bo_sigma.dic
        {(0,): -2, (1,): -4, (2,): 0, (3,): 2, (4,): 4}
        >>> rho_sigma = permutation_in_simplex_test(sigma, Permutation(0, 1))
        >>> bo_rho_sigma=boundary_op_n(rho_sigma)
        >>> bo_rho_sigma.dic
        {(1,): -4, (0,): -2, (2,): 0, (3,): 2, (4,): 4}
        >>> rho_bo_sigma == bo_rho_sigma
        True

        .. Note:: Then for this example the result is the same.

        And for the second property:

        >>> from pyreps.pchains import P_chains
        >>> from pyreps.utilities import boundary_op_n
        >>> sigma1 = P_chains([(0, 1, 2)], [1])
        >>> sigma2 = P_chains([(0, 1, 2)], [-1])
        >>> w1 = boundary_op_n(sigma1)
        >>> w2 = boundary_op_n(sigma2)
        >>> w1
        {(1, 2): 1, (0, 2): -1, (0, 1): 1}
        >>> w2
        {(1, 2): -1, (0, 2): 1, (0, 1): -1}
        >>> w1 == (-1)*w2 #Multiply by -1.
        True

        .. Note:: The simplices differ by the sign, and for this example is true
                  that partial_{p}(-sigma) = - partial_{p}(sigma)
                  like is wanted, and for all our cases the previous is true.

    """
    s = P_chains([], [])
    if (vec.dic != {}):
        v = list(vec.dic.keys())
        p = len(list(vec.dic.keys())[0]) - 1
        # ve = vec
        for a in v:
            if isinstance(a, int):
                return vec
            else:
                w = tuple_permutation(a, P)
                if orientation_function(tuple_sorted(w), w, p):
                    s = s + P_chains([tuple(tuple_sorted(w))], [vec.dic[a]])
                else:
                    s = s - P_chains([tuple(tuple_sorted(w))], [vec.dic[a]])
        return s
    else:
        return s


def partitions_list(n):
    """Returns a list of the partitions of n.

    Args:
        n (int): A integer that determine the partitions.

    Returns:
        w (list): A list of list that are the different partitions of n.

    Raises:
        ValueError: If ``n`` is not bigger than zero.

    Examples:
        To form all the partitions for the integer ``n``, use
        ``list_partitions``.

        >>> from pyreps.utilities import partitions_list
        >>> partitions_list(3)
        [[3], [1, 1, 1], [2, 1]]
        >>> partitions_list(4)
        [[4], [1, 1, 1, 1], [2, 1, 1], [2, 2], [3, 1]]


    """
    p = IntegerPartition([n])
    w = []
    while list(p.args[1]) not in w:
        w.append(list(p.args[1]))
        p = p.next_lex()
    return w


def form_matrix_yt(w):
    """Returns the a matrix that represent the character table of the symmetric
       group.

    Args:
        w (list): A list with the partitions for certain symmetric group.

    Returns:
        M (<class 'sympy.matrices.dense.MutableDenseMatrix'>): A matrix with
        the characters
        of the character table of the symmetric group.

    Examples:
        To form the matrix that represent the character table of the
        symmetric group, use ``form_matrix_yt``.

        >>> from pyreps.utilities import partitions_list, form_matrix_yt
        >>> v = partitions_list(3)
        >>> form_matrix_yt(v)
        Matrix([
        [ 1, 1,  1],
        [ 1, 1, -1],
        [-1, 2,  0]])

        .. Note:: The function need a list of the partitions for
                  ``n``, then is used the function ``partitions_list.``

    """
    M = zeros(len(w), len(w))
    for i in range(len(w)):
        for j in range(len(w)):
            M[i, j] = YoungTableaux(w[i], w[j]).CMNR()
    return M


def make_permutation(partition):
    """Given a partition returns the a representate of a conjugacy class.

    Args:
        partition (list): Represents the partitions of a symmetric group (n).

    Returns:
        sympy.combinatorics.permutations.Permutation: A representate of the
        conjugacy class.

    Examples:
        For find a representate of a conjugacy class of s simmetric group
        use ``make_permutation(partition)``.

        >>> from pyreps.utilities import make_permutation
        >>> make_permutation([5])
        Permutation(0, 1, 2, 3, 4)
        >>> make_permutation([1, 1, 1, 1, 1])
        Permutation()
        >>> make_permutation([2, 1, 1, 1])
        Permutation(4)(0, 1)
        >>> make_permutation([2, 2, 1])
        Permutation(4)(0, 1)(2, 3)
        >>> make_permutation([3, 1, 1])
        Permutation(4)(0, 1, 2)
        >>> make_permutation([3, 2])
        Permutation(0, 1, 2)(3, 4)
        >>> make_permutation([4, 1])
        Permutation(4)(0, 1, 2, 3)

    """
    P = Permutation()
    c = 0
    for j in range(len(partition)):
        a = []
        for h in range(partition[j]):
            a.append(c)
            c = c + 1
        if (c == 1):
            P1 = Permutation()
            c = 0
        else:
            P1 = Permutation([a])
        P = P*P1
    return P


def size_conjugacy_class(partition, n):
    """Returns the number of elements of a conjugacy class.

    Args:
        partition (list): Represents the partitions of a symmetric group (n).
        n (int): A integer to identify which the symmetric group.

    Returns:
        int: The number of elements of the conjugacy class.

    Examples:
        For find the number of elements of the conjugacy class of a
        symmetric group use ``size_conjugacy_class(partition,n)``.

        >>> from pyreps.utilities import size_conjugacy_class
        >>> size_conjugacy_class([4], 4)
        6
        >>> size_conjugacy_class([1,1,1,1], 4)
        1
        >>> size_conjugacy_class([2,1,1], 4)
        6
        >>> size_conjugacy_class([2,2], 4)
        3
        >>> size_conjugacy_class([3,1], 4)
        8

        .. Note:: The examples showed are all the partition for the case
                  4, and the sum of the results is 24 that is the cardinality of the
                  simmetric group of 4.

    """
    aux1 = 1
    c = 0
    aux = partition[0]
    flag = 1
    # c2 = 0
    for j in range(len(partition)):
        if (aux == partition[j]):
            c = c + 1
            flag = 1
        else:
            aux1 = aux1*(partition[j-1]**c)*(math.factorial(c))
            aux = partition[j]
            c = 1
            flag = 0
            # c2 = j
    if (flag == 1):
        aux1 = aux1*(partition[j-1]**c)*(math.factorial(c))
    else:
        aux1 = aux1*(partition[j]**c)*(math.factorial(c))
    card = (math.factorial(n))/aux1
    return int(card)
