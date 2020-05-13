import sympy as sp


class MatrixRepresentation:
    """Matricial representation of a group.

    Matrix representations are usually constructed using functions
    such as: `regular_representation`.

    Parameters
    ==========

    d : dictionary
        A dictionary that contains the map of the group into the
        group of nonsingular linear transformations of some finite
        dimensional vector space.

    G : A Sympy PermutationGroup
        The group.

    n : int
        The degree of the representation.

    """

    def __init__(self, d, G, n):

        self.map = d
        self.group = G
        self.degree = n

    def character(self):
        """Character of the representation

        Returns the character of a representation, as a dictionary
        with keys given by the elements of the group.

        Examples
        ========

        >>> from sympy.combinatorics.named_groups import SymmetricGroup
        >>> from pyreps.matrix import regular_representation
        >>> rr = regular_representation(SymmetricGroup(3))
        >>> rr.character()
        {Permutation(2): 6, Permutation(0, 1, 2): 0, Permutation(0, 2, 1): 0, Permutation(1, 2): 0, Permutation(2)(0, 1): 0, Permutation(0, 2): 0}

       """
        return {g: self.map[g].trace() for g in self.group.generate()}

    def is_unitary(self):
        """Tests if the matrix representation is unitary.

        Returns ``True`` if the matrix representation is unitary, and
        ``False`` otherwise.

        Examples
        ========

        >>> from sympy.combinatorics.named_groups import SymmetricGroup
        >>> from pyreps.matrix import regular_representation
        >>> rr = regular_representation(SymmetricGroup(3))
        >>> rr.is_unitary()
        True

        """
        for g in self.group.generate():
            if sp.expand(self.map[g].H*self.map[g]) != sp.eye(self.degree):
                return False
        else:
            return True


def _char_f(G, g, i, j):
    elems = list(G.generate())
    if elems[i]*g == elems[j]:
        return 1
    else:
        return 0


def regular_representation(G):
    """Builds the regular representation.

    Returns the regular representation of ``G`` as a ``MatrixRepresentation``.

    Parameters
    ==========

    G : A Sympy PermutationGroup
        The group.

    Examples
    ========

    >>> from sympy.combinatorics.named_groups import SymmetricGroup
    >>> from pyreps.matrix import regular_representation
    >>> type(regular_representation(SymmetricGroup(3)))
    <class 'pyreps.matrix.MatrixRepresentation'>

    """
    n = G.order()
    mydict = {g: sp.ImmutableMatrix(sp.Matrix(n, n,
                                              lambda i, j:
                                              _char_f(G, g, i, j)))
              for g in G.generate()}
    return MatrixRepresentation(mydict, G, n)
