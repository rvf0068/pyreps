import sympy as sp


class MatrixRepresentation:
    """A class of matricial representation of a group.

    ...

    Attributes:

        d (dict): A dict that must contains the map of a group into the group
            of nonsingular linear transformations of some finite dimensional
            vector space.
        G (sympy.combinatorics.perm_groups.PermutationGroup): The group.
            ..Note:: For our purposes we will work with this class of groups.
        n (int): The degree of the group.

    """
    def __init__(self, d, G, n):
        """Define a matrix representation.

        Args:
            d (dict): A dict that must contains the map of a group into the
                group of nonsingular linear transformations of some finite
                dimensional vector space.
            G (sympy.combinatorics.perm_groups.PermutationGroup): The group.

            n (int): The degree of the group.
        """
        self.map = d
        self.group = G
        self.degree = n

    def character(self):
        """Returns the character of a representation for every element in the group.

        Returns
            dict: A dictionary with the character of the matrix representation
            for every element in the group.

        """
        return dict([(g, self.map[g].trace()) for g in self.group.elements])

    def is_unitary(self):
        """Returns if the matrix representation is unitary.

        Returns
            bool: True if the matrix representation is unitary, False
            otherwise.

        Examples:
            To see if the representation is unitary use
            ``MatrixRepresentation.is_unitary()``, in this case
            we will help us of the ``regular representation(G)``.

            >>> from sympy.combinatorics.named_groups import SymmetricGroup
            >>> from pyreps.matrix import regular_representation
            >>> rr=regular_representation(SymmetricGroup(3))
            >>> print(rr.is_unitary())
            True

        """
        for g in self.group.elements:
            if sp.expand(self.map[g].H*self.map[g]) != sp.eye(self.degree):
                return False
        else:
            return True


def _char_f(G, g, i, j):
    elems = list(G.elements)
    if elems[i]*g == elems[j]:
        return 1
    else:
        return 0


def regular_representation(G):
    """Builds the regular representation.

    Args:
        G (sympy.combinatorics.perm_groups.PermutationGroup): A symmetric
        group.

    Returns:
        __main__.MatrixRepresentation: The matrix regular representation.

    """
    elems = list(G.elements)
    n = len(elems)
    mydict = {}
    for g in elems:
        mydict[g] = sp.ImmutableMatrix(sp.Matrix(n, n,
                                                 lambda i, j:
                                                 _char_f(G, g, i, j)))
    return MatrixRepresentation(mydict, G, n)
