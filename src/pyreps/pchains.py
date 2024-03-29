class P_chains:
    """Class for defining and operating p-chains.

    Attributes:
        keys (list): A list with the elements of the p-chains
        values (list): A list with the coefficient for every element in the
        p-chain

    """

    def __init__(self, keys, values):
        '''Defines a p-chain.

        Args:
            keys (list): A list with the elements of the p-chains.
            values (list): A list with the coefficient for every element in the
            p-chain.

        Raises:
            IndexError: If the number of p-simplexs is not equal to the number
            of coefficients.
            TypeError: If the p-simplex given is not immutable data types like
            a tuple.

        Examples:
            To make a p-chain, use the ``P_chains(list1,list2)`` class.
            A p-chain is constructed by providing a list1 of p-simplex and a
            list2 with their coefficients, i.e.

            >>> P = P_chains([(0, 1, 2, 3)], [1])
            >>> Q = P_chains([(0, 1, 2), (0, 1, 3)], [-1, 2])
            >>> P
            {(0, 1, 2, 3): 1}
            >>> Q
            {(0, 1, 2): -1, (0, 1, 3): 2}

            One important thing to note about P_chains is that the p-simplex
            must be a immutable data type like a tuple.

        '''
        self.keys = keys
        self.values = values
        self.dic = dict(zip(self.keys, self.values))

    def __repr__(self):
        return str(self.dic)

    def __add__(self, other):
        """Sums two p-chains.

        Args:
            other ( __main__.P_chains): Other p-chain.

        Returns
            __main__.P_chains: A new p-chain that is the sum of the two
            p-chains given.

        Examples:
            To sum two p-chains, use ``+``.

            >>> from pyreps.pchains import P_chains
            >>> P = P_chains([(0, 1, 2)], [2])
            >>> Q = P_chains([(0, 1, 3)], [5])
            >>> T = P_chains([(0, 1, 2)], [-2])
            >>> R = P + Q
            >>> L = P + T
            >>> R.dic
            {(0, 1, 2): 2, (0, 1, 3): 5}
            >>> L.dic
            {(0, 1, 2): 0}

            Note that the coefficient of L is zero because the P and T have the
            same simplex, then P + T is the sum of their respective
            coefficients.

        """
        D = {}
        for x in list(self.dic.keys()):
            D[x] = self.dic[x]
        for y in list(other.dic.keys()):
            if y not in list(D.keys()):
                D[y] = other.dic[y]
            else:
                D[y] = D[y] + other.dic[y]
        return P_chains(D.keys(), D.values())

    def __sub__(self, other):
        """Subtracts two p-chains.

        Args:
            other ( __main__.P_chains): Other p-chain.

        Returns
            __main__.P_chains: A new p-chains that is the difference of the two
            p-chains given.

        Examples:
            To subtract two p-chains, use ``-``.

            >>> from pyreps.pchains import P_chains
            >>> P = P_chains([(3,4,5)],[3])
            >>> Q = P_chains([(1,8,9)],[1])
            >>> T = P_chains([(0,1,2,3)],[-7])
            >>> R = P - Q
            >>> L = P - T
            >>> R.dic
            {(3, 4, 5): 3, (1, 8, 9): -1}
            >>> L.dic
            {(3, 4, 5): 3, (0, 1, 2, 3): 7}

        """
        D = {}
        for x in list(self.dic.keys()):
            D[x] = self.dic[x]
        for y in list(other.dic.keys()):
            if y not in list(D.keys()):
                D[y] = -other.dic[y]
            else:
                D[y] = D[y] - other.dic[y]
        return P_chains(D.keys(), D.values())

    def __eq__(self, other):
        '''Returns if the two P_chains are equal

        Args:
            other ( __main__.P_chains): Other p-chain.

        Returns
            bool: The return value. True for success, False otherwise.

        Examples:
            To find if two P_chains are equal use ``==``.

            >>> from pyreps.pchains import P_chains
            >>> P = P_chains([(0, 1, 2), (3, 4, 5)], [1, 1])
            >>> Q = P_chains([(3, 4, 5), (0, 1, 2)], [1, 1])
            >>> R = P_chains([(0, 1, 2)], [1])
            >>> L = P_chains([(1, 0, 2)], [1])
            >>> P == Q
            True
            >>> R == L
            False

            ..Note:: R and L are not equal even though they only are
            distint in orientation, moreover, in this class the
            orientation is not defined yet.

        '''
        return self.dic == other.dic

    def __ne__(self, other):
        '''Returns if the two P_chains are not equal

        Args:
            other ( __main__.P_chains): Other p-chain.

        Returns
            bool: The return value. True for success, False otherwise.

        Examples:
            To know if two P_chains are equal use ``!=``.

            >>> from pyreps.pchains import P_chains
            >>> P = P_chains([(0, 1, 2, 4)], [1])
            >>> Q = P_chains([(0, 1, 4, 2)], [-2])
            >>> R = P_chains([(5, 6, 7)], [1])
            >>> L = P_chains([(5, 6, 7)], [1])
            >>> P != Q
            True
            >>> R != L
            False

        '''
        return not self.__eq__(other)

    def __rmul__(self, other):
        '''Product for a scalar
        Args:
            other ( __main__.P_chains): Other p-chain.

        Returns
            __main__.P_chains: A new p-chain that is the product of the scalar by the 
            p-chain given.

        Examples:
            To multiply a scalar by a P_chain use ``*``.

            >>> from pyreps.pchains import P_chains
            >>> P = P_chains([(7, 8, 9), (10, 11, 12)], [3, 2])
            >>> Q = P_chains([(0, 1, 4, 2, 6)], [-5])
            >>> 3*P
            {(7, 8, 9): 9, (10, 11, 12): 6}
            >>> (-1)*Q
            {(0, 1, 4, 2, 6): 5}

        '''
        aux = [other*self.dic[x] for x in self.dic]
        return P_chains(self.dic.keys(), aux)
