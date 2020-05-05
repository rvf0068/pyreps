from itertools import permutations


class YoungTableaux:
    """A class to compute irreducible character values of a symmetric group.

    ...

    Attributes:
        p_lambda (list): A list that represent the first partition.
        p_rho (list): A list that represent the second partition.

    """

    def __init__(self, p_lambda, p_rho):
        '''
        Args:
            p_lambda (list): A list that represent the first partition.
            p_rho (list): A list that represent the second partition.

        '''
        self.p_lambda = p_lambda
        self.p_rho = p_rho

    def choose_tableaux(self, v):
        """A method to identify if a given list is a border-strip tableaux or not.

        Args:
            v (list): A list of list that is the candidate to be a border-strip
            tableaux.

        Returns:
            True if the list is a border-strip tableaux, False otherwise.

        Raises:
            IndexError: If the list given is not consistent with the partitions
            given.

        Examples:
            To see if a list given is a border-strip tableaux, you must
            use ``YoungTableaux.choose_tableaux(list)``.

            >>> from pyreps.young import YoungTableaux
            >>> YT = YoungTableaux([2,1],[1,1,1])
            >>> YT.choose_tableaux([[1,1,1],[1,2]])
            True
            >>> YT.choose_tableaux([[1,1,1],[2,1]])
            False
            >>> YT.choose_tableaux([[1,1,1],[1,1]])
            False

            .. Note:: The examples given must be interpreted
                      like:
                      Tableaux 1:

                      ::

                                  | 1 | 1 | 1 |
                                  | 1 | 2 |

                      Tableaux 2:

                      ::

                                  | 1 | 1 | 1 |
                                  | 2 | 1 |

                      Tableaux 1:

                      ::

                                  | 1 | 1 | 1 |
                                  | 1 | 1 |

                      In the two first examples we use two partitions of 3,
                      and the third example is done even though the second partition
                      is a partition of 2, for our purposes that mistake is not done,
                      because we give the correct partitions.

        """
        for i in v:
            for j in range(0, len(i)-1):
                if (i[j] > i[j+1]):
                    return False
        for i in range(1, len(v)):
            for j in range(0, len(v[i])):
                if (v[i][j] < v[i-1][j]):
                    return False
        for i in range(0, len(v)):
            for j in range(0, len(v[i])):
                c = 0
                c1 = 0
                if (j != 0):
                    if (v[i][j] == v[i][j-1]):
                        c = 1
                        c1 = c1 + 1
                if (j != (len(v[i])-1)):
                    if (v[i][j] == v[i][j+1]):
                        c = 1
                        c1 = c1 + 1
                if (i != 0):
                    if (v[i][j] == v[i-1][j]):
                        c = 1
                if (i != (len(v)-1)):
                    # if (len(v[i+1]) <= len(v[i])):
                    if (j < len(v[i+1])):
                        if (v[i][j] == v[i+1][j]):
                            c = 1
                            c1 = c1 + 1
                    if (j < (len(v[i+1]) - 1)):
                        if (v[i][j] == v[i+1][j+1]):
                            c1 = c1 + 1
                if ((c == 0) and (self.p_rho[v[i][j]-1] > 1)):
                    return False
                if (c1 == 3):
                    return False
        else:
            return True

    def MNR(self):
        """A method to generate all the border-strip tableauxs according to the
           partitions given.

        Returns:
            D1 (list): A list that contains the border-strip tableauxs.

        Examples:
            To see the border-strip tableauxs, you must
            use ``YoungTableaux.MNR()``.

            >>> from pyreps.young import YoungTableaux
            >>> YT1 = YoungTableaux([2, 2, 2, 1], [3, 3, 1])
            >>> YT2 = YoungTableaux([4, 1], [3, 2])
            >>> YT3 = YoungTableaux([2, 2, 1, 1], [6])
            >>> YT1.MNR()
            [[[1, 1], [1, 2], [2, 2], [3]], [[1, 2], [1, 2], [1, 2], [3]]]
            >>> YT2.MNR()
            [[[1, 1, 2, 2], [1]]]
            >>> YT3.MNR()
            []

            .. Note:: In the two first example the method found two list
                      that are border-strip tableaux, the which could be interpreted
                      like:

                      ::

                          | 1 | 1 |
                          | 1 | 2 |
                          | 2 | 2
                          | 3 |


                      and

                      ::

                          | 1 | 2 |
                          | 1 | 2 |
                          | 1 | 2 |
                          | 3 |.


                      And for the second example the method only found a list,
                      and his parallel border-strip tableaux looks like:

                      ::

                          | 1 | 1 | 2 | 2 |
                          | 1 |

                      And for the third example there are not any border-strip
                      tableaux.

        """
        p = []
        i = 1
        for h in self.p_rho:
            for j in range(0, h):
                p.append(i)
            i = i+1
        perm = permutations(p)
        D = []
        for i in list(perm):
            v = []
            for g in i:
                v.append(g)
            c = 0
            w = []
            for p in self.p_lambda:
                u = []
                for i in range(c, c+p):
                    u.append(v[i])
                w.append(u)
                c = c+p
            if self.choose_tableaux(w):
                D.append(w)
        D1 = []
        if (D != []):
            D1 = [D[0]]
            for k1 in D:
                if k1 not in D1:
                    D1.append(k1)
        return D1

    def heights(self):
        """A method to calculate the heights i.e the sum of the heights of the
           border strips.

        Returns:
            He (list): A list that contains heights.

        Examples:
            To see the heights of the border-strip tableauxs, you must
            use ``heights.MNR()``.

            >>> from pyreps.young import YoungTableaux
            >>> YT1 = YoungTableaux([5, 2, 1], [3, 3, 1, 1])
            >>> YT1.heights()
            [1, 1, 1, 2, 2, 3]

            .. Note:: The border-strip tableauxs generated by these partitions
                      are:
                      T1:

                      ::

                          | 1 | 1 | 1 | 3 | 4 |
                          | 2 | 2 |
                          | 2 |

                      T2:

                      ::

                          | 1 | 1 | 2 | 2 | 2 |
                          | 1 | 3 |
                          | 4 |

                      T3:

                      ::

                          | 1 | 1 | 2 | 2 | 2 |
                          | 1 | 4 |
                          | 3 |

                      T4:

                      ::

                          | 1 | 2 | 2 | 2 | 3 |
                          | 1 | 4 |
                          | 1 |

                      T5:

                      ::

                          | 1 | 2 | 2 | 2 | 4 |
                          | 1 | 3 |
                          | 1 |

                      T6:

                      ::

                          | 1 | 2 | 2 | 3 | 4 |
                          | 1 | 2 |
                          | 1 |.

                      And the last list are their respective heights, which can be
                      verified
                      like follows:
                      ht(T1) = 0 + 1 + 0 + 0 = 1
                      ht(T2) = 1 + 0 + 0 + 0 = 1
                      ht(T3) = 1 + 0 + 0 + 0 = 1
                      ht(T4) = 2 + 0 + 0 + 0 = 2
                      ht(T5) = 2 + 0 + 0 + 0 = 2
                      ht(T6) = 2 + 1 + 0 + 0 = 3

                      And the results above coincide with the list [1, 1, 1, 2, 2, 3].

        """
        H = self.MNR()
        He = []
        for h in H:
            he = []
            for i in range(0, len(self.p_rho)):
                c = 0
                for g in h:
                    if ((i+1) in g):
                        c = c+1
                he.append(c-1)
            He.append(sum(he))
        return He

    def CMNR(self):
        """A method to calculate irreducible character values ​​through the
           Murnaghan-Nakayama rule.

        Returns:
            s (int):  A irreducible character value of a symmetric group
            according to the given partition.

        Examples:
            Here are the last method of the class ``YoungTableaux``. To
            get the irreducible character value of a symmetric group according
            to the given partition use ``YoungTableaux.CMNR`` (We will use
            the same example that in the method ``heights``.

            >>> from pyreps.young import YoungTableaux
            >>> YT1 = YoungTableaux([5, 2, 1], [3, 3, 1, 1])
            >>> YT1.heights()
            [1, 1, 1, 2, 2, 3]
            >>> YT1.CMNR()
            -2

            .. Note:: In the method ``heights`` we saw that for the partitions
                      of 6 given
                      there are six such border-strip tableaux, and their heights are:
                      [1, 1, 1, 2, 2, 3]
          
                      According with the Murnaghan–Nakayama rule the character value is
                      therefore:
                      
                      .. math::

                      Chi_{3,3,1,1}^{5,2,1} = (-1)^(1) + (-1)^(1) + (-1)^(1) + (-1)^(2) +
                                      (-1)^(2) + (-1)^(3) =
                                            = - 1 - 1 - 1 + 1 + 1 - 1 = -2
          
                      Like the method ``CMNR`` got.
          
        """
        He = self.heights()
        s = 0
        for j in He:
            s = s + (-1)**(j)
        return s
