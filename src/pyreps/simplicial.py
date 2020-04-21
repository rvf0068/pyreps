import math
import networkx as nx
import numpy as np

from sympy import eye
from sympy.matrices import Matrix, zeros
from sympy.solvers.solveset import linsolve

from pyreps.pchains import P_chains
from pyreps.utilities import tuple_sorted, tuple_permutation, eq_elements, \
    orientation_function, boundary_op_n, nullspace, columnspace, Reduce, \
    permutation_in_simplex_test, partitions_list, form_matrix_yt, \
    make_permutation, size_conjugacy_class


class SimplicialComplex:
    """A class to make simplicial complex asociated with a graph.

    ...

    Attributes:
        G (networkx.classes.graph.Graph): A graph used to build a simplicial
        complex.

    """

    def __init__(self, G):
        '''Saves the graph and their nodes.

        Args:
            G (networkx.classes.graph.Graph): A graph used to build a
            simplicial complex.

        Raises:
            AttributeError: If G is not a graph.

        Examples:
            To make a simplicial complex asociated with a graph,
            use the ``SimplicialComplex`` class. We need a graph G.

            >>> import networkx as nx
            >>> from pyreps.simplicial import SimplicialComplex
            >>> G = nx.complete_graph(5)
            >>> sc = SimplicialComplex(G)
            >>> sc.vertices
            [0, 1, 2, 3, 4]


        '''
        self.G = G
        self.vertices = []
        for x in self.G.nodes():
            self.vertices.append(x)

    def faces(self):
        """Makes the faces of a simplicial complex.

        A simplicial complex must contains every face of a simplex
        and the intersection of any two simplexes of G is a face of
        each of them.

        Returns:
            list: A list of the faces of a simplex.

        Examples:
            To create the faces of a simplical complex use,
            ``SimplicialComplex.faces()``.

            >>> import networkx as nx
            >>> from pyreps.simplicial import SimplicialComplex
            >>> G = nx.complete_graph(3)
            >>> sc = SimplicialComplex(G)
            >>> sc.faces()
            [(0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2)]

            .. Note:: The faces are sorted by their dimension.

            """

        faceset = []
        for face in list(nx.enumerate_all_cliques(self.G)):
            faceset.append(tuple(face))
        return faceset

    def p_simplex(self, p):
        """Creates a list of the faces of a simplex with dimension p.

        Args:
            p (int): The dimension of the faces.

        Returns:
            list: A list of the faces of a simplex with dimension p.

        Examples:
            The p-simplices are done with
            "SimplicialComplex.p_simplex(p)".

            >>> import networkx as nx
            >>> from pyreps.simplicial import SimplicialComplex
            >>> G = nx.complete_graph(3)
            >>> sc = SimplicialComplex(G)
            >>> sc.faces()
            [(0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2)]
            >>> sc.p_simplex(0)
            [(0,), (1,), (2,)]
            >>> sc.p_simplex(2)
            [(0, 1, 2)]
            >>> sc.p_simplex(1)
            [(0, 1), (0, 2), (1, 2)]
            >>> sc.p_simplex(5)
            []

            .. Note:: If there are not faces of dimension p,
            the method return a empty list like in
            ``sc.p_simplex(5)``.

        """

        return list(filter(lambda face: (len(face) == p+1), self.faces()))

    def dimension(self):
        """Gives the dimension of a simplicial complex.

        Returns:
            a - 1 (int): The dimension of the simplicial complex.

        Raises:
            Return ``-1`` if the graph is empty.

        Examples:
            To use the method dimension write
            ``SimplicialComplex.dimension()``.

            >>> import networkx as nx
            >>> from pyreps.simplicial import SimplicialComplex
            >>> G = nx.petersen_graph()
            >>> sc = SimplicialComplex(G)
            >>> sc.dimension()
            1

        """
        a = 0
        for x in self.faces():
            if (len(x) > a):
                a = len(x)
        return a-1

#    def p_simplex(self, k):
#        p_simplex = []
#        for x in self.p_simplex(k):
#            p_simplex.append(x)
#        return p_simplex
#    def elementary_chain(self, simplex):
#        """Give the p-chains with their repective orientation.
#
#        Args:
#            simplex (tuple): A tuple which orientation is taken positive.
#
#        Returns:
#            __main__.P_chains: A new p-chains that is.
#
#
#        """
#        ec = P_chains([], [])
#        for x in set_oriented_p_simplices(simplex):
#            if (orientation_function(tuple_sorted(simplex), x, len(simplex)-1)
#                == True):
#                ec = ec + P_chains([x], [1])
#            else:
#                ec = ec - P_chains([x], [1])
#        return ec
#    def oriented_p_chains(self, k):
#        if ((k<0) or (k>self.dimension())):
#            return 0
#        else:
#            c_p = P_chains([], [1])
#            for x in self.p_simplex(k):
#                c_p = c_p + self.elementary_chain(tuple_sorted(x))
#            return c_p
    def basis_group_oriented_p_chains(self, p):
        """Gives a basis for the group of oriented p-chains.

        Args:
            p (int): Indicated the dimension of the p-simplex.

        Returns:
            __main__.P_chains: A new p-chains that contains the basis of the
            group.

            .. Note:: To every element is given the coefficiente ``1``
            by default, this is because the p-simplex are sorted with the
            lexicographical order, i.e, this orientation is taken positive.

        Raises:
            AttributeError: If p is lower than zero or bigger than the
            dimension of the simplicial complex dimension.

        Examples:
            To create a basis for the group of oriented p-chains, use
            ``SimplicialComplex.basis_group_oriented_p_chains(p)``.

            >>> import networkx as nx
            >>> from pyreps.simplicial import SimplicialComplex
            >>> from pyreps.graphs import matching_graph
            >>> G = matching_graph(3)
            >>> sc = SimplicialComplex(G)
            >>> sc.basis_group_oriented_p_chains(0).dic
            {((0, 1),): 1, ((0, 2),): 1, ((1, 2),): 1}

             .. Note:: We use the function ``matching_graph`` which
             will be explain after.

        """
        if ((p < 0) or (p > self.dimension())):
            return 0
        else:
            c_p = P_chains([], [])
            for x in self.p_simplex(p):
                c_p = c_p + P_chains([tuple(tuple_sorted(x))], [1])
            return c_p
#    def p_homology_group_dimention(self, k):
#        vk = self.simplex()[k]
#        vkf = self.n_faces(k-1)
#        M = zeros(len(vkf),len(vk.dic))
#        j=0
#        for u in list(vk.dic.keys()):
#            d={u: vk.dic[u]}
#            for a in list((boundary_op(d).dic).keys()):
#                i=0
#                for w in list(vkf):
#                    if (a == w):
#                        M[i,j]=(boundary_op(d).dic)[w]
#                    i=i+1
#            j=j+1
#        dimKe = len(M.rref()[1])
#        vk1 = self.simplex()[k+1]
#        vkf1 = self.n_faces(k)
#        N = zeros(len(vkf1),len(vk1.dic))
#        j=0
#        for u in list(vk1.dic.keys()):
#            d={u: vk1.dic[u]}
#            for a in list((boundary_op(d).dic).keys()):
#                i=0
#                for w in list(vkf1):
#                    if (a == w):
#                        N[i,j]=(boundary_op(d).dic)[w]
#                    i=i+1
#            j=j+1
#        dimIm = len((N.T).rref()[1])
#        dimH = dimKe - dimIm
#        return dimKe, dimIm, dimH

    def representate_in_simplex(self, vec, P):
        s = P_chains([], [])
        if (vec.dic != {}):
            v = list(vec.dic.keys())
            p = len(list(vec.dic.keys())[0]) - 1
            ve = self.basis_group_oriented_p_chains(p)
            for a in v:
                if isinstance(a, int):
                    return vec
                else:
                    w = tuple_permutation(a, P)
                    for b in list(ve.dic.keys()):
                        if eq_elements(b, w):
                            if orientation_function(b, w, p):
                                s = s + P_chains([b], [vec.dic[a]])
                            else:
                                s = s - P_chains([b], [vec.dic[a]])
            return s
        else:
            return s
    """
        Before see the documentation of the next methods, you must
        see the documentation of the functions.
    """
    def matrix_simmetric_representate(self, p):
        """Give the matrix associated to the boundary operator.

        Args:
            p (int): Determine with basis of the p-simplex will be used.
        Returns:
            A matrix if p is bigger than -1, and lower than the dimension
            of the simplicial complex, return False otherwise.

        Examples:
            To compute the matrix associated to the boundary operator, use
            ``SimplicialComplex.matrix_simmetric_representate(p)``.

            >>> import networkx as nx
            >>> from pyreps.simplicial import SimplicialComplex
            >>> from pyreps.graphs import matching_graph
            >>> G = matching_graph(3)
            >>> sc = SimplicialComplex(G)
            >>> sc.matrix_simmetric_representate(0)
            Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

            .. Note:: This matrix will be so useful to our purposes.

        """
        if (p > 0 and (p <= self.dimension())):
            v = self.basis_group_oriented_p_chains(p)
            p = p - 1
            ve = self.basis_group_oriented_p_chains(p)
            M = zeros(len(ve.dic), len(v.dic))
            j = 0
            for u1 in list(v.dic.keys()):
                d = P_chains([u1], [v.dic[u1]])
                for u2 in list(boundary_op_n(d).dic.keys()):
                    #                display(boundary_op(d, self.G).dic)
                    i = 0
                    for w in list(ve.dic.keys()):
                        if (w == u2):
                            # if (orientation_function(u2,tuple(w),p) == True):
                            M[i, j] = int((boundary_op_n(d).dic)[u2])
#                            else:
#                                M[i,j] = int((boundary_op_n(d).dic)[u2])*(-1)
                        i = i + 1
                j = j + 1
            return M
        else:
            if (p == 0):
                basisg = self.basis_group_oriented_p_chains(0).dic.keys()
                return eye(len(list(basisg)))
            else:
                return False

    def kernel_boundary_op(self, p):
        #        display(p,self.dimension())
        if ((p > 0) and (p <= self.dimension())):
            u = nullspace(self.matrix_simmetric_representate(p))
#            display(u)
            if (not u):
                v = self.basis_group_oriented_p_chains(p)
                w = []
                for i in range(len(u)):
                    s = P_chains([], [])
                    for j in range(len(u[i])):
                        if (u[i][j] != 0):
                            s = s + P_chains([list(v.dic.keys())[j]],
                                             [u[i][j]])
                    w.append(s)
                return w
            else:
                return False
        else:
            if (p == 0):
                return [self.basis_group_oriented_p_chains(p)]
            else:
                return False

    def image_boundary_op(self, p):
        #        display(p, self.dimension())
        if ((p > 0) and (p <= self.dimension())):
            u = columnspace(self.matrix_simmetric_representate(p))
            if (not u):
                v = self.basis_group_oriented_p_chains(p-1)
                w = []
                for i in range(len(u)):
                    s = P_chains([], [])
                    for j in range(len(u[i])):
                        if (u[i][j] != 0):
                            s = s + P_chains([list(v.dic.keys())[j]],
                                             [u[i][j]])
                    w.append(s)
                return w
            else:
                return False
        else:
            return False

    def character_kernel(self, p, P):
        """Gives character of a basis of the kernel under a Permutation.

        Args:
            p (int): Indicated the dimension of the p-simplex.
            P (sympy.combinatorics.permutations.Permutation): The permutation.

        Returns:
            int: The character of a basis of the kernel under a Permutation,
            in our case, we are interest in representates of the conjugacy
            class of the symmetric group.

        Examples:
            To get the character of a permutation acting on a basis
            of the kernel, use ``SimplicialComplex.character_kernel(p,
            Permutation)``.

            >>> import networkx as nx
            >>> from pyreps.simplicial import SimplicialComplex
            >>> from pyreps.graphs import matching_graph
            >>> G1 = matching_graph(4)
            >>> G = clique_graph(G1)
            >>> sc = SimplicialComplex(G)
            >>> print(sc.character_kernel(1,Permutation(0,1)))
            0
            >>> print(sc.character_kernel(1,Permutation(0,1,2,3)))
            0

        """
        A = self.matrix_simmetric_representate(p)
        if (p > 0 and (p <= self.dimension())):
            M = []
            null = nullspace(A)
            for i in range(len(null[0])):
                w = []
                for j in range(len(null)):
                    w.append(null[j][i])
                M.append(w)
            Q = Reduce(Matrix(M))
            M = Q[0]*Matrix(M)
        else:
            if (p == 0):
                M = A
                null = []
                for i in range(A.shape[0]):
                    aux = []
                    for j in range(A.shape[1]):
                        aux.append(M[i, j])
                    null.append(aux)
                Q = Reduce(Matrix(M))
                M = Q[0]*Matrix(M)
            else:
                return 0
        if (all(elem == null[0][0] for elem in null[0])):
            return 0
        else:
            w1 = []
            he = self.basis_group_oriented_p_chains(p)
            for a in range(len(null)):
                N = []
                v = P_chains([], [])
                c = 0
                for j in list(he.dic.keys()):
                    v = v + P_chains([j], [null[a][c]])
                    c = c+1
                v1 = permutation_in_simplex_test(v, P)
                u = []
                for i in list(he.dic.keys()):
                    for j in list(v1.dic.keys()):
                        if (i == j):
                            #  if (eq_elements(i, j) == True):
                            u.append(np.array([v1.dic[j]]))
                u = Q[0]*Matrix(u)
                N = np.append(M, u, axis=1)
                N = Matrix(N)
                w2 = []
                for i in tuple(linsolve(N)):
                    for j in i:
                        w2.append(j)
                w1.append(w2)
            N = Matrix(w1)
#            display(N.T)
            return np.trace(N.T)

    def character_image(self, p, P):
        """Gives character of a basis of the image under a Permutation.

        Args:
            p (int): Indicated the dimension of the p-simplex.
            P (sympy.combinatorics.permutations.Permutation): The permutation.

        Returns:
            int: The character of a basis of the image under a Permutation,
            in our case, we are interest in representates of the conjugacy
            class of the symmetric group.

        Examples:
            To get the character of a permutation acting on a basis
            of the image, use
            ``SimplicialComplex.character_image(p,Permutation)``.

            >>> n=5
            >>> G = matching_graph(n)
            >>> sc = SimplicialComplex(G)
            >>> print(sc.character_image(1,Permutation(0,1)))
            3
            >>> print(sc.character_image(1,Permutation()))
            9

        """
        if (p > 0 and (p <= self.dimension())):
            A = self.matrix_simmetric_representate(p)
            w1 = []
            M = []
            col = columnspace(A)
            for i in range(len(col[0])):
                w = []
                for j in range(len(col)):
                    w.append(col[j][i])
                M.append(w)
            Q = Reduce(Matrix(M))
            M = Q[0]*Matrix(M)
            he = self.basis_group_oriented_p_chains(p-1)
            for a in range(len(col)):
                N = []
                v = P_chains([], [])
                c = 0
                for j in list(he.dic.keys()):
                    v = v + P_chains([j], [col[a][c]])
                    c = c+1
                v1 = permutation_in_simplex_test(v, P)
                u = []
                for i in list(he.dic.keys()):
                    for j in list(v1.dic.keys()):
                        if (i == j):
                            # if (eq_elements(i, j) == True):
                            u.append(np.array([v1.dic[j]]))
                u = Q[0]*Matrix(u)
                N = np.append(M, u, axis=1)
                N = Matrix(N)
                w2 = []
                for i in tuple(linsolve(N)):
                    for j in i:
                        w2.append(j)
                w1.append(w2)
            N = Matrix(w1)
#            display(N.T)
            return np.trace(N.T)
        else:
            return 0

    def character_p_homology(self, p, P):
        """Gives character of the pth homology.

        Args:
            p (int): Indicated the dimension of the p-simplex.
            P (sympy.combinatorics.permutations.Permutation): The permutation.

        Returns:
            int: The character of the pth homology.

        Examples:
            To get the character of the pth homology use
            ``SimplicialComplex.character_p_homology(p, Permutation)``.

            >>> n=6
            >>> G1 = matching_graph(n)
            >>> G = clique_graph(G1)
            >>> sc = SimplicialComplex(G)
            >>> print(sc.character_p_homology(1,Permutation(0,1)))
            0
            >>> print(sc.character_p_homology(1,Permutation()))
            16

            ..Note:: The funcion only is the subtract of the
                character of the kernel (dimension p) and the
                character of the image (dimension p-1).

        """
        return self.character_kernel(p, P) - self.character_image(p + 1, P)

    def specific_function(self, n):
        """Returns a dictionary showing the descomposition into irreducibles.

        Args:
            n (int): Indicated what symmetric group act on the p-simplex.

        Returns:
            dict: A dictionary that contains show the descomposition
            into irreducibles, i.e the reduce homologies of the
            simplicial complex determinated by a graph.

        Examples:
            To get the reduce homologies of the simplicial complex by a
            graph, use ``SimplicialComplex.specific_function(n)``.

            >>> n=4
            >>> G1 = matching_graph(n)
            >>> G = clique_graph(G1)
            >>> sc1 = SimplicialComplex(G1)
            >>> print(sc1.specific_function(n))
            {0: {(4,): 1, (1, 1, 1, 1): 0, (2, 1, 1): 0, (2, 2): 1, (3, 1): 0}, 1: {(4,): 0, (1, 1, 1, 1): 0, (2, 1, 1): 0, (2, 2): 0, (3, 1): 0}}
            >>> sc2 = SimplicialComplex(G)
            >>> print(sc2.specific_function(n))
            {0: {(4,): 1, (1, 1, 1, 1): 0, (2, 1, 1): 0, (2, 2): 1, (3, 1): 0}}

        """
        w = partitions_list(n)
        M = form_matrix_yt(w)
        card = math.factorial(n)
        vec_dic = {}
        for k in range(self.dimension()+1):
            D = {}
            u = []
            v = []
            for h in w:
                u.append(self.character_p_homology(k, make_permutation(h)))
                v.append(size_conjugacy_class(h, n))
            for i in range(M.shape[0]):
                Ip = 0
                for j in range(M.shape[1]):
                    Ip = Ip + M[i, j]*u[j]*v[j]
                Ip = Ip/card
                D[tuple(w[i])] = Ip
            vec_dic[k] = D
        return vec_dic

    def character_matrix_permutation(self, p, P):
        #        display(1)
        v = self.basis_group_oriented_p_chains(p)
        #        display(2)
        v1 = permutation_in_simplex_test(v, P)
        #        display(3)
        M = zeros(len(v.dic.keys()), len(v.dic.keys()))
        i = 0
        for u1 in v.dic.keys():
            #            display(4)
            j = 0
            for u2 in v1.dic.keys():
                if (u1 == u2):
                    M[i, j] = (v1.dic)[u2]
                j = j + 1
            i = i + 1
        return np.trace(M)

    def specific_function_1(self, n):
        w = partitions_list(n)
        M = form_matrix_yt(w)
        card = math.factorial(n)
        vec_dic = {}
        for k in range(3):
            D = {}
            u = []
            v = []
            for h in w:
                mpk = self.character_matrix_permutation(k, make_permutation(h))
                u.append(mpk)
                v.append(size_conjugacy_class(h, n))
            for i in range(M.shape[0]):
                Ip = 0
                for j in range(M.shape[1]):
                    Ip = Ip + M[i, j]*u[j]*v[j]
                Ip = Ip/card
                D[tuple(w[i])] = Ip
            vec_dic[k] = D
        return vec_dic

    def specific_function_2(self, p, n):
        w = partitions_list(n)
        print(1)
        M = form_matrix_yt(w)
        card = math.factorial(n)
        print(2)
        v = self.basis_group_oriented_p_chains(p)
        print(3)
        leng = len(v.dic.keys())
        print(4)
        M1 = zeros(leng, leng)
        print(5)
        vec_dic = {}
        D = {}
        au = []
        av = []
        for h in w:
            print(6)
            M2 = M1.copy()
            print(7)
            v1 = permutation_in_simplex_test(v, make_permutation(h))
            print(8)
            i = 0
            for u1 in v.dic.keys():
                #            display(4)
                j = 0
                for u2 in v1.dic.keys():
                    if (u1 == u2):
                        M2[i, j] = (v1.dic)[u2]
                    j = j + 1
                i = i + 1
            traza = np.trace(M2)
            au.append(traza)
            av.append(size_conjugacy_class(h, n))
            print(traza, h)
        for i in range(M.shape[0]):
            Ip = 0
            for j in range(M.shape[1]):
                Ip = Ip + M[i, j]*au[j]*av[j]
            Ip = Ip/card
            D[tuple(w[i])] = Ip
        vec_dic[p] = D
        return vec_dic

    def specific_function_5(self, p, n, valores):
        w = [[7], [1, 1, 1, 1, 1, 1, 1], [2, 1, 1, 1, 1, 1], [2, 2, 1, 1, 1],
             [2, 2, 2, 1], [3, 1, 1, 1, 1], [3, 2, 1, 1], [3, 2, 2], [3, 3, 1],
             [4, 1, 1, 1], [4, 2, 1], [4, 3], [5, 1, 1], [5, 2], [6, 1]]
        #        w = partitions_list(n)
        M = form_matrix_yt(w)
        card = math.factorial(n)
        vec_dic = {}
        D = {}
        au = valores
        av = []
        for h in w:
            av.append(size_conjugacy_class(h, n))
        for i in range(M.shape[0]):
            Ip = 0
            for j in range(M.shape[1]):
                Ip = Ip + M[i, j]*au[j]*av[j]
            Ip = Ip/card
            D[tuple(w[i])] = Ip
        vec_dic[p] = D
        return vec_dic
