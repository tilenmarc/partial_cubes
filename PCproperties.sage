from sage.graphs.convexity_properties import ConvexityProperties
from sage.all_cmdline import *


class PCproperties:
    def __init__(self, G, precompute=True):
        self.G = G.copy()

        check, self.labels = G.is_partial_cube(certificate=True)
        assert check, 'G is not a partial cube'
        self.G.relabel(self.labels)
        self.dimension = len(self.labels.values()[0])
        self.COM = None
        self.OM = None
        self.LOP = None
        self.median = None
        self.almost_median = None
        self.pasch = None
        self.tree_zone = None
        self.cellular = None
        self.hypercellular = None
        self.well_embedded = None
        self.netlike = None
        self.polat = None
        self.peano = None
        self.affine = None
        if precompute:
            self.precompute_values()

    def precompute_values(self):
        self.distances = self.G.distance_all_pairs()
        self.CP = ConvexityProperties(self.G)

    def is_antipodal_subgraph(self, C):
        if len(C) == 1:
            return True
        for x in C:
            max_dist = 0
            for y in C:
                if self.distances[x][y] > max_dist:
                    max_dist = self.distances[x][y]
                    antipode = y
            if Set(C) != Set(self.CP.hull([x, antipode])):
                return False
        else:
            return True

    def is_gated(self, C):
        for a in self.G.vertices():
            if a in C:
                continue
            min_dist = self.G.order()
            for x in C:
                if self.distances[x][a] < min_dist:
                    gate = x
                    min_dist = self.distances[x][a]
            for x in C:
                if self.distances[x][a] != self.distances[x][gate] + self.distances[gate][a]:
                    return False
        return True

    def is_COM(self):
        for x, y in Combinations(self.G.vertices(), 2):
            C = self.CP.hull([x, y])
            if self.is_antipodal_subgraph(C) and not self.is_gated(C):
                self.COM = False
                return False
        self.COM = True
        return True

    def is_LOP(self):
        for x, y in Combinations(self.G.vertices(), 2):
            C = self.CP.hull([x, y])
            if self.is_antipodal_subgraph(C) and len(C) != 2 ** (self.distances[x][y]):
                self.LOP = False
                return False
        self.LOP = True
        return True

    def is_OM(self):
        if self.COM == None:
            self.is_COM()
        return self.COM and self.is_antipodal_subgraph(self.G.vertices())

    def has_minor(self, s, i):
        # s is a list of canonical labels of minors, i is a list of their isometric dimensions
        for j in range(len(s)):
            H = s[j]
            for S in Combinations(range(self.dimension), self.dimension - i[j]):
                for contract in Combinations(S):
                    restrict = [x for x in S if x not in contract]
                    for positive in Combinations(restrict):
                        negative = [x for x in restrict if x not in positive]
                        H = self.G.copy()
                        H = contract_and_restrict(H, contract, negative, positive)
                        if H.canonical_label().graph6_string() in s:
                            return True
        return False

    def is_median(self):
        s = [graphs.CycleGraph(6).canonical_label().graph6_string()]
        Q = graphs.CubeGraph(3)
        Q.delete_vertex('000')
        s.append(Q.canonical_label().graph6_string())
        self.median = not self.has_minor(s, [3, 3])
        return self.median

    def is_median_check(self):
        for u, v, w in Combinations(self.G.vertices(), 3):
            m = ''
            for i in range(len(u)):
                if u[i] in [v[i], w[i]]:
                    m += u[i]
                else:
                    m += v[i]
            if m not in self.G.vertices():
                return False
        return True

    def is_loopsided(self):
        s = [graphs.CycleGraph(6).canonical_label().graph6_string()]
        Q = graphs.CubeGraph(4)
        Q.delete_vertex('0000')
        Q.delete_vertex('1111')
        s.append(Q.canonical_label().graph6_string())
        self.cellular = not self.has_minor(s, [3, 4])
        return self.cellular

    def is_hypercellular(self):
        s = []
        Q = graphs.CubeGraph(3)
        Q.delete_vertex('000')
        s.append(Q.canonical_label().graph6_string())
        self.hypercellular = not self.has_minor(s, [3])
        return self.hypercellular

    def is_cellular(self):
        s = [graphs.CubeGraph(3).canonical_label().graph6_string()]
        Q = graphs.CubeGraph(3)
        Q.delete_vertex('000')
        s.append(Q.canonical_label().graph6_string())
        self.cellular = not self.has_minor(s, [3, 3])
        return self.cellular

    def is_almost_median(self):
        s = [graphs.CycleGraph(6).canonical_label().graph6_string()]
        self.almost_median = not self.has_minor(s, [3])
        return self.almost_median

    def is_pasch(self):
        s = []
        Q = graphs.CubeGraph(4)
        Q.delete_vertex('0000')
        Q.delete_vertex('1110')
        s.append(Q.canonical_label().graph6_string())
        Q.delete_vertex('1111')
        s.append(Q.canonical_label().graph6_string())
        Q.delete_vertex('1101')
        s.append(Q.canonical_label().graph6_string())
        Q.delete_vertex('1011')
        s.append(Q.canonical_label().graph6_string())
        Q.delete_vertex('0111')
        s.append(Q.canonical_label().graph6_string())
        Q = graphs.CubeGraph(4)
        Q.delete_vertex('0000')
        s.append(Q.canonical_label().graph6_string())
        Q.delete_vertex('1111')
        s.append(Q.canonical_label().graph6_string())
        self.pasch = not self.has_minor(s, 7 * [4])
        return self.pasch

    def is_well_embedded(self):
        s = []
        Q = graphs.CubeGraph(4)
        Q.delete_vertex('0000')
        Q.delete_vertex('1110')
        s.append(Q.canonical_label().graph6_string())
        Q.delete_vertex('1111')
        s.append(Q.canonical_label().graph6_string())
        Q.delete_vertex('1101')
        s.append(Q.canonical_label().graph6_string())
        Q.delete_vertex('1011')
        s.append(Q.canonical_label().graph6_string())
        Q.delete_vertex('0111')
        s.append(Q.canonical_label().graph6_string())
        self.well_embedded = not self.has_minor(s, 5 * [4])
        return self.well_embedded

    def is_tree_zone(self):
        Q = graphs.CubeGraph(3)
        s = [Q.canonical_label().graph6_string()]
        Q = graphs.CubeGraph(4)
        Q.delete_vertex('0000')
        Q.delete_vertex('1110')
        Q.delete_vertex('1111')
        Q.delete_vertex('1101')
        Q.delete_vertex('1011')
        Q.delete_vertex('0111')
        s.append(Q.canonical_label().graph6_string())
        self.tree_zone = not self.has_minor(s, [3, 4])
        return self.tree_zone

    def is_polat(self):
        Q = graphs.CycleGraph(6).cartesian_product(graphs.CubeGraph(1))
        s = [Q.canonical_label().graph6_string()]
        Q = graphs.CubeGraph(3)
        Q.delete_vertex('000')
        s.append(Q.canonical_label().graph6_string())
        self.polat = not self.has_minor(s,  [4, 3])
        return self.polat

    def is_peano(self):
        for u in self.G.vertices():
            for v in self.G.vertices():
                for w in self.G.vertices():
                    for x in self.CP.hull([u, w]):
                        for y in self.CP.hull([v, x]):
                            for z in self.CP.hull([v, w]):
                                if y in self.CP.hull([u, z]):
                                    break
                            else:
                                self.peano = False
                                return False
        self.peano = True
        return True

    def is_netlike(self):
        if not self.peano:
            return False
        theta_classes = {i: [] for i in xrange(self.dimension)}
        for edge in self.G.edges():
            for i in xrange(self.dimension):
                if edge[0][i] != edge[1][i]:
                    theta_classes[i].append(edge)
                    break
        for i in xrange(self.dimension):
            U0, U1 = Set(), Set()
            u0, u1 = Set(), Set()
            for x, y in Combinations(theta_classes[i], 2):
                u0 = u0.union(Set([x[0], y[0]]))
                u1 = u1.union(Set([x[1], y[1]]))
                U0 = U0.union(Set(self.CP.hull([x[0], y[0]])))
                U1 = U1.union(Set(self.CP.hull([x[1], y[1]])))
            G0 = self.G.subgraph(U0)
            G1 = self.G.subgraph(U1)
            for v in G0.vertices():
                if G0.degree(v) >= 3 and v not in u0:
                    self.netlike = False
                    return False
            for v in G1.vertices():
                if G1.degree(v) >= 3 and v not in u1:
                    self.netlike = False
                    return False
        self.netlike = True
        return True

    def isometric_cycles_gated(self):
        for i in [6]:
            for edges in Combinations(self.G.edges(), i):
                H = Graph(edges)
                if H.is_connected():
                    if H.order() == i:
                        if H.canonical_label() ==\
                                graphs.CycleGraph(i).canonical_label():
                            distH = H.distance_all_pairs()
                            for x, y in Combinations(H.vertices(), 2):
                                if self.distances[x][y] != distH[x][y]:
                                    break
                            else:
                                hull = self.CP.hull(H.vertices())
                                if not self.is_gated(hull):
                                    return False
        else:
            return True


    def is_netlike_check(self):
        if self.polat:
            return True
        if self.tree_zone and self.peano:
            return True
        return False

    def is_affine(self):
        Q = graphs.CubeGraph(self.dimension + 1)
        Qvertices = [x + '0' for x in self.G.vertices()]
        Qvertices += [''.join(['0' if x == '1' else '1' for x in y]) + '1' for
            y in self.G.vertices()]
        Q = Q.subgraph(Qvertices)
        return Q.is_partial_cube(), Q

    def find_counterexample(self):
        for e in range(self.dimension):
            contracted = contract_and_restrict(self.G, [e], [], [])
            min_deg = min(contracted.degree_sequence())
            print min_deg, contracted.degree_sequence().count(min_deg)
        for e in Combinations(range(self.dimension), 2):
            contracted = contract_and_restrict(self.G, e, [], [])
            min_deg = min(contracted.degree_sequence())
            print min_deg, contracted.degree_sequence().count(min_deg)

    # works only for loopsided
    def euclidean(self, max_dim):
        cocircuits = []
        for x, y in Combinations(self.G.vertices(), 2):
            if self.distances[x][y] != max_dim:
                continue
            H_vertices = self.CP.hull([x, y])
            if len(H_vertices) != 2 ** max_dim:
                continue
            cocircuits.append(tuple(H_vertices))
        cograph = Graph()
        for x, y in Combinations(cocircuits, 2):
            intersect = [i for i in x if i in y]
            if len(intersect) == 2 ** (max_dim - 1):
                cograph.add_edge(x, y)
        cograph.plot(vertex_labels=False).save('cograph.png')
        class_to_cube = {}
        for x, y, w in Combinations(range(self.dimension), 3):
            intersect = False
            for z in cograph.vertices():
                if len(set([v[x] for v in z])) == 2:
                    if len(set([v[y] for v in z])) == 2:
                        if len(set([v[w] for v in z])) == 2:
                            intersect = True
                            class_to_cube[tuple(sorted([x, y, w]))] = z
                            break
            if not intersect:
                print x, y, w, 'not_intersect'

        def covector_dist(x, y):
            return min([self.distances[i][j] for i in x for j in y])

        for x in range(self.dimension):
            dicograph = DiGraph()
            for e in cograph.edges():
                # print e[0]
                dif1 = [i for i in range(self.dimension) if len(set([v[i] for v in e[0]])) == 2]
                dif2 = [i for i in range(self.dimension) if len(set([v[i] for v in e[1]])) == 2]
                dif = [y for y in dif1 if y in dif2]
                # print dif
                if x in dif:
                    continue
                if covector_dist(e[0], class_to_cube[tuple(sorted([dif[0], dif[1], x]))]) <\
                        covector_dist(e[1], class_to_cube[tuple(sorted([dif[0], dif[1], x]))]):
                    dicograph.add_edge([e[0], e[1]])
                else:
                    dicograph.add_edge([e[1], e[0]])
            dicograph.plot(vertex_labels=False).save('dicograph.png')
            if not dicograph.is_directed_acyclic():
                dicograph.plot(vertex_labels=False).save('dicograph.png')
                print 'not acyclic', x

    def is_appiculate(self):
        for v in self.G.vertices():
            for u, w in Combinations(self.G.vertices(), 2):
                inter = set(self.CP.hull([v, u])).intersection(set(self.CP.hull([v, w])))
                upper_v = v
                dist = 0
                for x in inter:
                    current_dist = self.distances[x][v]
                    if current_dist > dist:
                        dist = current_dist
                        upper_v = x
                if inter != set(self.CP.hull([v, upper_v])):
                    return False
        return True

    def is_appiculate_local(self, v):
        for u, w in Combinations(self.G.vertices(), 2):
            inter = set(self.CP.hull([v, u])).intersection(set(self.CP.hull([v, w])))
            upper_v = v
            dist = 0
            for x in inter:
                current_dist = self.distances[x][v]
                if current_dist > dist:
                    dist = current_dist
                    upper_v = x
            if inter != set(self.CP.hull([v, upper_v])):
                return False
        return True

    # check if all its antipodal subgraphs are prisms, works to
    # isometric dimension 6 of subgraphs
    def is_COP(self):
        # prisms are the Cartesian product of edges and cycles up to
        # isometric dimension 6
        one_vertex = Graph()
        one_vertex.add_vertex(0)
        hypercubes6 = [one_vertex] + [graphs.CubeGraph(i) for i in range(1, 7)]
        prisms = [x.canonical_label().graph6_string()
                                for x in hypercubes6]
        prisms += [graphs.CycleGraph(6).cartesian_product(x).canonical_label().graph6_string()\
                                 for x in hypercubes6[:4]]
        prisms += [graphs.CycleGraph(8).cartesian_product(x).canonical_label().graph6_string()\
                                 for x in hypercubes6[:3]]
        prisms += [graphs.CycleGraph(10).cartesian_product(x).canonical_label().graph6_string()\
                                 for x in hypercubes6[:2]]
        prisms += [graphs.CycleGraph(12).canonical_label().graph6_string()]
        prisms += [graphs.CycleGraph(6).cartesian_product(graphs.CycleGraph(6)).canonical_label().graph6_string()]
        prisms = set(prisms)

        for x, y in Combinations(self.G.vertices(), 2):
            C = self.CP.hull([x, y])
            if self.is_antipodal_subgraph(C):
                name = self.G.subgraph(C).canonical_label().graph6_string()
                if name not in PRISMS:
                    return False
        return True

    def get_rank(self):
        self.rank = 1
        hypercubes = [graphs.CubeGraph(i).canonical_label().graph6_string()
                        for i in range(1, self.dimension + 1)]
        for i in range(2, self.dimension + 1):
            if self.has_minor([hypercubes[i - 1]], [i]):
                self.rank = i
            else:
                break
        return self.rank

    # finds all the mutations of self
    def find_mutations(self, rank):
        ret = []
        cube = graphs.CubeGraph(self.dimension)

        for v in self.G.vertices():
            if self.G.degree(v) > rank:
                continue
            indices = []
            for u in self.G[v]:
                for i in range(len(u)):
                    if u[i] != v[i]:
                        indices.append(i)
                        break
            v1 = flip(v, indices)
            H_vertices = self.G.vertices()
            H_vertices.remove(v)
            H_vertices.remove(flip(v, range(self.dimension)))
            H_vertices.append(v1)
            H_vertices.append(flip(v1, range(self.dimension)))
            H = cube.subgraph(H_vertices)
            mut = H.canonical_label().graph6_string()
            if mut not in ret:
                ret.append(mut)
        return ret

def flip(v, indices):
    return ''.join([v[i] if i not in indices else ('0' if v[i] == '1' else '1')
                    for i in range(len(v))])

def antipode(x):
    return ''.join(['0' if s == '1' else '1' for s in x])

# helping function: assumes vertices are string of zeros and ones
def contract_and_restrict(H, contracted_dimensions, restricted_dimensions0, restricted_dimensions1):
    I = Graph()
    for edge in H.edges():
        check = True
        for i in restricted_dimensions0:
            if edge[0][i] != "0" or  edge[1][i] != '0':
                check = False
                break
        if not check:
            continue
        for i in restricted_dimensions1:
            if edge[0][i] != "1" or  edge[1][i] != '1':
                check = False
                break
        if not check:
            continue
        for i in contracted_dimensions:
            if edge[0][i] != edge[1][i]:
                check = False
                break
        if not check:
            continue
        u = ''.join([edge[0][i] for i in range(len(edge[0])) if i not in contracted_dimensions + restricted_dimensions0 + restricted_dimensions1])
        v = ''.join([edge[1][i] for i in range(len(edge[0])) if i not in contracted_dimensions + restricted_dimensions0 + restricted_dimensions1])
        I.add_edge(u, v)
    return I

