"""
Created using PyCharm
@author: Madalina Giurgiu
@date: 5:30 PM 9/19/22

Definition and methods of the underlying graph.
"""
import math
from collections import defaultdict
from copy import deepcopy
import numpy as np

import decoil.search.cycles as cycles
from decoil.utils import FRAGMENT as fg
from decoil.utils import GRAPH_EDGE_PROP as gep
from decoil.utils import GRAPH_NODE_PROP as gnp
from decoil.utils import GRAPH_PROP as gp
from decoil.utils import QUAL
from intervaltree import IntervalTree
from overrides import overrides

NODE_REF = 0
STATE = 1
VISITED = 2
WEIGHT = 3


class FragmentsIntervals(object):
        """
        Tracks the collection of fragments
        """

        def __init__(self):
                self.tree_fragments = defaultdict(IntervalTree)

        def add_fragment(self, chr, start, stop, data):
                """
                Add fragment to interval tree + data information
                """
                # if chr not in self.tree_fragments:
                #       self.tree_fragments[chr] = IntervalTree()
                self.tree_fragments[chr].addi(int(start), int(stop), data)

        # def remove_fragment(self, chr, start, stop):
        #       self.tree_fragments[chr].remove_overlap(int(start) + 2 * QUAL.DISTANCE, int(stop) - 2 * QUAL.DISTANCE)

        def remove_fragment(self, chr, pos):
                self.tree_fragments[chr].remove_overlap(pos)

        def overlap_interval(self, chr, start, stop):
                if chr not in self.tree_fragments:
                        return None
                return self.tree_fragments[chr].overlap(int(start), int(stop))

        def overlap_breakpoint(self, chr, pos):
                if chr not in self.tree_fragments:
                        return None
                return self.tree_fragments[chr].at(int(pos))

        def find_tail_node(self, chr, pos):
                """
                Return the node id for fragment tail
                """
                results = list(self.overlap_breakpoint(chr, int(pos) + 2 * QUAL.DISTANCE))
                if len(results) == 1:
                        return results[0].data[fg.NODE_TAIL]
                elif len(results) > 1:
                        raise Exception("Overlapping fragments", results)
                return None

        def find_head_node(self, chr, pos):
                """
                Return the node id for fragment head
                """
                results = list(self.overlap_breakpoint(chr, int(pos) - 2 * QUAL.DISTANCE))
                if len(results) == 1:
                        return results[0].data[fg.NODE_HEAD]
                elif len(results) > 1:
                        raise Exception("Overlapping fragments", results)
                return None

        def __str__(self):
                return self.tree_fragments.__str__()


class Fragment(object):
        """
        Genomic fragment representation
        """

        def __init__(self, _id, _id_tail, _id_head, _id_edge, _props):
                self._id = _id
                self._id_tail = _id_tail
                self._id_head = _id_head
                self._id_edge = _id_edge
                self._props = _props
                self._len = int(abs(self._props[fg.START] - self._props[fg.END]))

        def set_norm_cov(self, meancov):
                """
                Compute normalised coverage = (cov / (wgs_meancov + 1)) * (fragment_len / 10^6)
                """
                self._props[fg.NORM_COV] = float(self._props[fg.COVERAGE] / (meancov + 1) * (fg.MB / self._len))

        def get_norm_cov(self):
                return self._props[fg.NORM_COV]

        def is_overlap_with(self, fobj):
                """
                Check if current fragment overlaps with another

                Arguments:
                        fobj (decoil.graph.Fragment)
                """
                if self.chr == fobj.chr:
                        if self.props[fg.STRAND] == "+" and fobj.props[fg.STRAND] == "+" and abs(self.end - fobj.start) < 10:
                                return True
                        if self.props[fg.STRAND] == "-" and fobj.props[fg.STRAND] == "-" and abs(self.start - fobj.end) < 10:
                                return True

                return False


        @property
        def id(self):
                return self._id

        @property
        def head(self):
                return self._id_head

        @head.setter
        def head(self, id):
                self._id_head = int(id)

        @property
        def tail(self):
                return self._id_tail

        @tail.setter
        def tail(self, id):
                self._id_tail = int(id)

        @property
        def edge(self):
                return self._id_edge

        @edge.setter
        def edge(self, id):
                self._id_edge = int(id)

        @property
        def coverage(self):
                return self._props[fg.COVERAGE]

        @coverage.setter
        def coverage(self, cov):
                self._props[fg.COVERAGE] = float(cov)

        @property
        def len(self):
                return self._len

        @len.setter
        def len(self, size):
                self._len = float(size)


        @property
        def coords(self):
                return """{}:{}-{}""".format(self._props[fg.CHR],
                                             str(self._props[fg.START]),
                                             str(self._props[fg.END]),
                                             )

        @property
        def chr(self):
                return self._props[fg.CHR]

        @property
        def start(self):
                """
                Returns:
                        (int) start of the genomic fragment
                """
                return int(self._props[fg.START])

        @start.setter
        def start(self, _start):
                self._props[fg.START] = int(_start)

        @property
        def end(self):
                """
                Returns:
                        (int) end of the genomic fragment
                """
                return int(self._props[fg.END])

        @end.setter
        def end(self, _end):
                self._props[fg.END] = int(_end)

        @property
        def props(self):
                return self._props


class Node(object):

        def __init__(self, _id, _part, _parent_fragment, _data):
                self._id = _id
                self._parent_fragment = _parent_fragment
                self._part = _part  # 0 /1
                self._props = _data
                self._idpair = None
                self._start_time = None
                self._end_time = None

        @property
        def id(self):
                return self._id

        @property
        def parent_fragment(self):
                return self._parent_fragment

        @property
        def part(self):
                return self._part

        @property
        def props(self):
                return self._props

        @property
        def idpair(self):
                return self._idpair

        @idpair.setter
        def idpair(self, id):
                """
                Set the node pair corresponding to same fragment;
                A fragment is composed of 2 paired nodes
                """
                self._idpair = id

        @property
        def start_time(self):
                """
                Return iteration time
                """
                return self._start_time

        @start_time.setter
        def start_time(self, t):
                """
                Set iteration time
                """
                self._start_time = t

        @property
        def end_time(self):
                """
                Return iteration time
                """
                return self._end_time

        @end_time.setter
        def end_time(self, t):
                """
                Set iteration timne
                """
                self._end_time = t

        @overrides
        def __str__(self):
                return """id:{} part:{} parent_fragment:{} idpair:{} data:{}""".format(self._id,
                                                                                       self._part,
                                                                                       self._parent_fragment,
                                                                                       self._idpair,
                                                                                       str(self._props))


class Edge(object):

        def __init__(self, _id, _u, _v, _svtype, _data):
                self._id = _id
                self._u = _u
                self._v = _v
                self._svtype = _svtype
                self._props = _data

        @overrides
        def __str__(self):
                return """id:{} u:{} v:{} svtype:{} data:{}""".format(self._id,
                                                                      self._u,
                                                                      self._v,
                                                                      self._svtype,
                                                                      str(self._props))

        @property
        def u(self):
                return self._u

        @property
        def v(self):
                return self._v

        @property
        def svtype(self):
                return self._svtype

        @property
        def props(self):
                return self._props

        @property
        def weight(self):
                return self._props[gep.WEIGHT]


class MultiGraph(object):
        """
        Multigraph with special properties.
        Fragments are meta-nodes containing two simple nodes (head and tail)
        """

        # this is an unique number for every node and edge
        id_generator = 0

        def __init__(self):

                self.graph = defaultdict(dict)  # edges list for every u,v are ordered by weight (descendent)
                self.fragments = defaultdict(Fragment)
                self.fragments_intervals = FragmentsIntervals()
                self.nodes = defaultdict(Node)
                self.edges = defaultdict(Edge)
                self.paths = []
                self.start_time = 0  # used to generate DFS tree
                self.end_time = 0
                self.hash_circles = defaultdict(int) # hash individual cycles
                self.hash_canonical = defaultdict(int)
                # self.metahash_circles = defaultdict(int) # hash equivalent cycles

        def adjacency_matrix(self):
                """
                Return adjacency matrix for the multi-graph.
                """
                size = len(self.nodes.keys())
                m = np.zeros((size,size),dtype=int)
                cols = list(self.nodes.keys())
                rows = list(self.nodes.keys())
                
                print(cols)
                for u in self.graph:
                        print(u, self.graph[u])
                        for v in self.graph[u]:
                                i = rows.index(u)
                                j = cols.index(v)
                                m[i][j] += 1

                return m

        def add_fragment(self, id, id_tail, id_head, id_edge, data):
                """
                Add fragment to the graph.

                Arguments:
                        id (int): UID
                        id_tail (int): node id of the fragment's tail
                        id_head (int): node id of the fragment's head
                        id_edge (int): edge id of the connecting edge
                        data (dict): fragment information {"coverage":x,
                                                                                                "strand": 1,
                                                                                                "chr": x,
                                                                                                "start": 111,
                                                                                                "end": 123}
                """
                if id not in self.fragments:
                        self.fragments[id] = Fragment(id, id_tail, id_head, id_edge, data)

        def add_node(self, id, part, parent_fragment, data):
                """
                Add node to graph
                
                Arguments:
                        id (int): unique identifier of the node
                        part (int): 0-tail, 1-head of the fragment
                        parent_fragment (str): name of the parent fragment
                        data (dict): properties of the node
                """
                if id not in self.nodes:
                        self.nodes[id] = Node(id, part, parent_fragment, data)
                        self.graph[id] = {}

        def add_edge(self, id, u, v, svtype, data):

                if id not in self.edges:
                        self.edges[id] = Edge(id, u, v, svtype, data)
                        if v not in self.graph[u]:
                                self.graph[u][v] = []
                        if u not in self.graph[v]:
                                self.graph[v][u] = []

                        self.graph[u][v].append(id)
                        self.graph[v][u].append(id)

        def remove_node(self, u):
                """
                Remove node and all adjacent edges
                """

                # delete adjacent edges
                toremove = []
                for v in self.neighbors(u):
                        for e in self.graph[u][v]:
                                if e in self.edges:
                                        self.edges.pop(e)
                        toremove.append(v)

                # remove u entries from the adjancency matrix
                for key in toremove:
                        self.graph[key].pop(u)
                self.graph.pop(u)

                # delete node itself
                self.nodes.pop(u)

        def remove_edge(self, edgeid, u, v):
                """
                Remove edge;
                (1) simple removal if edge represents a structural variant
                (2) remove edge and its parent fragment if edge connects fragment's head to tail

                Arguments:
                        edgeid (int): edge id connecting u und v
                        u (int): left connecting node id
                        v (int): right connecting node id
                """
                self.graph[u][v].remove(edgeid)
                self.graph[v][u].remove(edgeid)

                # remove fragment if edge was connecting fragment;s head to tail
                if self.edges.get(edgeid).svtype == gp.FRAGMENT:
                        fid = self.nodes.get(u).parent_fragment
                        self.remove_fragment(fid)

                # remove edge from edges dict
                self.edges.pop(edgeid)

        def remove_fragment(self, fid):
                """
                Remove genomic fragments from the data structure
                """
                if fid in self.fragments:
                        frag_obj = self.fragments[fid]
                        u = frag_obj.tail
                        v = frag_obj.head

                        # remove fragment from interval tree
                        fragment_middle = (self.nodes[u].props[gnp.POS] + self.nodes[v].props[gnp.POS]) / 2
                        self.fragments_intervals.remove_fragment(self.nodes[u].props[gnp.CHR], fragment_middle)

                        # remove nodes (this automatically removes also the edge)
                        self.remove_node(u)
                        self.remove_node(v)

                        # remove the fragment object
                        self.fragments.pop(fid)

        def change_fragment_boundry(self, fid, part, offset):
                """
                Change fragment boundries, i.e. elongate
                """
                if part == gp.TAIL:
                        # elongate to left
                        self.fragments[fid].start = self.fragments[fid].start - offset
                        self.fragments[fid].len = self.fragments[fid].len - offset
                        # change props for the node
                        node = self.fragments[fid].tail
                        self.nodes[node].props[gnp.POS] = self.nodes[node].props[gnp.POS] - offset

                elif part == gp.HEAD:
                        # elongate to right
                        self.fragments[fid].end = self.fragments[fid].end + offset
                        self.fragments[fid].len = self.fragments[fid].len + offset
                        # change props for the node
                        node = self.fragments[fid].head
                        self.nodes[node].props[gnp.POS] = self.nodes[node].props[gnp.POS] + offset


        def rewire_fragment_neighbors(self, fid, dist_neighbors = 50):
                """
                Connect neighbors of a fragment f. Assume
                       fa -- f -- fb
                           should harbor fa --- fb

                Arguments:
                                fid (int): Fragment id
                """
                to_rewire = []
                if fid in self.fragments:
                        frag_obj = self.fragments[fid]
                        u = frag_obj.tail
                        v = frag_obj.head

                        neighbors_u = self.neighbors(u)
                        neighbors_v = self.neighbors(v)

                        # connect these neighbors (if edge type does not exist)
                        for ui in neighbors_u:
                                for vj in neighbors_v:
                                        # ignore edges connected to fragment fid
                                        if ui == v or vj == u:
                                                continue
                                        # add only if edge does not exist
                                        if vj not in self.graph[ui]:
                                                props = {gep.COLOR: "gray",
                                                          gep.SVTYPE: gp.REWIRE,
                                                          gep.WEIGHT: math.inf}
                                                to_rewire.append((ui,vj,props))

                        # try to extend neighbor fragment (left)
                        for ui in neighbors_u:
                                fid_ui = self.nodes[ui].parent_fragment
                                fid_ui_obj = self.fragments[fid_ui]
                                # if ui is not the pair of fid
                                # if ui is the head of the fragment fid_ui
                                # if ui and u are in close proximity
                                if ui != v and fid_ui_obj.head == ui and abs(fid_ui_obj.end - frag_obj.start) < dist_neighbors:
                                        self.change_fragment_boundry(fid_ui, gp.HEAD, frag_obj.len / 2)

                        # try to extend neighbor fragment (right)
                        for vj in neighbors_v:
                                fid_vj = self.nodes[vj].parent_fragment
                                fid_vj_obj = self.fragments[fid_vj]
                                # if vj is not the pair of fid
                                # if vj is the tail of the fragment fid_ui
                                # if vj and u are in close proximity
                                if vj != u and fid_vj_obj.tail == vj and abs(fid_vj_obj.start - frag_obj.end) < dist_neighbors:
                                        self.change_fragment_boundry(fid_vj, gp.TAIL, frag_obj.len / 2)

                for (u,v,p) in to_rewire:
                        self.add_edge(MultiGraph.generate_id(),u,v,gp.REWIRE,p)


        # print(
        #       """Fragment id:(tail,head,edge) {}:({},{}.{}) with coverage {} was removed""".format(fid, u, v, e, cov))

        def DFSUtil(self, temp, v, visited):

                # Mark the current vertex as visited
                visited[v] = True

                # Store the vertex to list
                temp.append(v)

                # Repeat for all vertices adjacent
                # to this vertex v
                for n in self.neighbors(v):
                        if not visited[n]:
                                # Update the list
                                temp = self.DFSUtil(temp, n, visited)
                return temp

        def degree(self, u):
                """
                Degree of a node (does not include the edge to its pair (head or tail)
                """
                count = 0
                for v in self.graph[u]:
                        for e in self.graph[u][v]:
                                if self.edges[e].svtype != gp.FRAGMENT:
                                        count += 1
                return count


        def neighbors(self, u):
                """
                Return neighbors for the node
                """
                return self.graph[u]

        def get_next_node(self, u, e):
                """
                Get next node starting from u, connection e
                """
                for v in self.graph[u]:
                        if e in self.graph[u][v]:
                                return v
                return None

        def __collect_edge_info__(self, e, u):
                """
                Collect edge information required for __insert_edges_sorted__
                """
                v = self.get_next_node(u, e)
                v_fid = self.nodes[v].parent_fragment

                data = {"rank1": self.edges[e].weight,
                        "rank2": self.fragments[v_fid].coverage,
                        "rank3": self.edges[e].svtype,
                        "v": v,
                        "u": u
                        }

                return data

        def convert_path2fragments(self, path):
                """
                Convert paths (returned by the DFS search) to fragments list.
                Include orientation of fragments.
                By convention negative fid means inverted fragment fid.

                Arguments:
                        path (list): List of nodes [(v1, e1), (v2, e2)] where e2 represents trasition v1 -> v2

                Returns:
                        List of paths
                        [ [fid1, fid2],
                          [fid1, -fid2],
                        ]
                """
                chain_fragments = []

                for i in range(len(path) - 1):
                        (n1, e1) = path[i]
                        (n2, e2) = path[i + 1]

                        # edge connects head to tail of a fragment
                        if self.edges[e2].svtype == gp.FRAGMENT:

                                fid = self.nodes[n1].parent_fragment

                                # positive strand
                                if self.fragments[fid].tail == n1 and self.fragments[fid].head == n2:
                                        fid = fid
                                elif self.fragments[fid].tail == n2 and self.fragments[fid].head == n1:
                                        # negative strand
                                        fid = -fid
                                else:
                                        raise Exception("""Wrong fragment conformation fid:{} head:{} tail:{} but discovered n1:{} n2:{}""".
                                                        format(fid, self.fragments[fid].head, self.fragments[fid].tail, n1, n2))

                                chain_fragments.append(fid)

                return chain_fragments

        def adjacent_edges_sorted(self, u):
                """
                Get all adjancent edges to a node u.
                Sort the output:
                        - edges weight (descendent)
                        - for same weight consider sorting:
                                - connecting fragment coverage (descending)
                                - if 2 edges connect 2 same nodes, with same weight then priority DEL, DUP, TRA, INV, INS

                Arguments:
                        u (int): None id

                Returns:
                        list edges ids which neighbor vertex u and fulfill the filtering criteria
                """
                available_edges = []

                for n in self.neighbors(u):
                        for k in self.graph[u][n]:
                                # ignore edges connecting head and tail of a fragment
                                if self.edges[k].svtype != gp.FRAGMENT:
                                        data = self.__collect_edge_info__(k, u)
                                        cycles.insert_edges_sorted(available_edges, k, data)

                return available_edges

        def connected_components(self):
                """
                Compute connected components
                """
                visited = {v: False for v in self.nodes}
                cc = []
                for v in self.nodes:
                        if not visited[v]:
                                temp = []
                                arr = self.DFSUtil(temp, v, visited)
                                cc.append(arr)
                return cc

        def subgraph(self, list_vertices):
                """
                Extracts the subgraph based on a list of nodes
                Returns:
                        the subgraph which contains the nodes listed in vertices
                """
                new_subgraph = deepcopy(self)
                # nodes
                for v in self.get_nodes():
                        if v not in list_vertices:
                                new_subgraph.remove_node(v)
                return new_subgraph

        def update_pairs(self):
                """
                Go over each fragment and match the node pairs
                """
                for fid, frag in self.fragments.items():
                        n1 = frag.tail
                        n2 = frag.head
                        self.nodes[n1].idpair = n2
                        self.nodes[n2].idpair = n1

        def reset_time(self):
                """
                Reset time for DFS tree search
                """
                self.start_time = self.end_time = 0

                for n in self.nodes:
                        self.nodes[n].start_time = self.fnodes[n].end_time = 0

        def __find_cycle__(self, node, stack):
                """
                Find circular path by popping out the stack

                Arguments:
                        node (decoil.graph.Node):
                        stack (list): Stack
                """
                stack_c = deepcopy(stack)
                cycle = []

                #         print("Find cycle for " , node)
                #         print("Current stack: ", stack)

                while stack_c:
                        # pop always 2xnodes (which represent one fragment)
                        v, k = stack_c.pop()  # v=node_id, k=edge_id
                        cycle.append((v, k))
                        v, k = stack_c.pop()
                        cycle.append((v, k))
                        if v == node and len(cycle) > 1:
                                #                 print("Found cycle",cycle)
                                #                 print()
                                # reverse ordering
                                return cycle[::-1]
                return None

        def __stack_to_list__(self, stack):
                stack_c = deepcopy(stack)
                cycle = []

                while stack_c:
                        # pop always 2xnodes (which represent one fragment)
                        v, k = stack_c.pop()  # v=node_id, k=edge_id
                        cycle.append((v, k))
                        v, k = stack_c.pop()
                        cycle.append((v, k))

                return cycle[::-1]

        def __visit_fragment__(self, u, u_pair, e, fid, visited_nodes, visited_edges, visited_fragments, stack):
                """
                Arguments:
                        u (int): Next node id
                        u_pair (int): Pair of u belonging to same fragment
                        e (int): Edge connecting fragments
                        fid (int): Fragment ID
                        visited_nodes (dict): records the visited nodes
                        visited_edges (dict): records the visited edges
                        visited_fragments (dict): records the visited fragments
                        stack (list): stack tracking current path (node, incoming_edge)
                """
                self.nodes[u].start_time = self.start_time
                self.start_time += 1
                self.nodes[u_pair].start_time = self.start_time
                self.start_time += 1

                # mark edge as visited
                if e:
                        visited_edges[e] = visited_edges[e] + 1

                # mark fragment as visited
                visited_fragments[fid] = visited_fragments[fid] + 1

                # mark head tail of fragment as visited
                visited_nodes[u] = visited_nodes[u] + 1
                visited_nodes[u_pair] = visited_nodes[u_pair] + 1

                # mark edge connecting head to tail as visited
                e_con = self.fragments[fid].edge
                visited_edges[e_con] = visited_edges[e_con] + 1

                stack.append((u, e))
                stack.append((u_pair, e_con))

        def __unvisit_fragment__(self, fid, visited_edges, visited_nodes, visited_fragments, stack):
                """
                Undo states

                Arguments:
                        fid (int): Fragment ID
                        visited_nodes (dict): records the visited nodes
                        visited_edges (dict): records the visited edges
                        visited_fragments (dict): records the visited fragments
                        stack (list): stack tracking current path
                """

                # pop one part of the fragment
                _, _ = stack.pop()
                # pot the other part of the fragment
                # e - incoming edge for node u
                u, e = stack.pop()

                # reset edge e (check if edge is not none)
                if e:
                        visited_edges[e] = visited_edges[e] - 1

                # reset fragment
                u, v, e = self.fragments[fid].tail, self.fragments[fid].head, self.fragments[fid].edge
                visited_nodes[v] = visited_nodes[v] - 1
                visited_nodes[u] = visited_nodes[u] - 1
                visited_edges[e] = visited_edges[e] - 1
                visited_fragments[fid] = visited_fragments[fid] - 1

                self.nodes[u].start_time = self.nodes[u].end_time = 0
                self.nodes[v].start_time = self.nodes[v].end_time = 0

        def __find3__(self, node, visited_edges, visited_nodes, visited_fragments, stack, paths):
                """
                Find DFS tree, consider unique use of edges but allow multiple use of fragments
                """
                # DFS recursion
                available_edges = self.adjacent_edges_sorted(node)
                # print("forward", node, stack)

                for e in available_edges:
                        k = e[0]
                        v = e[1]["v"]
                        if visited_edges[k] == 0 and visited_nodes[v] > 0:
                                newpath = self.__find_cycle__(v, stack)
                                paths = self.__add_cycle__(paths, newpath)
                                # print("Circular paths", paths)

                for e in available_edges:
                        k = e[0]
                        v = e[1]["v"]
                        fid = self.nodes[v].parent_fragment

                        if visited_edges[k] == 0:
                                # mark node visited
                                v_pair = self.nodes[v].idpair
                                self.__visit_fragment__(v,
                                                        v_pair,
                                                        k,
                                                        fid,
                                                        visited_nodes,
                                                        visited_edges,
                                                        visited_fragments,
                                                        stack
                                                        )
                                self.__find3__(v_pair, visited_edges, visited_nodes, visited_fragments, stack, paths)


                (k, e) = stack.pop()
                self.nodes[k].end_time = self.end_time
                self.end_time += 1
                (k, e) = stack.pop()
                self.nodes[k].end_time = self.end_time
                self.end_time += 1
                # print("backtrack", node, stack)

        def __find__(self, node, visited_edges, visited_nodes, visited_fragments, stack, paths):
                """
                Find DFS tree
                """
                # DFS recursion
                available_edges = self.adjacent_edges_sorted(node)
                # print("forward", node, stack)

                for e in available_edges:
                        k = e[0]
                        v = e[1]["v"]
                        fid = self.nodes[v].parent_fragment
                        if visited_edges[k] == 0 and visited_nodes[v] == 0:
                                # mark node visited
                                v_pair = self.nodes[v].idpair
                                self.__visit_fragment__(v,
                                                        v_pair,
                                                        k,
                                                        fid,
                                                        visited_nodes,
                                                        visited_edges,
                                                        visited_fragments,
                                                        stack
                                                        )
                                self.__find__(v_pair, visited_edges, visited_nodes, visited_fragments, stack, paths)

                # detect backedges
                for e in available_edges:
                        k = e[0]
                        v = e[1]["v"]
                        if visited_edges[k] == 0 and visited_nodes[v] > 0:
                                newpath = self.__find_cycle__(v, stack)
                                paths = self.__add_cycle__(paths, newpath)
                                # if newpath:
                                #       paths.append(newpath)

                # pop fragment=2xnodes
                (k, e) = stack.pop()
                self.nodes[k].end_time = self.end_time
                self.end_time += 1
                (k, e) = stack.pop()
                self.nodes[k].end_time = self.end_time
                self.end_time += 1
                # print("backtrack", node, stack)

        def __add_cycle__(self, listpaths, newpath):
                """
                Append cycle to list of found cycles if not exists.
                Does not take circularity into account
                """

                # 0. ignore entry if empty
                if not newpath or len(newpath) == 0:
                        return listpaths

                # 1. convert path to fragments chain (keep orientation)
                #    negative ids mean inverted fragment
                chain_fragments = self.convert_path2fragments(newpath)

                # 2. create and add hash
                added = self.add_to_hash(chain_fragments)

                # 3. new entry
                if added:
                        listpaths.append(chain_fragments)
                return listpaths

        def __find2__(self, node, visited_edges, visited_nodes, visited_fragments, stack, paths):
                """
                Find DFS tree
                """
                # DFS recursion
                available_edges = self.adjacent_edges_sorted(node)
                fid_node = self.nodes[node].parent_fragment
                # print("forward", node, stack)

                # detect backedges
                for e in available_edges:
                        k = e[0]
                        v = e[1]["v"]
                        if visited_edges[k] == 0 and visited_nodes[v] > 0:
                                newpath = self.__find_cycle__(v, stack)
                                paths = self.__add_cycle__(paths, newpath)

                for e in available_edges:
                        # edge id
                        k = e[0]
                        # next vertex
                        v = e[1]["v"]
                        fid = self.nodes[v].parent_fragment
                        if visited_edges[k] == 0 and visited_nodes[v] == 0:
                                # mark node visited
                                v_pair = self.nodes[v].idpair
                                self.__visit_fragment__(v,
                                                        v_pair,
                                                        k,
                                                        fid,
                                                        visited_nodes,
                                                        visited_edges,
                                                        visited_fragments,
                                                        stack
                                                        )
                                self.__find2__(v_pair, visited_edges, visited_nodes, visited_fragments, stack, paths)

                self.__unvisit_fragment__(fid_node,
                                          visited_edges,
                                          visited_nodes,
                                          visited_fragments,
                                          stack)
                # print("backtrack", node, stack)

        def find_simple_circles(self, start):
                """
                Find all simple circles and add the unique candidates to self.paths
                
                Returns:
                        [path1 [(n1, None),(n2, e1)], which means a transition from n1 to n2 over edge e1
                        path2 []]
                """
                stack = []
                paths = []
                visited_nodes = {n: 0 for n in self.nodes}
                visited_edges = {e: 0 for e in self.edges}
                visited_fragments = {f: 0 for f in self.fragments}

                # add head and tail of the fragment
                self.__visit_fragment__(start,
                                        self.nodes[start].idpair,
                                        None,
                                        self.nodes[start].parent_fragment,
                                        visited_nodes,
                                        visited_edges,
                                        visited_fragments,
                                        stack)

                # find simple circles
                self.__find__(self.nodes[start].idpair, visited_edges, visited_nodes, visited_fragments, stack, paths)

                # add deduplicated candidates to the larger list
                return paths

        def find_simple_circles2(self, start):
                """
                Find all simple circles and add the unique candidates to self.paths

                Returns:
                        [path1 [(n1, None),(n2, e1)], which means a transition from n1 to n2 over edge e1
                        path2 []]
                """
                stack = []
                paths = []
                visited_nodes = {n: 0 for n in self.nodes}
                visited_edges = {e: 0 for e in self.edges}
                visited_fragments = {f: 0 for f in self.fragments}

                # add head and tail of the fragment
                self.__visit_fragment__(start,
                                        self.nodes[start].idpair,
                                        None,
                                        self.nodes[start].parent_fragment,
                                        visited_nodes,
                                        visited_edges,
                                        visited_fragments,
                                        stack)

                # find simple circles
                self.__find2__(self.nodes[start].idpair, visited_edges, visited_nodes, visited_fragments, stack, paths)

                return paths

        def find_simple_circles3(self, start):
                """
                Find all simple circles and add the unique candidates to self.paths

                Returns:
                        [path1 [(n1, None),(n2, e1)], which means a transition from n1 to n2 over edge e1
                        path2 []]
                """
                stack = []
                paths = []
                visited_nodes = {n: 0 for n in self.nodes}
                visited_edges = {e: 0 for e in self.edges}
                visited_fragments = {f: 0 for f in self.fragments}

                # add head and tail of the fragment
                self.__visit_fragment__(start,
                                        self.nodes[start].idpair,
                                        None,
                                        self.nodes[start].parent_fragment,
                                        visited_nodes,
                                        visited_edges,
                                        visited_fragments,
                                        stack)

                # find simple circles
                self.__find3__(self.nodes[start].idpair, visited_edges, visited_nodes, visited_fragments, stack, paths)

                # add deduplicated candidates to the larger list
                return paths


        def add_to_hash(self, path):
                """
                Add to metahash the path

                Arguments:
                        path (list): List of fragments composing a cycle

                Returns:
                        True/False (if added or not to self.hash_circles)
                """

                h = cycles.get_hash(path)

                cpath = self.get_canonical(path)
                ch = cycles.get_hash(cpath)

                added = False
                if h not in self.hash_circles:
                        added = True
                        # map unique hash to canonical hash of the cycle
                        self.hash_circles[h] = ch
                        if ch not in self.hash_canonical:
                                self.hash_canonical[ch] = cpath

                return added

        def remove_from_hash(self, id):
                del self.hash_canonical[id]

        def get_canonical(self, path):
                """
                Find canonical for a cycle
                """
                return cycles.convert2canonical_path(path)

        def __str__(self):
                _str = "\nNodes:\n"
                for node in self.nodes:
                        _str += self.nodes[node].__str__() + "\n"

                _str += "\nEdges:\n"
                for edges in self.edges:
                        _str += self.edges[edges].__str__() + "\n"

                _str += "\nEdges links:\n"
                for u in dict(self.graph.items()):
                        for v in self.graph[u]:
                                for edgeid in self.graph[u][v]:
                                        _str += str(edgeid) + ":" + str(u) + "-" + str(v) + "\n"

                _str += "\nFragments:\n"
                for f in self.fragments:
                        _str += str(f) + ":" + str(self.fragments[f]) + "\n"

                return _str

        def print_paths(self):
                """
                Print all circular paths
                """
                for p in self.paths:
                        print(p)

        def get_nodes(self):
                return dict(self.nodes.items())

        def get_edges(self):
                return dict(self.edges.items())

        def get_fragment_intervals(self):
                return self.fragments_intervals

        def get_fragment_intervals_bychr(self, chr):
                return self.fragments_intervals.tree_fragments[chr]

        def get_fragment(self, fid):
                return dict(self.fragments.items())[fid]

        def get_fragments(self):
                return self.fragments

        def set_nodes(self, nodes):
                self.nodes = nodes

        def set_edges(self, edges):
                self.edges = edges

        def set_fragment_intervals(self, fi):
                self.fragments_intervals = fi

        def get_hash_circles(self):
                return self.hash_circles

        # def get_metahash_circles(self):
        #       return self.metahash_circles

        def get_hash_canonical(self):
                return self.hash_canonical

        @staticmethod
        def load_graph(file):
                """
                Load a graph based on the specification file. Only for test
                Example:
                        # test file for ABCBCBBA
                        # part -> {head=1, tail=0}
                        # svtype -> {frag=0, DUP=1, DEL=2, INV=3, BND=4, INS=5, spatial=6}
                        #node1  node2   frag1   frag2   part1   part2   edgeid  svtype cov
                        1       2       A       A       0       1       0       0   100
                        2       3       A       B       1       0       1       6   40
                        3       4       B       B       0       1       2       0
                        4       5       B       C       1       0       3       6
                        5       6       C       C       0       1       4       0
                        1       2       A       A       0       1       5       1
                        3       4       B       B       0       1       6       1
                        1       4       A       B       0       1       7       1
                        3       6       B       C       0       1       8       1
                        
                Arguments:
                        file (str): graph specification file
                
                Returns:
                        graph (MultiGraph)
                """
                graph = MultiGraph()

                with open(file, "r") as f:
                        for l in f:
                                if l.startswith("#"):
                                        continue
                                else:
                                        (node1, node2, frag1, frag2, part1, part2, edgeid, svtype, cov) = l.strip().split("\t")
                                        node1 = int(node1)
                                        node2 = int(node2)
                                        part1 = int(part1)
                                        part2 = int(part2)
                                        edgeid = int(edgeid)
                                        svtype = int(svtype)
                                        cov = float(cov)

                                        # insert fragment
                                        if frag1 == frag2 and frag1 not in graph.get_fragments():

                                                id1 = node1
                                                id2 = node2
                                                id3 = edgeid

                                                props1 = {gnp.CHR: "chr1", gnp.POS: 1, gnp.COLOR: "black"}
                                                props2 = {gnp.CHR: "chr1", gnp.POS: 1, gnp.COLOR: "black"}

                                                graph.add_node(id1, part1, frag1, props1)
                                                graph.add_node(id2, part2, frag1, props2)

                                                frag_len = (int(2) - int(1))
                                                props3 = {gep.CHR: chr,
                                                          gep.START: 1,
                                                          gep.END: 2,
                                                          gep.COLOR: "gray",
                                                          gep.SVTYPE: gp.FRAGMENT,
                                                          gep.FRAGMENT_LEN: frag_len,
                                                          gep.FRAGMENT_NAME: frag1,
                                                          gep.WEIGHT: cov}
                                                graph.add_edge(id3, id1, id2, gp.FRAGMENT, props3)

                                                graph.add_fragment(frag1, id1, id2, id3, {fg.CHR: "chr1",
                                                                                          fg.START: 1,
                                                                                          fg.END: 2,
                                                                                          fg.COVERAGE: cov
                                                                                          })
                                                graph.get_fragment_intervals().add_fragment("chr1", 1, 2, {fg.FRAGID: frag1,
                                                                                                           fg.EDGEID: id3,
                                                                                                           fg.NODE_TAIL: id1,
                                                                                                           fg.NODE_HEAD: id2
                                                                                                           })

                                        else:
                                                graph.add_edge(edgeid, node1, node2, svtype, {"total_cov": cov,
                                                                                              "len": 10,
                                                                                              "color": "black",
                                                                                              "weight": cov})

                        # match nodes belonging to the same fragment
                        graph.update_pairs()

                return graph

        @staticmethod
        def generate_id():
                MultiGraph.id_generator += 1
                return MultiGraph.id_generator
