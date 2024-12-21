class UnionFind:
    def __init__(self, size):
        self.parent = list(range(size))
        self.rank = [0] * size

    def find(self, p):
        if self.parent[p] != p:
            self.parent[p] = self.find(self.parent[p])  # Path compression
        return self.parent[p]

    def union(self, p, q):
        rootP = self.find(p)
        rootQ = self.find(q)
        if rootP != rootQ:
            # Union by rank
            if self.rank[rootP] > self.rank[rootQ]:
                self.parent[rootQ] = rootP
            elif self.rank[rootP] < self.rank[rootQ]:
                self.parent[rootP] = rootQ
            else:
                self.parent[rootQ] = rootP
                self.rank[rootP] += 1

    def connected(self, p, q):
        return self.find(p) == self.find(q)

    def component_count(self):
        return len(set(self.find(x) for x in range(len(self.parent))))


def find_connected_components(graph):
    # Initialize Union-Find for the graph's nodes
    uf = UnionFind(max(graph) + 1)  # Assumes nodes are labeled from 0 to max

    # Union nodes based on edges
    for node, neighbors in graph.items():
        for neighbor in neighbors:
            uf.union(node, neighbor)

    # Extract the components
    component_map = {}
    for node in graph:
        root = uf.find(node)
        if root not in component_map:
            component_map[root] = []
        component_map[root].append(node)

    return list(component_map.values())

# Define the graph
graph = {
    4: [75, 103],
    40: [171],
    53: [],
    75: [4],
    103: [4],
    171: [40],
    220: [],
    229: [],
    261: [327],
    327: [261],
    340: [343],
    343: [340],
    388: []
}


# Calculate connected components
components = find_connected_components(graph)
print("Connected components:", components)

