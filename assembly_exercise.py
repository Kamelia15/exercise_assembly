def build_adjacency_list(edges, directed=False):
    graph = {}
    for u, v in edges:
        graph.setdefault(u, [])
        graph.setdefault(v, [])
        graph[u].append(v)
        if not directed:
            graph[v].append(u)
    return graph


def degrees(graph):
    return {v: len(neighbors) for v, neighbors in graph.items()}


def in_out_degree(graph):
    out_degree = degrees(graph)
    in_degree = {v: 0 for v in graph}
    for neighbors in graph.values():
        for v in neighbors:
            in_degree[v] += 1
    return in_degree, out_degree


def dfs(graph, start):
    visited = [start]
    seen = {start}
    stack = [(start, iter(sorted(graph.get(start, []))))]
    while stack:
        v, neighbors = stack[-1]
        try:
            n = next(neighbors)
            if n not in seen:
                seen.add(n)
                visited.append(n)
                stack.append((n, iter(sorted(graph.get(n, [])))))
        except StopIteration:
            stack.pop()
    return visited


def connected_components(graph):
    components = []
    visited_all = set()
    for v in sorted(graph):
        if v not in visited_all:
            component = set(dfs(graph, v))
            components.append(component)
            visited_all.update(component)
    return components


def is_weakly_connected(graph):
    undirected = {v: set() for v in graph}
    for u, neighbors in graph.items():
        for v in neighbors:
            undirected.setdefault(v, set())
            undirected[u].add(v)
            undirected[v].add(u)
    graph_undirected = {v: list(n) for v, n in undirected.items()}
    return len(connected_components(graph_undirected)) <= 1


def eulerian_type(graph, directed=False):
    active = {u for u, neighbors in graph.items() if neighbors}
    for neighbors in graph.values():
        active.update(neighbors)

    if not active:
        return "circuit"

    undirected = {v: set() for v in active}
    for u, neighbors in graph.items():
        for v in neighbors:
            undirected[u].add(v)
            undirected[v].add(u)

    visited = set()
    stack = [next(iter(active))]
    while stack:
        v = stack.pop()
        if v not in visited:
            visited.add(v)
            stack.extend(undirected[v] - visited)

    if visited != active:
        return "none"

    if directed:
        in_deg, out_deg = in_out_degree(graph)
        starts = ends = 0
        for v in active:
            diff = out_deg[v] - in_deg[v]
            if diff == 1:
                starts += 1
            elif diff == -1:
                ends += 1
            elif diff != 0:
                return "none"
        if starts == 0 and ends == 0:
            return "circuit"
        if starts == 1 and ends == 1:
            return "path"
        return "none"

    deg = degrees(graph)
    odd = sum(deg[v] % 2 for v in active)
    if odd == 0:
        return "circuit"
    if odd == 2:
        return "path"
    return "none"


def kmer_composition(text, k):
    return [text[i:i + k] for i in range(len(text) - k + 1)]


def debruijn_graph(kmers):
    graph = {}
    for km in kmers:
        prefix, suffix = km[:-1], km[1:]
        graph.setdefault(prefix, []).append(suffix)
        graph.setdefault(suffix, [])
    return graph


def spell_path(nodes):
    return nodes[0] + "".join(node[-1] for node in nodes[1:])


genome = "CTCGGACATTAATGC"
kmers = kmer_composition(genome, 3)
graph = debruijn_graph(kmers)

assert len(kmers) == len(genome) - 3 + 1
assert sum(len(v) for v in graph.values()) == len(kmers)
assert set(graph) == {km[:-1] for km in kmers} | {km[1:] for km in kmers}
assert spell_path(kmers) == genome
assert eulerian_type(graph, directed=True) == "path"

for text in ["GATTACA", "CATCATCAT", "ATATAT", genome]:
    for k in range(1, len(text) + 1):
        assert spell_path(kmer_composition(text, k)) == text

repeat_kmers = kmer_composition("CATCATCAT", 3)
repeat_graph = debruijn_graph(repeat_kmers)

assert repeat_graph["CA"].count("AT") == 3
assert sum(len(v) for v in repeat_graph.values()) == len(repeat_kmers)
assert eulerian_type(repeat_graph, directed=True) == "path"

circuit_graph = debruijn_graph(kmer_composition("ATATAT", 3))
assert eulerian_type(circuit_graph, directed=True) == "circuit"

print("Genome:", genome)
print("3-mers:", kmers)
print("De Bruijn graph:", graph)
print("In-degrees, out-degrees:", in_out_degree(graph))
print("Eulerian type:", eulerian_type(graph, directed=True))
print("Spelled genome:", spell_path(kmers))
