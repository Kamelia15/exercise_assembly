from collections import Counter
from heapq import heapify, heappop
from Bio import SeqIO


#________________________TASK_1________________________

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
    vertices = set(graph)
    for neighbors in graph.values():
        vertices.update(neighbors)

    out_degree = {v: len(graph.get(v, [])) for v in vertices}
    in_degree = dict.fromkeys(vertices, 0)

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

    components = connected_components(
        {v: list(n) for v, n in undirected.items()}
    )
    return len(components) <= 1


def eulerian_type(graph, directed=False):
    active = set()

    for u, neighbors in graph.items():
        if neighbors:
            active.add(u)
            active.update(neighbors)

    if not active:
        return "circuit"

    undirected = {v: set() for v in active}

    for u, neighbors in graph.items():
        for v in neighbors:
            undirected[u].add(v)
            undirected[v].add(u)

    visited = set()
    stack = [min(active)]

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

        if starts == ends == 0:
            return "circuit"
        if starts == ends == 1:
            return "path"
        return "none"

    deg = degrees(graph)
    odd = sum(deg.get(v, 0) % 2 for v in active)

    if odd == 0:
        return "circuit"
    if odd == 2:
        return "path"
    return "none"


#________________________TASK_2________________________

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

for text in ["GATTACA", "CATCATCAT", genome]:
    for k in range(1, len(text) + 1):
        assert spell_path(kmer_composition(text, k)) == text

repeat_kmers = kmer_composition("CATCATCAT", 3)
repeat_graph = debruijn_graph(repeat_kmers)

assert repeat_graph["CA"].count("AT") == 3
assert sum(len(v) for v in repeat_graph.values()) == len(repeat_kmers)


#________________________TASK_3________________________

def edge_counter(graph):
    return Counter(
        (u, v)
        for u, neighbors in graph.items()
        for v in neighbors
    )


def hierholzer(graph, start, verbose=False):
    nodes = set(graph)
    for neighbors in graph.values():
        nodes.update(neighbors)

    unused = {v: list(graph.get(v, [])) for v in nodes}
    for neighbors in unused.values():
        heapify(neighbors)

    remaining = sum(len(neighbors) for neighbors in graph.values())
    walk_number = 1

    def walk(first):
        nonlocal remaining
        path = [first]
        v = first

        while unused[v]:
            v = heappop(unused[v])
            remaining -= 1
            path.append(v)

        return path

    path = walk(start)

    if verbose:
        print("Walk 1:", path)

    while remaining:
        i = next(i for i, v in enumerate(path) if unused[v])
        restart = path[i]
        new_walk = walk(restart)
        walk_number += 1

        if verbose:
            print(f"Walk {walk_number}, restart at {restart}:", new_walk)

        path = path[:i] + new_walk + path[i + 1:]

        if verbose:
            print("Spliced:", path)

    return path


def eulerian_path(graph, verbose=False):
    if not graph:
        return []

    in_deg, out_deg = in_out_degree(graph)
    starts = [
        v for v in graph
        if out_deg[v] - in_deg[v] == 1
    ]

    start = starts[0] if starts else min(graph)
    return hierholzer(graph, start, verbose)


original_graph = {u: ns.copy() for u, ns in graph.items()}
walk = eulerian_path(graph)
spelled = spell_path(walk)

assert len(walk) == len(kmers) + 1
assert Counter(zip(walk, walk[1:])) == edge_counter(graph)
assert Counter(kmer_composition(spelled, 3)) == Counter(kmers)
assert graph == original_graph
assert spelled == genome
assert eulerian_path(debruijn_graph(kmers[::-1])) == walk

assert eulerian_type(graph, directed=True) == "path"
assert eulerian_type(repeat_graph, directed=True) == "path"

print("Genome:", genome)
print("3-mers:", kmers)
print("De Bruijn graph:", graph)
print("In/out degrees:", in_out_degree(graph))
print("Eulerian type:", eulerian_type(graph, directed=True))
print("Eulerian path:", walk)
print("Spelled genome:", spelled)


#________________________TASK_4________________________

def contigs(graph):
    in_deg, out_deg = in_out_degree(graph)
    vertices = set(in_deg) | set(out_deg)
    one = {
        v for v in vertices
        if in_deg[v] == 1 and out_deg[v] == 1
    }

    adj = {v: sorted(graph.get(v, [])) for v in vertices}
    used = {v: [False] * len(adj[v]) for v in vertices}
    result = []

    for u in sorted(vertices):
        if u in one:
            continue

        for i, v in enumerate(adj[u]):
            if used[u][i]:
                continue

            used[u][i] = True
            path = [u, v]

            while v in one:
                j = used[v].index(False)
                used[v][j] = True
                v = adj[v][j]
                path.append(v)

            result.append(spell_path(path))

    for u in sorted(vertices):
        for i in range(len(adj[u])):
            if used[u][i]:
                continue

            path = [u]
            v = u

            while True:
                j = used[v].index(False)
                used[v][j] = True
                v = adj[v][j]
                path.append(v)

                if v == u:
                    break

            result.append(spell_path(path))

    return sorted(result)


def trim_cycle(text, k):
    overlap = k - 1
    if text[:overlap] == text[-overlap:]:
        return text[:-overlap]
    return text


for text in [genome, "CATCATCAT", "ACGTGCA", "ACGTTGCAAC"]:
    g = debruijn_graph(kmer_composition(text, 3))
    cs = contigs(g)
    combined = Counter()

    for c in cs:
        combined.update(
            edge_counter(debruijn_graph(kmer_composition(c, 3)))
        )

    assert combined == edge_counter(g)

unique_text = "ACGTGCA"
assert contigs(debruijn_graph(kmer_composition(unique_text, 3))) == [unique_text]

cycle_text = "ACGTTGCAAC"
cycle_contigs = contigs(debruijn_graph(kmer_composition(cycle_text, 3)))
assert len(cycle_contigs) == 1

trimmed_cycle = trim_cycle(cycle_contigs[0], 3)
cycle_core = cycle_text[:-2]

assert len(trimmed_cycle) == len(cycle_core)
assert trimmed_cycle in cycle_core + cycle_core

print("Contigs:", contigs(graph))

#________________________TASK_5________________________

reads = [
    str(record.seq)
    for record in SeqIO.parse("reads.fasta", "fasta")
]

assemblies = {}

print("k |V| |E| type contigs longest_contig spelled_length")

for k in [9, 12, 13, 14, 21]:
    kmers = sorted({
        read[i:i + k]
        for read in reads
        for i in range(len(read) - k + 1)
    })

    graph = debruijn_graph(kmers)
    typ = eulerian_type(graph, directed=True)
    cs = contigs(graph)

    combined = Counter()
    for c in cs:
        combined.update(
            edge_counter(debruijn_graph(kmer_composition(c, k)))
        )

    assert combined == edge_counter(graph)

    trimmed_contigs = [trim_cycle(c, k) for c in cs]
    longest = max(map(len, trimmed_contigs), default=0)
    assembled = None

    if typ != "none":
        walk = eulerian_path(graph)
        raw = spell_path(walk)

        assert len(walk) == len(kmers) + 1
        assert len(raw) == len(kmers) + k - 1
        assert Counter(zip(walk, walk[1:])) == edge_counter(graph)
        assert set(kmer_composition(raw, k)) == set(kmers)

        assembled = trim_cycle(raw, k) if walk[0] == walk[-1] else raw
        assemblies[k] = assembled

    print(
        k, len(graph), len(kmers), typ, len(cs),
        longest, len(assembled) if assembled is not None else "—"
    )

    if assembled is not None:
        print(assembled)

keys = sorted(assemblies)

for i, k1 in enumerate(keys):
    for k2 in keys[i + 1:]:
        s, t = assemblies[k1], assemblies[k2]
        same = len(s) == len(t) and s in t + t
        print("Same circle", k1, k2, same)
