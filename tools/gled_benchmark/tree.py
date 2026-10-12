"""UPGMA on GLED cognate distances; Newick and forms for clepsydra."""


def upgma(names, dist):
    """names: list; dist(a, b) -> float. Returns nested (left, right, height)
    tuples with leaves as names; ultrametric. Missing distances take the mean."""
    pairs = [(a, b) for i, a in enumerate(names) for b in names[i + 1:]]
    vals = [dist(a, b) for a, b in pairs]
    known = [v for v in vals if v is not None]
    fill = sum(known) / len(known) if known else 1.0
    d = {}
    for (a, b), v in zip(pairs, vals):
        d[(a, b)] = d[(b, a)] = fill if v is None else v
    clusters = {n: (n, 1, 0.0) for n in names}      # key -> (tree, size, height)
    keys = list(names)
    D = dict(d)
    k = 0
    while len(keys) > 1:
        best = min(((a, b) for i, a in enumerate(keys) for b in keys[i + 1:]),
                   key=lambda p: (D[p], p))
        a, b = best
        h = D[best] / 2.0
        ta, sa, _ = clusters[a]
        tb, sb, _ = clusters[b]
        k += 1
        new = f"_u{k}"
        for c in keys:
            if c in (a, b):
                continue
            D[(new, c)] = D[(c, new)] = (D[(a, c)] * sa + D[(b, c)] * sb) / (sa + sb)
        clusters[new] = ((ta, tb, h), sa + sb, h)
        keys = [c for c in keys if c not in (a, b)] + [new]
    return clusters[keys[0]][0]


def leaves(t):
    return [t] if isinstance(t, str) else leaves(t[0]) + leaves(t[1])


def newick(tree, root="root", clades=None):
    """Labelled Newick with unit branch lengths implied by heights (relative)."""
    counter = [0]

    def h(t):
        return 0.0 if isinstance(t, str) else t[2]

    def rec(t, parent_h):
        if isinstance(t, str):
            return f"{t}:{parent_h:.6f}"
        counter[0] += 1
        lab = root if parent_h is None else f"n{counter[0]}"
        l = rec(t[0], t[2] - h(t[0]))
        r = rec(t[1], t[2] - h(t[1]))
        s = f"({l},{r}){lab}"
        if clades is not None:
            clades[lab] = frozenset(leaves(t))
        return s if parent_h is None else f"{s}:{parent_h:.6f}"
    # label the root last so n-numbers are stable
    body = rec(tree, None)
    return body + ";"
