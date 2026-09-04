from pathlib import Path
from collections import defaultdict, Counter

RBH_DIR = Path("work/orthology/rbh_conservative10")
OUTDIR = Path("work/orthology/strict_groups_conservative10")
OUTDIR.mkdir(parents=True, exist_ok=True)

parent = {}
rank = {}

def make(x):
    if x not in parent:
        parent[x] = x
        rank[x] = 0

def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

def union(a, b):
    make(a)
    make(b)
    ra, rb = find(a), find(b)
    if ra == rb:
        return
    if rank[ra] < rank[rb]:
        ra, rb = rb, ra
    parent[rb] = ra
    if rank[ra] == rank[rb]:
        rank[ra] += 1

rbh_files = sorted(RBH_DIR.glob("*.tsv"))
print(f"RBH files: {len(rbh_files)}")

edge_count = 0

for path in rbh_files:
    stem = path.stem
    acc1, acc2 = stem.split("__", 1)

    with path.open() as f:
        for line in f:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            query = fields[0]
            target = fields[1]

            node1 = f"{acc1}|{query}"
            node2 = f"{acc2}|{target}"

            union(node1, node2)
            edge_count += 1

print(f"RBH edges: {edge_count}")

components = defaultdict(list)
for node in parent:
    components[find(node)].append(node)

print(f"Connected components with RBH edges: {len(components)}")

strict_groups = []

for nodes in components.values():
    species = [node.split("|", 1)[0] for node in nodes]
    counts = Counter(species)

    # strict one-to-one:
    # no species may contribute more than one protein
    if max(counts.values()) == 1:
        strict_groups.append(sorted(nodes))

coverage_counts = Counter(len(group) for group in strict_groups)

print(f"Strict one-to-one groups: {len(strict_groups)}")
print("Species coverage distribution:")
for n in sorted(coverage_counts, reverse=True):
    print(f"  {n}/10 species: {coverage_counts[n]}")

for cutoff in (8, 7, 6):
    n = sum(1 for g in strict_groups if len(g) >= cutoff)
    print(f"Strict groups >= {cutoff}/10 species: {n}")

outfile = OUTDIR / "strict_one_to_one_groups.tsv"

with outfile.open("w") as out:
    out.write("orthogroup\tspecies_count\tmembers\n")

    groups_sorted = sorted(
        strict_groups,
        key=lambda g: (-len(g), g)
    )

    for i, group in enumerate(groups_sorted, start=1):
        out.write(
            f"orthogroup_{i}\t{len(group)}\t"
            + ";".join(group)
            + "\n"
        )

print(f"Output: {outfile}")
