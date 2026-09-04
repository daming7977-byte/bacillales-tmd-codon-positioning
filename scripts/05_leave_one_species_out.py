from pathlib import Path
from collections import defaultdict, Counter

RBH_DIR = Path("work/orthology/rbh")

rbh_files = sorted(RBH_DIR.glob("*.tsv"))

species = sorted({
    x
    for p in rbh_files
    for x in p.stem.split("__")
})

print("Species/accessions:")
for s in species:
    print(" ", s)

def analyze(excluded=None):

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

        ra = find(a)
        rb = find(b)

        if ra == rb:
            return

        if rank[ra] < rank[rb]:
            ra, rb = rb, ra

        parent[rb] = ra

        if rank[ra] == rank[rb]:
            rank[ra] += 1

    for path in rbh_files:

        acc1, acc2 = path.stem.split("__", 1)

        if excluded is not None and (
            acc1 == excluded or acc2 == excluded
        ):
            continue

        with path.open() as f:
            for line in f:

                if not line.strip():
                    continue

                fields = line.rstrip("\n").split("\t")

                q = fields[0]
                t = fields[1]

                union(
                    f"{acc1}|{q}",
                    f"{acc2}|{t}"
                )

    components = defaultdict(list)

    for node in parent:
        components[find(node)].append(node)

    strict = []

    for nodes in components.values():

        spp = [
            node.split("|", 1)[0]
            for node in nodes
        ]

        counts = Counter(spp)

        if max(counts.values()) == 1:
            strict.append(nodes)

    n_species = len(species) - (1 if excluded else 0)

    full = sum(
        len(g) == n_species
        for g in strict
    )

    minus1 = sum(
        len(g) >= n_species - 1
        for g in strict
    )

    return len(strict), full, minus1


total, full, near = analyze()

print()
print("FULL 8-species panel")
print(f"Strict groups total: {total}")
print(f"8/8 groups: {full}")
print(f">=7/8 groups: {near}")

print()
print("LEAVE-ONE-SPECIES-OUT")

for s in species:

    total, full, near = analyze(s)

    print(
        f"{s}\t"
        f"strict_total={total}\t"
        f"7/7={full}\t"
        f">=6/7={near}"
    )
