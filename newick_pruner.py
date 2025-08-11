from Bio import Phylo
inp="../agora/diptera.pruned.newick"
out="../agora/bibio_trio.newick"
keep={"Bibio_marci","Dilophus_febrilis","Plecia_longiforceps"}

t = Phylo.read(inp, "newick")          # parse the tree
for tip in list(t.get_terminals()):    # iterate all leaves
    if tip.name not in keep:           # if not one of the trio
        t.prune(tip)                   # remove it (updates tree)

Phylo.write(t, out, "newick")          # write the pruned tree


