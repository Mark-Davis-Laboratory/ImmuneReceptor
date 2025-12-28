using Graphs: SimpleGraph, edges, induced_subgraph, ne, nv, vertices

using MetaGraphsNext: MetaGraph, edge_labels, labels

using ProgressMeter: @showprogress

using Random: seed!

using Nucleus

using ImmuneReceptor

# ---- #

# loading data
df = load_cdr3s(/path/to/your/data.csv)
hla_df = CSV.read(/path/to/your/hla_data.csv, DataFrame) # if choosing to analyse hla also

# find motifs
all_motifs = get_motifs(list_of_cdrs, min_length, max_length)
sig_motifs = find_significant_motifs(all_motifs, sample_cdrs, reference_cdrs)

# make graph
graph = make_edges(df, sig_motifs, isglobal = true, islocal = true)

# score
length_pvals = find_length_pvals(graph, reference_cdrs, sim_depth)
vgene_pvals = score_vgene(graph, sim_depth)
hla_pvals = score_hla(graph, hla_df) # see above for hla_df
