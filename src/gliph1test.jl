using ImmuneReceptor

using CSV
using DataFrames

using Graphs
using MetaGraphsNext

using StatsBase: sample, mode
using Statistics: mean
using Random
using Distributions

using ProgressMeter: @showprogress
using Nucleus


# ---- #

# to do
# load practice data and test
# load cdr3s from paper
# do analysis on paper cdr3s and compare

# loading data
reference_cdrs = load_cdr3_fasta("in/gliph1/gliph/db/warren-naive.fa")
reference_cdrs = load_cdr3_fasta("in/gliph1/gliph/db/warren-naive.fa")

df = load_cdr3s("in/gliph1/gliph/db/test_data_originalgliph.csv")
#hla_df = CSV.read(/path/to/your/hla_data.csv, DataFrame) # if choosing to analyse hla also

# find motifs
#all_motifs = get_motifs(cd.cdr3, min_length, max_length)
all_motifs = get_motifs(df.cdr3, 3, 5)
#sig_motifs = find_significant_motifs(all_motifs, sample_cdrs, reference_cdrs, nsim, ove_cutoff)
sig_motifs = find_significant_motifs(all_motifs, df.cdr3, reference_cdrs, 1000, true)

# make graph
graph = make_edges(df, sig_motifs, true, true)

# score
length_pvals = find_length_pvals(graph, reference_cdrs, sim_depth)
vgene_pvals = score_vgene(graph, sim_depth)
#hla_pvals = score_hla(graph, hla_df) # see above for hla_df
cluster_motifs = motif_scoring(graph)

summarize_data(graph, vgene_pvals, length_pvals, cluster_motifs)
