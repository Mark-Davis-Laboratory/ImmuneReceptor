#0 load functions
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "ImmuneReceptor.jl"))

# add postional ID to original df
# in graph when tracking and when u make cdrs unique, use this positional ID, call it a barcode

const PK = pkgdir(ImmuneReceptor)

const IN = joinpath(PK, "in")

const OU = joinpath(PK, "ou")

PK = @__DIR__
const IN = joinpath(PK, "in")
const OU = joinpath(PK, "ou")

# packages
using .ImmuneReceptor
using CSV
using DataFrames
using ProgressMeter
using Statistics
using HypothesisTests
using StatsBase
using Distributions
using Random
using Graphs
using MetaGraphsNext
using CairoMakie
using GraphMakie
using NetworkLayout
using FASTX

#1. load data
df = CSV.read("/Users/dr/Downloads/Dagny_Reese/a1bf881a/RNA_TCR/b.clones_TRB.tsv", DataFrame; delim='\t')
cdrs = load_cdr3s("/Users/dr/Downloads/Dagny_Reese/a1bf881a/RNA_TCR/b.clones_TRB.tsv")
fasta_path = joinpath(IN, "tcrab-naive-refdb.fa")
cdrs_ref = load_cdr3_fasta(fasta_path)

cdrs = cdrs[.!occursin.(r"[*_]", cdrs.cdr3), :]
println("after cleaning: $(nrow(cdrs)) CDR3s")

const VGENE_FREQ_PATH = joinpath(IN, "tcrb-human.v-freq.txt")
const LEN_FREQ_PATH = joinpath(IN, "tcrb-human.cdr3len-freq.txt")

##################################################

# helpful functions specific to current dataset:
#
function fix_tsv(path::AbstractString)

    df = CSV.read(path, DataFrame; delim='\t')

    # rename columns to standard names
    rename!(df,
        :allVHitsWithScore => :vgene,
        :allDHitsWithScore => :dgene,
        :allJHitsWithScore => :jgene,
        :allCHitsWithScore => :cgene,
        :aaSeqCDR3 => :cdr3,
        :nSeqCDR3 => :nucleotide_cdr3,
        :readCount => :reads,
        :uniqueMoleculeCount => :uniqueMoleculeCount,
    )

    # normalise gene names — take only the part before the '*'
    for col in (:vgene, :dgene, :jgene, :cgene)
        col in propertynames(df) || continue
        df[!, col] = map(df[:, col]) do val
            ismissing(val) && return val
            s = String(val)
            isempty(s) && return val
            contains(s, '*') ? split(s, '*')[1] : s
        end
    end

    return df

end

function normalise_gene_col!(df::DataFrame, col::Symbol)

    col in propertynames(df) || return df

    df[!, col] = map(df[:, col]) do val

        if ismissing(val) || val === nothing
            return val
        end

        if val isa AbstractVector
            return [contains(String(v), '*') ? split(String(v), '*')[1] : String(v) for v in val]
        else
            s = String(val)
            return contains(s, '*') ? split(s, '*')[1] : s
        end

    end

    return df
end

##################################################

# add row index — required for make_edges
cdrs[!, :row_index] = 1:nrow(cdrs)

println("loaded $(nrow(cdrs)) unique cdr3s")

# make qc plots
fig_vgene = plot_distribution(summarize_vgene(raw_df); title="V Gene Usage")
fig_jgene = plot_distribution(summarize_jgene(raw_df); title="J Gene Usage")
fig_dgene = plot_distribution(summarize_dgene(raw_df); title="D Gene Usage")
fig_len = plot_distribution(summarize_cdr3_lengths(raw_df); title="CDR3 Length Distribution")

save(joinpath(OU, "sample1_trb_vgene.png"), fig_vgene)
save(joinpath(OU, "qc_jgene.png"), fig_jgene)
save(joinpath(OU, "qc_dgene.png"), fig_dgene)
save(joinpath(OU, "qc_lengths.png"), fig_len)

##################################################

# find motifs
all_motifs = get_motifs(
    cdrs.cdr3,
    3,      # min motif length
    5,      # max motif length
    false,  # discontiguous — set true to include gapped motifs
)

println("found $(length(all_motifs)) candidate motifs")
# significant motifs

sig_motifs = find_significant_motifs(
    all_motifs,
    cdrs.cdr3,
    cdrs_ref;
    pvalue_cutoff=0.01,
    min_fold=(1000.0, 100.0, 10.0),
    fold_lengths=(2, 3, 4),
    min_count=3,
    use_fold_cutoff=true,
)

println("found $(length(sig_motifs)) significant motifs")

##################################################

cdrs[!, :row_index] = 1:nrow(cdrs)

# build graph
g = make_edges(
    cdrs,
    sig_motifs,
    true,   # isglobal
    true,   # islocal
)

println("graph: $(nv(g)) vertices, $(ne(g)) edges")

##################################################

# filter graph

min_cluster_size = 10

g_filtered = remove_small_clusters(g, min_cluster_size)

clusters, sizes = get_clusters_and_sizes(g_filtered)
println("After filtering: $(length(clusters)) clusters, $(nv(g_filtered)) vertices")

##################################################

# scoring

cluster_summary = summarize_clusters(
    g_filtered,
    min_cluster_size;
    sim_depth=1000,
    cdrs_ref=cdrs_ref,
    motifs=sig_motifs,
    isglobal=true,
    islocal=true,
    nsim=1000,
    run_length=true,
    run_vgene=true,
    run_jgene=false,
    run_clone=false,
    run_size=false,
)

# cdrs in each cluster w/ row id from orig table
cdr3_summary = summarize_cluster_cdr3s(g_filtered, min_cluster_size)

println("scored $(nrow(cluster_summary)) clusters")

##################################################

#o utputs

CSV.write(joinpath(OU, "cluster_summary.csv"), cluster_summary)
CSV.write(joinpath(OU, "cluster_cdr3s.csv"), cdr3_summary)

println("saved cluster_summary.csv and cluster_cdr3s.csv to $OU")
