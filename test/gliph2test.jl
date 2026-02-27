##################################################
### GLIPH2 WORKFLOW
##################################################

### 1. LOAD DATA

raw_df = CSV.read("path/to/your/input.csv", DataFrame)
cdrs = load_cdr3s("path/to/your/input.csv")
cdrs_ref = load_cdr3_fasta("path/to/your/reference.fasta")  # or just a String vector

# add row index — required for make_edges
cdrs[!, :row_index] = 1:nrow(cdrs)

println("Loaded $(nrow(cdrs)) unique CDR3s")

##################################################

### 2. INPUT QC PLOTS

fig_vgene = plot_distribution(summarize_vgene(raw_df); title="V Gene Usage")
fig_jgene = plot_distribution(summarize_jgene(raw_df); title="J Gene Usage")
fig_dgene = plot_distribution(summarize_dgene(raw_df); title="D Gene Usage")
fig_len = plot_distribution(summarize_cdr3_lengths(raw_df); title="CDR3 Length Distribution")

save(joinpath(OU, "qc_vgene.png"), fig_vgene)
save(joinpath(OU, "qc_jgene.png"), fig_jgene)
save(joinpath(OU, "qc_dgene.png"), fig_dgene)
save(joinpath(OU, "qc_lengths.png"), fig_len)

##################################################

### 3. FIND MOTIFS

# get all motifs from the sample CDR3s
all_motifs = get_motifs(
    cdrs.cdr3,
    2,      # min motif length
    4,      # max motif length
    false,  # discontiguous — set true to include gapped motifs
)

println("Found $(length(all_motifs)) candidate motifs")

##################################################

### 4. FIND SIGNIFICANT MOTIFS vs REFERENCE

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

println("Found $(length(sig_motifs)) significant motifs")

##################################################

### 5. BUILD GRAPH

g = make_edges(
    cdrs,
    sig_motifs,
    true,   # isglobal — hamming distance edges
    true,   # islocal  — shared motif edges
)

println("Graph: $(nv(g)) vertices, $(ne(g)) edges")

##################################################

### 6. FILTER GRAPH — remove small clusters

min_cluster_size = 2

g_filtered = remove_small_clusters(g, min_cluster_size)

clusters, sizes = get_clusters_and_sizes(g_filtered)
println("After filtering: $(length(clusters)) clusters, $(nv(g_filtered)) vertices")

##################################################

### 7. SCORE GRAPH

cluster_summary = summarize_clusters(
    g_filtered,
    min_cluster_size;
    sim_depth=1000,
    cdrs_ref=cdrs,       # used for size pval resampling
    motifs=sig_motifs,
    isglobal=true,
    islocal=true,
    nsim=1000,
    run_length=true,
    run_vgene=true,
    run_jgene=true,
    run_clone=true,
    run_size=true,
)

# also get per-cluster CDR3 membership
cdr3_summary = summarize_cluster_cdr3s(g_filtered, min_cluster_size)

println("Scored $(nrow(cluster_summary)) clusters")

##################################################

### 8. SAVE OUTPUTS

CSV.write(joinpath(OU, "cluster_summary.csv"), cluster_summary)
CSV.write(joinpath(OU, "cluster_cdr3s.csv"), cdr3_summary)

println("Saved cluster_summary.csv and cluster_cdr3s.csv to $OU")

##################################################

### 9. CLUSTER PLOTS

# plot all clusters in a grid
fig_clusters = plot_all_clusters(g_filtered, min_cluster_size; ncols=3)
save(joinpath(OU, "clusters.png"), fig_clusters)

# optionally plot individual clusters of interest — e.g. top 5 by size
top_clusters = first(sort(cluster_summary, :cluster_size; rev=true), 5)

for row in eachrow(top_clusters)
    fig = plot_cluster(g_filtered, row.cluster_id)
    save(joinpath(OU, "cluster_$(row.cluster_id).png"), fig)
end

println("Done — all outputs saved to $OU")
