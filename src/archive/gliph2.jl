# gliph 2

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
using FASTX

#################
### PATHS
#################

const PK = @__DIR__
const IN = joinpath(PK, "in")
const OU = joinpath(PK, "ou")

const VGENE_FREQ_PATH = joinpath(IN, "tcrb-human_v-freq.txt")
const LEN_FREQ_PATH = joinpath(IN, "tcrb-human_cdr3len-freq.txt")

#################
### READING
#################

function is_t_gene(::Any)

    false

end

function is_t_gene(st::AbstractString)

    startswith(st, 'T')

end

function is_cdr3(st)

    st[1] == 'C' && st[end] == 'F'

end

using CSV, DataFrames

function load_cdr3s(csvpath)

    df = CSV.read(csvpath, DataFrame)

    cdr3 = String.(coalesce.(df.cdr3, ""))
    cdr3 = replace.(cdr3, r"\s+" => "")
    cdr3 = uppercase.(cdr3)

    bad = Set(["", "NA", "NONE", "NAN", "NULL"])
    keep = .!(cdr3 .∈ Ref(bad))
    keep .&= startswith.(cdr3, "C")
    keep .&= endswith.(cdr3, "F")

    df = df[keep, :]
    df.cdr3 = cdr3[keep]
    df.cdr3_length = length.(df.cdr3)

    new = groupby(df, :cdr3)

    new_cols = Any[nrow=>:duplicate_count, :cdr3_length=>first=>:cdr3_length]

    if "barcode" in names(df)
        push!(new_cols, :barcode => (x -> collect(unique(x))) => :barcodes)
    end

    if "sample" in names(df)
        push!(new_cols, :sample => (x -> collect(unique(x))) => :samples)
    end

    if "vgene" in names(df)
        push!(new_cols, :vgene => (x -> collect(unique(x))) => :vgenes)
    end
    if "jgene" in names(df)
        push!(new_cols, :jgene => (x -> collect(unique(x))) => :jgenes)
    end
    if "dgene" in names(df)
        push!(new_cols, :dgene => (x -> collect(unique(x))) => :dgenes)
    end

    if "clones" in names(df)
        df.clones = Float64.(coalesce.(df.clones, 0))
        push!(new_cols, :clones => (x -> collect(x)) => :clones)
        push!(new_cols, :clones => sum => :clones_sum)
    end

    if "umis" in names(df)
        df.umis = Float64.(coalesce.(df.umis, 0))
        push!(new_cols, :umis => (x -> collect(x)) => :clones)
        push!(new_cols, :umis => sum => :clones_sum)
    end

    out = combine(new, new_cols...)

    return out
end

function load_cdr3_fasta(path::AbstractString)
    cdr3s = String[]
    open(FASTA.Reader, path) do reader
        for record in reader
            push!(cdr3s, String(sequence(record)))
        end
    end
    return cdr3s
end

#################
### CDR3 and Scoring
#################

using ProgressMeter

function make_hamming_distance(s1, s2)

    sum(s1[nd] != s2[nd] for nd in eachindex(s1))

end

const blosum62 = Dict{Tuple{Char,Char},Int}(
    ('A', 'A') => 4, ('A', 'R') => -1, ('A', 'N') => -2, ('A', 'D') => -2, ('A', 'C') => 0,
    ('A', 'Q') => -1, ('A', 'E') => -1, ('A', 'G') => 0, ('A', 'H') => -2, ('A', 'I') => -1,
    ('A', 'L') => -1, ('A', 'K') => -1, ('A', 'M') => -1, ('A', 'F') => -2, ('A', 'P') => -1,
    ('A', 'S') => 1, ('A', 'T') => 0, ('A', 'W') => -3, ('A', 'Y') => -2, ('A', 'V') => 0, ('R', 'R') => 5, ('R', 'N') => 0, ('R', 'D') => -2, ('R', 'C') => -3, ('R', 'Q') => 1,
    ('R', 'E') => 0, ('R', 'G') => -2, ('R', 'H') => 0, ('R', 'I') => -3, ('R', 'L') => -2,
    ('R', 'K') => 2, ('R', 'M') => -1, ('R', 'F') => -3, ('R', 'P') => -2, ('R', 'S') => -1,
    ('R', 'T') => -1, ('R', 'W') => -3, ('R', 'Y') => -2, ('R', 'V') => -3, ('N', 'N') => 6, ('N', 'D') => 1, ('N', 'C') => -3, ('N', 'Q') => 0, ('N', 'E') => 0,
    ('N', 'G') => 0, ('N', 'H') => 1, ('N', 'I') => -3, ('N', 'L') => -3, ('N', 'K') => 0,
    ('N', 'M') => -2, ('N', 'F') => -3, ('N', 'P') => -2, ('N', 'S') => 1, ('N', 'T') => 0,
    ('N', 'W') => -4, ('N', 'Y') => -2, ('N', 'V') => -3, ('D', 'D') => 6, ('D', 'C') => -3, ('D', 'Q') => 0, ('D', 'E') => 2, ('D', 'G') => -1,
    ('D', 'H') => -1, ('D', 'I') => -3, ('D', 'L') => -4, ('D', 'K') => -1, ('D', 'M') => -3,
    ('D', 'F') => -3, ('D', 'P') => -1, ('D', 'S') => 0, ('D', 'T') => -1, ('D', 'W') => -4,
    ('D', 'Y') => -3, ('D', 'V') => -3, ('C', 'C') => 9, ('C', 'Q') => -3, ('C', 'E') => -4, ('C', 'G') => -3, ('C', 'H') => -3,
    ('C', 'I') => -1, ('C', 'L') => -1, ('C', 'K') => -3, ('C', 'M') => -1, ('C', 'F') => -2,
    ('C', 'P') => -3, ('C', 'S') => -1, ('C', 'T') => -1, ('C', 'W') => -2, ('C', 'Y') => -2,
    ('C', 'V') => -1, ('Q', 'Q') => 5, ('Q', 'E') => 2, ('Q', 'G') => -2, ('Q', 'H') => 0, ('Q', 'I') => -3,
    ('Q', 'L') => -2, ('Q', 'K') => 1, ('Q', 'M') => 0, ('Q', 'F') => -3, ('Q', 'P') => -1,
    ('Q', 'S') => 0, ('Q', 'T') => -1, ('Q', 'W') => -2, ('Q', 'Y') => -1, ('Q', 'V') => -2, ('E', 'E') => 5, ('E', 'G') => -2, ('E', 'H') => 0, ('E', 'I') => -3, ('E', 'L') => -3,
    ('E', 'K') => 1, ('E', 'M') => -2, ('E', 'F') => -3, ('E', 'P') => -1, ('E', 'S') => 0,
    ('E', 'T') => -1, ('E', 'W') => -3, ('E', 'Y') => -2, ('E', 'V') => -2, ('G', 'G') => 6, ('G', 'H') => -2, ('G', 'I') => -4, ('G', 'L') => -4, ('G', 'K') => -2,
    ('G', 'M') => -3, ('G', 'F') => -3, ('G', 'P') => -2, ('G', 'S') => 0, ('G', 'T') => -2,
    ('G', 'W') => -2, ('G', 'Y') => -3, ('G', 'V') => -3, ('H', 'H') => 8, ('H', 'I') => -3, ('H', 'L') => -3, ('H', 'K') => -1, ('H', 'M') => -2,
    ('H', 'F') => -1, ('H', 'P') => -2, ('H', 'S') => -1, ('H', 'T') => -2, ('H', 'W') => -2,
    ('H', 'Y') => 2, ('H', 'V') => -3, ('I', 'I') => 4, ('I', 'L') => 2, ('I', 'K') => -3, ('I', 'M') => 1, ('I', 'F') => 0,
    ('I', 'P') => -3, ('I', 'S') => -2, ('I', 'T') => -1, ('I', 'W') => -3, ('I', 'Y') => -1,
    ('I', 'V') => 3, ('L', 'L') => 4, ('L', 'K') => -2, ('L', 'M') => 2, ('L', 'F') => 0, ('L', 'P') => -3,
    ('L', 'S') => -2, ('L', 'T') => -1, ('L', 'W') => -2, ('L', 'Y') => -1, ('L', 'V') => 1, ('K', 'K') => 5, ('K', 'M') => -1, ('K', 'F') => -3, ('K', 'P') => -1, ('K', 'S') => 0,
    ('K', 'T') => -1, ('K', 'W') => -3, ('K', 'Y') => -2, ('K', 'V') => -2, ('M', 'M') => 5, ('M', 'F') => 0, ('M', 'P') => -2, ('M', 'S') => -1, ('M', 'T') => -1,
    ('M', 'W') => -1, ('M', 'Y') => -1, ('M', 'V') => 1, ('F', 'F') => 6, ('F', 'P') => -4, ('F', 'S') => -2, ('F', 'T') => -2, ('F', 'W') => 1,
    ('F', 'Y') => 3, ('F', 'V') => -1, ('P', 'P') => 7, ('P', 'S') => -1, ('P', 'T') => -1, ('P', 'W') => -4, ('P', 'Y') => -3,
    ('P', 'V') => -2, ('S', 'S') => 4, ('S', 'T') => 1, ('S', 'W') => -3, ('S', 'Y') => -2, ('S', 'V') => -2, ('T', 'T') => 5, ('T', 'W') => -2, ('T', 'Y') => -2, ('T', 'V') => 0, ('W', 'W') => 11, ('W', 'Y') => 2, ('W', 'V') => -3, ('Y', 'Y') => 7, ('Y', 'V') => -1, ('V', 'V') => 4
)

function make_symmetric!(d::Dict{Tuple{Char,Char},Int})
    for ((a, b), score) in collect(d)
        d[(b, a)] = score
    end
    return d
end

const blosum62_sym = make_symmetric!(blosum62)

function make_blosum_score(s1, s2)

    @assert length(s1) == length(s2)

    mismatch_indices = Int[]

    for nd in eachindex(s1)
        if s1[nd] != s2[nd]
            push!(mismatch_indices, nd)
        end
    end

    if isempty(mismatch_indices)
        return 0
    elseif length(mismatch_indices) > 1
        return -Inf
    end

    nd = mismatch_indices[1]

    a1 = s1[nd]
    a2 = s2[nd]

    return get(blosum62_sym, (a1, a2), -Inf)

end

function make_distance(st_, ids)

    u1 = lastindex(st_)

    u2 = div(u1 * (u1 - 1), 2)

    in__ = Vector{Tuple{Int,Int}}(undef, u2)

    po_ = Vector{Int}(undef, u2)

    i1 = 0

    @showprogress for i2 in 1:u1, i3 in (i2+1):u1

        s1 = st_[i2]
        s2 = st_[i3]

        if lastindex(s1) != lastindex(s2)

            continue

        end

        in__[i1+=1] = (Int(ids[i2]), Int(ids[i3]))
        po_[i1] = make_hamming_distance(s1, s2)

    end

    resize!(in__, i1)
    resize!(po_, i1)

    return in__, po_
end

#################
### Motifs
#################

# count whether motif is present at least once per cdr3 for all cdr3s

function get_motif_counts(motif::AbstractVector{<:AbstractString}, cdrs::AbstractVector{<:AbstractString},)

    motif_counts = Dict{String,Int}(m => 0 for m in motif)

    cdrs = unique(cdrs)
    cores = [c[4:(end-3)] for c in cdrs if lastindex(c) >= 7]

    for cdr in cores

        for m in motif

            if occursin('.', m)

                hit = occursin(Regex(m), cdr)

            else

                hit = occursin(m, cdr)

            end

            if hit
                motif_counts[m] += 1
            end

        end

    end

    return motif_counts

end

# makes list of motifs in cdr3s, counts them, filters based on a cutoff, then provides set count
function get_motifs(st_::AbstractVector{<:AbstractString}, min::Int, max::Int, discontiguous=false)

    mo_ = Dict{String,Int}() # make dictionary

    st_ = unique(st_)
    cores = [c[4:(end-3)] for c in st_ if lastindex(c) >= 7]
    #cores = unique(cores)

    @showprogress desc = "Finding and counting motifs..." for s1 in cores

        lastindex(s1) < 7 && continue # only keep cdr3 7 or longer

        #s2 = s1[4:(end-3)] # remove first and last 3 aa
        s2 = s1
        lind = lastindex(s2)

        seen = Set{String}()

        for um in min:max

            um > lastindex(s2) && continue # make sure cdr3 is longer than the motif size

            i1 = 0
            i2 = i1 + um - 1

            for i1 in 1:(lind-um+1)

                i2 = i1 + um - 1

                m = s2[i1:i2]

                push!(seen, m)

            end

        end

        if discontiguous == true
            for span in (4, 5)

                lastindex(s2) < span && continue

                for i in 1:(lastindex(s2)-span+1)

                    w = s2[i:i+span-1]

                    for g in 2:(span-1)
                        p = string(w[1:g-1], ".", w[g+1:end])
                        push!(seen, p)

                    end
                end
            end
        end

        for m in seen
            mo_[m] = get!(mo_, m, 0) + 1
        end

    end

    mo2_ = Dict(m => num for (m, num) in mo_ if num >= 3)
    return mo2_

end

function find_significant_motifs(motifs, cdrs1, cdrs2; pvalue_cutoff::Float64=0.01, min_fold=(1000.0, 100.0, 10.0), fold_lengths=(2, 3, 4), min_count::Int=3, use_fold_cutoff::Bool=true, pseudocount::Float64=0.01,)

    motifs_list = collect(keys(motifs))
    L = length(motifs_list)

    pat = Vector{Union{Nothing,Regex}}(undef, L)
    for j in 1:L
        m = motifs_list[j]
        pat[j] = occursin('.', m) ? Regex(m) : nothing
    end

    cdrs1 = unique(cdrs1)
    cdrs2 = unique(cdrs2)

    cores1 = unique([c[4:(end-3)] for c in cdrs1 if lastindex(c) >= 7])
    cores2 = unique([c[4:(end-3)] for c in cdrs2 if lastindex(c) >= 7])

    n1 = length(cores1)  # sample unique
    n2 = length(cores2)  # reference unique

    # no data => no motifs
    if n1 == 0 || n2 == 0 || L == 0
        return Dict{String,Float64}()
    end

    # per-motif fold threshold (GLIPH2: depends on motif length; discontinuous motifs count as length-1)
    minfold_for = Vector{Float64}(undef, L)
    if length(min_fold) == 1
        fill!(minfold_for, float(min_fold))
    else

        @assert length(min_fold) == length(fold_lengths)

        for j in 1:L
            m = motifs_list[j]
            eff_len = length(m) - (occursin('.', m) ? 1 : 0)

            thr = 0.0
            for (k, Lk) in pairs(fold_lengths)
                if eff_len == Lk
                    thr = float(min_fold[k])
                    break
                end
            end

            # if it doesn't match any provided fold_lengths, set to 0 (i.e., don't filter by fold)
            minfold_for[j] = thr

        end
    end

    significant = Dict{String,Float64}()

    for j in 1:L
        m = motifs_list[j]

        a = 0
        for c in cores1
            hit = pat[j] === nothing ? occursin(m, c) : occursin(pat[j], c)
            a += hit ? 1 : 0
        end

        c = 0
        for cdr in cores2
            hit = pat[j] === nothing ? occursin(m, cdr) : occursin(pat[j], cdr)
            c += hit ? 1 : 0
        end

        a < min_count && continue

        b = n1 - a
        d = n2 - c

        ove = (a / n1) / ((c + pseudocount) / n2)

        if use_fold_cutoff
            ove >= minfold_for[j] || continue
        end

        p = pvalue(FisherExactTest(a, b, c, d); tail=:right)
        p <= pvalue_cutoff || continue

        significant[m] = p  # raw p-value, like GLIPH2’s “motif enrichment” filtering step

    end

    return significant

end

#_________#

# takes list of cdrs and checks for a motif, making pairs

function make_motif_pairs(st_::AbstractVector{<:AbstractString}, motif::AbstractString)
    u1 = lastindex(st_)
    u2 = div(u1 * (u1 - 1), 2)

    in__ = Vector{Tuple{Int,Int}}(undef, u2)
    po_ = Vector{Int}(undef, u2)

    i1 = 0

    booler = nothing

    if occursin('.', motif)
        booler = Regex(motif)
    end

    for i2 in 1:u1, i3 in (i2+1):u1

        s1, s2 = st_[i2], st_[i3]
        core1 = lastindex(s1) < 7 ? "" : s1[4:(end-3)]
        core2 = lastindex(s2) < 7 ? "" : s2[4:(end-3)]

        has1 = false
        has2 = false

        if booler === nothing
            has1 = occursin(motif, core1)
            has2 = occursin(motif, core2)
        else
            has1 = occursin(booler, core1)
            has2 = occursin(booler, core2)
        end

        if has1 && has2

            nd1 = nothing
            nd2 = nothing

            if booler === nothing
                nd1 = findfirst(motif, core1)
                nd2 = findfirst(motif, core2)
            else
                nd1 = findfirst(booler, core1)
                nd2 = findfirst(booler, core2)
            end

            a1 = first(nd1)
            b1 = first(nd2)

            if abs(a1 - b1) ≤ 3
                i1 += 1
                in__[i1] = (i2, i3)
                po_[i1] = 1
            end

        end

    end

    resize!(in__, i1)
    resize!(po_, i1)

    return in__, po_
end

#################
### Making Edges
#################

function make_vertices!(g, cdrs)
    isblank(x) = x === nothing || x === missing || x in ("None", "none", "NA", "na", "Na", " ", "  ", "")

    for row in 1:nrow(cdrs)
        label = String(cdrs.cdr3[row])

        add_vertex!(g, label, Dict(:index => row))

        g[label][:cdr3_length] = cdrs.cdr3_length[row]
        g[label][:duplicate_count] = cdrs.duplicate_count[row]

        if "vgenes" in names(cdrs) && cdrs.vgenes[row] !== missing && cdrs.vgenes[row] !== nothing
            v = cdrs.vgenes[row]
            if !(v isa AbstractVector && isempty(v))
                g[label][:vgene] = v
            end
        end

        if "jgenes" in names(cdrs) && cdrs.jgenes[row] !== missing && cdrs.jgenes[row] !== nothing
            v = cdrs.jgenes[row]
            if !(v isa AbstractVector && isempty(v))
                g[label][:jgene] = v
            end
        end

        if "dgenes" in names(cdrs) && cdrs.dgenes[row] !== missing && cdrs.dgenes[row] !== nothing
            v = cdrs.dgenes[row]
            if !(v isa AbstractVector && isempty(v))
                g[label][:dgene] = v
            end
        end

        if "barcodes" in names(cdrs) && cdrs.barcodes[row] !== missing && cdrs.barcodes[row] !== nothing
            v = cdrs.barcodes[row]
            if !(v isa AbstractVector && isempty(v))
                g[label][:barcodes] = v
            end
        end

        if "samples" in names(cdrs) && cdrs.samples[row] !== missing && cdrs.samples[row] !== nothing
            v = cdrs.samples[row]
            if !(v isa AbstractVector && isempty(v))
                g[label][:samples] = v
            end
        end

        if "clones" in names(cdrs) && cdrs.clones[row] !== missing && cdrs.clones[row] !== nothing
            g[label][:clones] = cdrs.clones[row]
        end

    end

    return g
end

# make a global edge w/ annotation of hamming disrance
function add_global_edge!(g, u, v, distance)
    add_edge!(g, u, v, Dict(:distance => distance))
end

# add a local edge w/ annotation of motif and the pval (from get_significant_motifs)
function add_local_edge!(g, u, v, motif::String, pval)
    add_edge!(g, u, v, Dict(
        :motifs => [motif],
        :motif_pvals => [pval],
    ))
end

# initialises the graph and makes edges
function make_edges(cdrs, motifs, isglobal, islocal)

    cdrs.cdr3 = String.(cdrs.cdr3) # ensure is string

    g = MetaGraph(
        Graph();
        label_type=Int,
        vertex_data_type=Dict{Symbol,Any},
        edge_data_type=Dict{Symbol,Any},
    )

    g = make_vertices!(g, cdrs)

    if isglobal == true

        pairs, dists = make_distance(cdrs.cdr3, cdrs.row_index)

        @showprogress desc = "Making edges..." for (index, (u, v)) in enumerate(pairs)

            d = dists[index]

            # check distance and calculate blosum score if needed
            if d == 1

                s1 = cdrs.cdr3[u]
                s2 = cdrs.cdr3[v]

                make_blosum_score(s1, s2) >= 0 || continue

            else
                d <= 1 || continue
            end

            # add edge
            if haskey(g, u, v)
                g[u, v][:distance] = d
            else
                add_global_edge!(g, u, v, d)
            end

        end
    end


    if islocal == true

        @showprogress desc = "Making edges..." for (motif, pval) in motifs

            pairs, _ = make_motif_pairs(cdrs.cdr3, motif)

            for (u, v) in pairs

                if haskey(g, u, v)

                    if !haskey(g[u, v], :motifs)
                        g[u, v][:motifs] = String[]

                    end

                    if !haskey(g[u, v], :motif_pvals)
                        g[u, v][:motif_pvals] = Float64[]
                    end

                    push!(g[u, v][:motifs], motif)
                    push!(g[u, v][:motif_pvals], pval)

                else
                    add_local_edge!(g, u, v, motif, pval)
                end


            end
        end
    end # end of local edges

    return g

end

function get_clusters_and_sizes(g)
    clusters = connected_components(g)
    sizes = [length(c) for c in clusters]
    return clusters, sizes
end

function remove_small_clusters(g, min_size::Int)

    # get clusters and filter to only those meeting the minimum size
    clusters = connected_components(g)
    keep_labels = Set{Int}()

    for cluster in clusters
        if length(cluster) >= min_size
            for v in cluster
                push!(keep_labels, label_for(g, v))
            end
        end
    end

    # initialise new graph of the same type
    g_new = MetaGraph(
        Graph();
        label_type=Int,
        vertex_data_type=Dict{Symbol,Any},
        edge_data_type=Dict{Symbol,Any},
    )

    # copy kept vertices with their metadata
    for lbl in keep_labels
        add_vertex!(g_new, lbl, deepcopy(g[lbl]))
    end

    # copy edges where both endpoints are kept
    for lbl in keep_labels
        v = code_for(g, lbl)
        for u in neighbors(g, v)
            lu = label_for(g, u)
            if lu in keep_labels && lbl < lu
                add_edge!(g_new, lbl, lu, deepcopy(g[lbl, lu]))
            end
        end
    end

    n_removed = nv(g) - nv(g_new)
    println("Removed $n_removed vertices across $(length(clusters) - count(c -> length(c) >= min_size, clusters)) clusters.")

    return g_new
end

#################
### Scoring / P-Vals
#################

# length score
# for each group, resample 1000 groups of group size n and measure cdr3 length distribution.

function load_cdr3len_freqs(path::AbstractString)

    freqs = Dict{Int,Float64}()

    open(path) do f

        for line in eachline(f)
            isempty(strip(line)) && continue
            parts = split(strip(line), '\t')
            length(parts) == 2 || continue
            len = parse(Int, replace(String(parts[1]), "Len" => ""))
            freq = parse(Float64, parts[2])
            freqs[len] = freq
        end

    end

    return freqs
end

function find_length_pvals(g, sim_depth)

    len_freq_path = joinpath(IN, "tcrb-human_cdr3len-freq.txt")
    ref_freqs = load_cdr3len_freqs(len_freq_path)

    len_keys = collect(keys(ref_freqs))
    len_probs = Weights([ref_freqs[l] for l in len_keys])

    clusters = connected_components(g)

    sample_scores = Float64[]
    ns = Int[]

    for cluster in clusters

        labels = label_for.(Ref(g), cluster)

        cluster_lengths = Int[]

        for lbl in labels
            push!(cluster_lengths, g[lbl][:cdr3_length])
        end

        n = length(cluster_lengths)

        if n == 0
            push!(sample_scores, 1.0)
            push!(ns, 0)
            continue
        end

        cl_tab = countmap(cluster_lengths)
        p_len = [cnt / n for cnt in values(cl_tab)]
        sample_score = prod(p_len)

        push!(sample_scores, sample_score)
        push!(ns, n)

    end

    p_vals = Float64[]

    for (index, n) in enumerate(ns)

        if n == 0
            push!(p_vals, 1.0)
            continue
        end

        sim_scores = Float64[]

        for i in 1:sim_depth
            sim_lengths = sample(len_keys, len_probs, n; replace=true)
            sim_tab = countmap(sim_lengths)
            p_len = [cnt / n for cnt in values(sim_tab)]
            push!(sim_scores, prod(p_len))
        end

        p_val = (count(>=(sample_scores[index]), sim_scores) + 1) / (sim_depth + 1)
        push!(p_vals, p_val)

    end

    return p_vals
end

function find_length_pvals_sampling(g, cdrs2, sim_depth)

    clusters = connected_components(g)

    sample_scores = Float64[]
    ns = Int[]

    # get product-of-frequencies score for each cluster
    for cluster in clusters
        labels = label_for.(Ref(g), cluster)

        cluster_lengths = Int[]
        for lbl in labels
            push!(cluster_lengths, g[lbl][:cdr3_length])
        end

        n = length(cluster_lengths)
        if n == 0 # skip if empty
            push!(sample_scores, 1.0)
            push!(ns, 0)
            continue
        end

        cl_tab = countmap(cluster_lengths)
        p_len = [cnt / n for cnt in values(cl_tab)]

        sample_score = prod(p_len)

        push!(sample_scores, sample_score)
        push!(ns, n)
    end

    p_vals = Float64[]
    for (index, n) in enumerate(ns)
        if n == 0
            push!(p_vals, 1.0)
            continue
        end

        sim_scores = Float64[]
        for i in 1:sim_depth
            random_cdrs = sample(cdrs2, n; replace=true, ordered=false)
            sim_lengths = [length(s) for s in random_cdrs]

            sim_tab = countmap(sim_lengths)
            p_len = [cnt / n for cnt in values(sim_tab)]
            push!(sim_scores, prod(p_len))
        end

        p_val = (count(>=(sample_scores[index]), sim_scores) + 1) / (sim_depth + 1)
        push!(p_vals, p_val)

    end

    return p_vals

end

# network size pval

function find_size_pval(g, cdrs_ref, motifs; nsim::Int=1000, isglobal::Bool=true, islocal::Bool=true, rng::AbstractRNG=Random.default_rng())

    n_sample = nv(g)
    n_ref = nrow(cdrs_ref)

    clusters, sizes_obs = get_clusters_and_sizes(g)

    sim_all_sizes = Int[]

    for x in 1:nsim
        idx = sample(rng, 1:n_ref, n_sample; replace=true)
        sub = cdrs_ref[idx, :]
        sub[!, :row_index] = 1:nrow(sub)

        g_ref = make_edges(sub, motifs, isglobal, islocal)
        _, sim_sizes = get_clusters_and_sizes(g_ref)

        append!(sim_all_sizes, sim_sizes)
    end

    pvals = Float64[]
    for s_obs in sizes_obs
        wins = count(x -> x >= s_obs, sim_all_sizes)
        p = wins / length(sim_all_sizes)
        push!(pvals, p == 0.0 ? 1.0 / nsim : p)
    end

    return pvals
end

# v gene and j gene enrichment score

using Distributions

using StatsBase, Statistics

using StatsBase

function load_vgene_freqs(path::AbstractString)

    freqs = Dict{String,Float64}()
    open(path) do f

        for line in eachline(f)
            isempty(strip(line)) && continue
            parts = split(strip(line), '\t')
            length(parts) == 2 || continue
            gene = String(parts[1])
            freq = parse(Float64, parts[2])
            freqs[gene] = freq
        end

    end

    return freqs
end

function find_vgene_pval(g, sim_depth)

    vgene_freq_path = joinpath(IN, "tcrb-human_v-freq.txt")
    ref_freqs = load_vgene_freqs(vgene_freq_path)

    v_keys = collect(keys(ref_freqs))
    v_probs = Weights([ref_freqs[v] for v in v_keys])

    clusters = connected_components(g)

    counts = Vector{Dict{String,Int}}(undef, length(clusters))
    entry_sizes = Vector{Int}(undef, length(clusters))

    for (index, cluster) in enumerate(clusters)

        di = Dict{String,Int}()
        entry_count = 0

        for vertex in cluster

            lbl = label_for(g, vertex)
            haskey(g[lbl], :vgene) || continue

            vgenes = g[lbl][:vgene]
            vgenes = vgenes isa AbstractVector ? vgenes : [vgenes]

            for vg in vgenes
                s = String(vg)
                di[s] = get!(di, s, 0) + 1
                entry_count += 1
            end

        end

        counts[index] = di
        entry_sizes[index] = entry_count

    end

    cluster_pvals = Vector{Dict{String,Float64}}(undef, length(clusters))

    for (index, di) in enumerate(counts)

        n = entry_sizes[index]

        if isempty(di) || n == 0
            cluster_pvals[index] = Dict("NA" => 1.0)
            continue
        end

        sample_score = sum(log(cnt / n) for cnt in values(di))

        wins = 0
        for i in 1:sim_depth
            random_vgenes = sample(v_keys, Weights(v_probs), n; replace=true)
            sim_tab = countmap(random_vgenes)
            sim_score = sum(log(cnt / n) for cnt in values(sim_tab))

            if sim_score >= sample_score
                wins += 1
            end
        end

        p = wins / sim_depth
        p = (p == 0.0) ? (1.0 / sim_depth) : p

        cluster_pvals[index] = Dict("vgene_spectratype_p" => p)

    end

    return cluster_pvals
end

function find_vgene_pval_sampling(g, sim_depth)

    clusters = connected_components(g)

    # collect all vgenes from the graph for resampling
    totals = Dict{String,Int}()

    for v in vertices(g)
        lbl = label_for(g, v)
        haskey(g[lbl], :vgene) || continue

        vgenes = g[lbl][:vgene]
        vgenes = vgenes isa AbstractVector ? vgenes : [vgenes]

        for vg in vgenes
            s = String(vg)
            totals[s] = get!(totals, s, 0) + 1
        end
    end

    v_keys = collect(keys(totals))
    v_probs = Weights([totals[v] for v in v_keys])
    v_probs = v_probs ./ sum(v_probs)

    # collect per cluster vgene counts
    counts = Vector{Dict{String,Int}}(undef, length(clusters))
    entry_sizes = Vector{Int}(undef, length(clusters))

    for (index, cluster) in enumerate(clusters)
        di = Dict{String,Int}()
        entry_count = 0

        for vertex in cluster
            lbl = label_for(g, vertex)
            haskey(g[lbl], :vgene) || continue

            vgenes = g[lbl][:vgene]
            vgenes = vgenes isa AbstractVector ? vgenes : [vgenes]

            for vg in vgenes
                s = String(vg)
                di[s] = get!(di, s, 0) + 1
                entry_count += 1
            end
        end

        counts[index] = di
        entry_sizes[index] = entry_count
    end

    cluster_pvals = Vector{Dict{String,Float64}}(undef, length(clusters))

    for (index, di) in enumerate(counts)
        n = entry_sizes[index]

        if isempty(di) || n == 0
            cluster_pvals[index] = Dict("NA" => 1.0)
            continue
        end

        sample_score = sum(log(cnt / n) for cnt in values(di))

        wins = 0
        for i in 1:sim_depth
            random_vgenes = sample(v_keys, Weights(v_probs), n; replace=true)
            sim_tab = countmap(random_vgenes)
            sim_score = sum(log(cnt / n) for cnt in values(sim_tab))

            if sim_score >= sample_score
                wins += 1
            end
        end

        p = wins / sim_depth
        p = (p == 0.0) ? (1.0 / sim_depth) : p

        cluster_pvals[index] = Dict("vgene_spectratype_p" => p)
    end

    return cluster_pvals
end

function find_jgene_pval_sampling(g, sim_depth)
    clusters = connected_components(g)

    counts = Vector{Dict{String,Int}}(undef, length(clusters))
    totals = Dict{String,Int}()

    entry_sizes = Vector{Int}(undef, length(clusters))

    for (index, cluster) in enumerate(clusters)
        di = Dict{String,Int}()
        entry_count = 0

        for vertex in cluster
            lbl = label_for(g, vertex)
            haskey(g[lbl], :jgene) || continue

            jgenes = g[lbl][:jgene]

            for jg in jgenes
                j = String(jg)
                di[j] = get!(di, j, 0) + 1
                totals[j] = get!(totals, j, 0) + 1
                entry_count += 1
            end
        end

        counts[index] = di
        entry_sizes[index] = entry_count
    end

    j_keys = collect(keys(totals))
    j_probs = [totals[j] for j in j_keys]
    j_probs = j_probs ./ sum(j_probs)

    cluster_pvals = Vector{Dict{String,Float64}}(undef, length(clusters))

    for (index, di) in enumerate(counts)

        n = entry_sizes[index]

        if isempty(di) || n == 0
            cluster_pvals[index] = Dict("NA" => 1.0)
            continue
        end

        p_j = [cnt / n for cnt in values(di)]
        sample_score = prod(p_j)

        wins = 0
        for i in 1:sim_depth
            random_jgenes = sample(j_keys, Weights(j_probs), n; replace=true)
            sim_tab = countmap(random_jgenes)
            p_sim = [cnt / n for cnt in values(sim_tab)]
            sim_score = prod(p_sim)

            if sim_score >= sample_score
                wins += 1
            end
        end

        p = wins / sim_depth
        p = (p == 0.0) ? (1.0 / sim_depth) : p

        cluster_pvals[index] = Dict("jgene_spectratype_p" => p)
    end

    # write pval back to vertices
    #for (index, cluster) in enumerate(clusters)
    #for vertex in cluster
    #lbl = label_for(g, vertex)
    #g[lbl][:jgene_pval] = cluster_pvals[index]
    #end
    #end

    return cluster_pvals
end

##################### motif scoring
#
function motif_hits_by_cdr3(cdr3s::AbstractVector{<:AbstractString}, motifs::AbstractVector{<:AbstractString})

    hits = Dict{String,Dict{Int,Vector{Int}}}()


    matchers = Dict{String,Union{Nothing,Regex}}()

    for m in motifs

        hits[m] = Dict{Int,Vector{Int}}()
        matchers[m] = occursin('.', m) ? Regex(m) : nothing

    end

    for (i, s_) in pairs(cdr3s)
        s = uppercase(replace(String(s_), r"\s+" => ""))
        core = (lastindex(s) < 7) ? "" : s[4:(end-3)]
        isempty(core) && continue


        for motif in motifs

            rx = matchers[motif]
            start = firstindex(core)
            pos = Int[]

            searcher = rx === nothing ? motif : rx
            start = firstindex(core)

            while start <= lastindex(core)

                match = findnext(searcher, core, start)
                match === nothing && break

                pos_start = first(match)
                push!(pos, pos_start)

                start = pos_start + 1
            end

            if !isempty(pos)
                hits[motif][i] = pos
            end

        end
    end

    return hits

end

function find_cluster_motif_pval(g; ref_hits::Dict{String,<:Any}, N_ref::Int, tail::Symbol=:right)

    clusters = connected_components(g)
    N_sample = nv(g)

    cluster_motifs = Vector{Vector{Tuple{String,Float64}}}(undef, length(clusters))

    for (index, cluster) in enumerate(clusters)

        in_cluster = falses(N_sample)
        in_cluster[cluster] .= true

        motif_vertices = Dict{String,Set{Int}}()

        for u in cluster
            for v in neighbors(g, u)

                if in_cluster[v] && u < v

                    lu = label_for(g, u)
                    lv = label_for(g, v)

                    if haskey(g, lu, lv) && haskey(g[lu, lv], :motifs)

                        motifs = g[lu, lv][:motifs]::Vector{String}

                        for motif in motifs
                            s = get!(motif_vertices, motif, Set{Int}())
                            push!(s, u)
                            push!(s, v)
                        end

                    end
                end
            end
        end

        pairs = Tuple{String,Float64}[]

        for (motif, vertex_set) in motif_vertices

            k_obs = length(vertex_set)
            k_ref = haskey(ref_hits, motif) ? length(ref_hits[motif]) : 0

            if N_ref <= 0 || k_ref <= 0
                push!(pairs, (motif, 1.0))
                continue
            end

            t = FisherExactTest(k_obs, N_sample - k_obs, k_ref, N_ref - k_ref)

            p = pvalue(t, tail=tail)

            push!(pairs, (motif, p))

        end

        cluster_motifs[index] = pairs
    end

    return cluster_motifs
end


function find_global_motif_pval(g)

    clusters = connected_components(g)
    cluster_motifs = Vector{Vector{Tuple{String,Float64}}}(undef, length(clusters))

    for (index, cluster) in enumerate(clusters)

        in_cluster = falses(nv(g))
        in_cluster[cluster] .= true

        pairs = Tuple{String,Float64}[]

        for u in cluster

            for v in neighbors(g, u)

                if in_cluster[v] && u < v

                    lu = label_for(g, u)
                    lv = label_for(g, v)

                    if haskey(g, lu, lv) && haskey(g[lu, lv], :motifs) && haskey(g[lu, lv], :motif_pvals)

                        motifs = g[lu, lv][:motifs]
                        pvals = g[lu, lv][:motif_pvals]

                        for ind in eachindex(motifs)
                            push!(pairs, (motifs[ind], pvals[ind]))

                        end

                    end
                end
            end
        end

        cluster_motifs[index] = pairs

    end

    return cluster_motifs
end


##############################

# clonal enrichment pval
function find_clone_pval(g, sim_depth)
    clusters = connected_components(g)

    all_clones = Float64[]

    for v in vertices(g)

        lbl = label_for(g, v)
        haskey(g[lbl], :clones) || continue
        append!(all_clones, Float64.(g[lbl][:clones]))

    end

    sample_scores = Float64[]
    cluster_sizes = Int[]

    for cluster in clusters

        labels = label_for.(Ref(g), cluster)
        clones = Float64[]

        for lbl in labels
            haskey(g[lbl], :clones) || continue
            append!(clones, Float64.(g[lbl][:clones]))
        end

        n = length(clones)
        push!(cluster_sizes, n)
        push!(sample_scores, n == 0 ? 0.0 : sum(clones) / n)

    end

    p_vals = Float64[]

    for (index, n) in enumerate(cluster_sizes)

        if n == 0
            push!(p_vals, 1.0)
            continue
        end

        sample_score = sample_scores[index]
        counter = 0

        for i in 1:sim_depth
            sim_clones = sample(all_clones, n; replace=true)
            test_score = sum(sim_clones) / n
            counter += test_score >= sample_score ? 1 : 0
        end

        push!(p_vals, counter == 0 ? 1.0 / sim_depth : counter / sim_depth)

    end

    return p_vals

end

function find_cluster_type(g)
    clusters = connected_components(g)
    types = Vector{Symbol}(undef, length(clusters))

    for (i, cluster) in enumerate(clusters)
        has_distance = false
        has_motifs = false
        saw_edge = false

        # iterate edges induced by the cluster
        for a in 1:length(cluster)-1
            u = cluster[a]
            for b in a+1:length(cluster)
                v = cluster[b]

                has_edge(g, u, v) || continue
                saw_edge = true

                lu = label_for(g, u)
                lv = label_for(g, v)

                if haskey(g[lu, lv], :distance)
                    has_distance = true
                end

                if haskey(g[lu, lv], :motifs)
                    has_motifs = true
                end

                if has_distance && has_motifs
                    types[i] = :mixed
                    @goto done_cluster
                end

            end
        end

        if !saw_edge
            types[i] = :none
        elseif has_distance && !has_motifs
            types[i] = :global
        elseif has_motifs && !has_distance
            types[i] = :local
        else
            # edges exist but neither key found anywhere
            types[i] = :none
        end

        @label done_cluster
    end

    return types
end

function summarize_clusters(
    g,
    min_size::Int;
    sim_depth::Int=1000,
    cdrs_ref=nothing,
    motifs=nothing,
    isglobal::Bool=true,
    islocal::Bool=true,
    nsim::Int=1000,
    rng::AbstractRNG=Random.default_rng(),
    run_length::Bool=true,
    run_vgene::Bool=true,
    run_jgene::Bool=true,
    run_clone::Bool=true,
    run_size::Bool=true,
    cdrs2=nothing,
)

    g = remove_small_clusters(g, min_size)
    clusters, sizes = get_clusters_and_sizes(g)
    n_clusters = length(clusters)

    types = find_cluster_type(g)

    length_pvals = run_length && cdrs2 !== nothing ? find_length_pvals(g, cdrs2, sim_depth) : fill(missing, n_clusters)
    vgene_pvals = run_vgene ? find_vgene_pval(g, sim_depth) : fill(missing, n_clusters)
    jgene_pvals = run_jgene ? find_jgene_pval(g, sim_depth) : fill(missing, n_clusters)
    clone_pvals = run_clone ? find_clone_pval(g, sim_depth) : fill(missing, n_clusters)
    size_pvals = run_size && cdrs_ref !== nothing && motifs !== nothing ?
                 find_size_pval(g, cdrs_ref, motifs; nsim, isglobal, islocal, rng) :
                 fill(missing, n_clusters)

    motif_results = find_global_motif_pval(g)

    rows = []

    for (i, cluster) in enumerate(clusters)

        motif_strs = [pair[1] for pair in motif_results[i]]
        motif_pvals = [pair[2] for pair in motif_results[i]]

        lp = length_pvals[i]
        sp = size_pvals[i]
        cp = clone_pvals[i]

        vp = if vgene_pvals[i] isa Dict
            get(vgene_pvals[i], "vgene_spectratype_p", missing)
        else
            missing
        end

        jp = if jgene_pvals[i] isa Dict
            get(jgene_pvals[i], "jgene_spectratype_p", missing)
        else
            missing
        end

        # combined score — product of all available pvals
        available = filter(!ismissing, [lp, sp, cp, vp, jp])
        if !isempty(motif_pvals)
            append!(available, motif_pvals)
        end
        combined = isempty(available) ? missing : prod(Float64.(available))

        push!(rows, (
            cluster_id=i,
            cluster_size=sizes[i],
            cluster_type=types[i],
            length_pval=lp,
            size_pval=sp,
            clone_pval=cp,
            vgene_pval=vp,
            jgene_pval=jp,
            motifs=isempty(motif_strs) ? missing : motif_strs,
            motif_pvals=isempty(motif_pvals) ? missing : motif_pvals,
            combined_score=combined,
        ))
    end

    return DataFrame(rows)

end

function summarize_cluster_cdr3s(g, min_size::Int)

    g = remove_small_clusters(g, min_size)
    clusters, _ = get_clusters_and_sizes(g)

    rows = []

    for (i, cluster) in enumerate(clusters)

        labels = label_for.(Ref(g), cluster)
        cdr3s = String[]
        indices = Int[]
        barcodes = String[]

        for lbl in labels
            meta = g[lbl]

            push!(cdr3s, String(lbl))

            if haskey(meta, :index)
                push!(indices, meta[:index])
            end

            if haskey(meta, :barcodes)
                bc = meta[:barcodes]
                bc = bc isa AbstractVector ? bc : [bc]
                append!(barcodes, String.(bc))
            end
        end

        push!(rows, (
            cluster_id=i,
            cdr3s=cdr3s,
            row_ids=isempty(indices) ? missing : indices,
            barcodes=isempty(barcodes) ? missing : barcodes,
        ))
    end

    return DataFrame(rows)

end


# plotting

function summarize_cdr3_lengths(df::DataFrame)

    col = nothing
    if "cdr3" in names(df)
        col = :cdr3
    end

    col === nothing && error("No cdr3 column found in DataFrame")

    seqs = filter(!ismissing, df[:, col])
    seqs = filter(x -> x ∉ ("", "NA", "None", "none"), String.(seqs))
    seqs = unique(seqs)

    lengths = length.(seqs)

    tab = sort(DataFrame(countmap(lengths)), :first)
    rename!(tab, :first => :cdr3_length, :second => :count)
    tab[!, :freq] = tab.count ./ sum(tab.count)

    return tab
end

function summarize_umi_distribution(df::DataFrame)

    col = nothing
    if "umis" in names(df)
        col = :umis
    elseif "clones" in names(df)
        col = :clones
    end

    col === nothing && error("No UMI or clone count column found in DataFrame")

    vals = filter(!ismissing, df[:, col])
    vals = Float64.(vals)

    isempty(vals) && error("UMI column is empty")

    result = DataFrame(
        metric=[
            "n",
            "mean",
            "median",
            "std",
            "min",
            "p25",
            "p75",
            "p90",
            "p95",
            "max",
            "total",
        ],
        value=[
            Float64(length(vals)),
            mean(vals),
            median(vals),
            std(vals),
            minimum(vals),
            quantile(vals, 0.25),
            quantile(vals, 0.75),
            quantile(vals, 0.90),
            quantile(vals, 0.95),
            maximum(vals),
            sum(vals),
        ]
    )

    return result
end

function summarize_vgene(df::DataFrame)

    col = nothing
    if "vgene" in names(df)
        col = :vgene
    elseif "v_gene" in names(df)
        col = :v_gene
    else
        return nothing
    end

    vals = filter(!ismissing, df[:, col])
    vals = filter(x -> x ∉ ("", "NA", "None", "none"), String.(vals))
    isempty(vals) && return nothing

    tab = sort(DataFrame(countmap(vals)), :second => rev = true)
    rename!(tab, :first => :vgene, :second => :count)
    tab[!, :freq] = tab.count ./ sum(tab.count)

    return tab
end


function summarize_jgene(df::DataFrame)

    col = nothing
    if "jgene" in names(df)
        col = :jgene
    elseif "j_gene" in names(df)
        col = :j_gene
    else
        return nothing
    end

    vals = filter(!ismissing, df[:, col])
    vals = filter(x -> x ∉ ("", "NA", "None", "none"), String.(vals))
    isempty(vals) && return nothing

    tab = sort(DataFrame(countmap(vals)), :second => rev = true)
    rename!(tab, :first => :jgene, :second => :count)
    tab[!, :freq] = tab.count ./ sum(tab.count)

    return tab
end


function summarize_dgene(df::DataFrame)

    col = nothing
    if "dgene" in names(df)
        col = :dgene
    elseif "d_gene" in names(df)
        col = :d_gene
    else
        return nothing
    end

    vals = filter(!ismissing, df[:, col])
    vals = filter(x -> x ∉ ("", "NA", "None", "none"), String.(vals))
    isempty(vals) && return nothing

    tab = sort(DataFrame(countmap(vals)), :second => rev = true)
    rename!(tab, :first => :dgene, :second => :count)
    tab[!, :freq] = tab.count ./ sum(tab.count)

    return tab
end


function plot_distribution(df::DataFrame; title::String="", top_n::Int=20)

    cols = names(df)

    x_col = nothing
    y_col = :freq

    for c in cols
        if c in ("vgene", "jgene", "dgene", "chain", "cdr3_length", "clonotype_id")
            x_col = Symbol(c)
            break
        end
    end

    x_col === nothing && error("Could not detect category column in DataFrame. Columns: $cols")

    plot_df = nrow(df) > top_n ? first(df, top_n) : df

    plot_df = sort(plot_df, y_col; rev=true)

    x_vals = string.(plot_df[:, x_col])
    y_vals = plot_df[:, y_col]

    fig = Figure(resolution=(max(600, length(x_vals) * 35), 500))
    ax = Axis(
        fig[1, 1];
        title=isempty(title) ? string(x_col) * " distribution" : title,
        xlabel=string(x_col),
        ylabel="frequency",
        xticks=(1:length(x_vals), x_vals),
        xticklabelrotation=π / 3,
        xticklabelsize=11,
    )

    barplot!(ax, 1:length(x_vals), y_vals; color=:steelblue)

    return fig
end

#plot_distribution(summarize_vgene(df); title="V Gene Usage")
#plot_distribution(summarize_jgene(df); title="J Gene Usage")
#plot_distribution(summarize_dgene(df); title="D Gene Usage")
#plot_distribution(summarize_cdr3_lengths(df); title="CDR3 Length Distribution")
#plot_distribution(summarize_chain_usage(df); title="Chain Usage")
