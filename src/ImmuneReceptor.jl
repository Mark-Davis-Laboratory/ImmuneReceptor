module ImmuneReceptor

const PK = pkgdir(ImmuneReceptor)

const IN = joinpath(PK, "in")

const OU = joinpath(PK, "ou")

# ----------------------------------------------------------------------------------------------- #

using ProgressMeter: @showprogress

using Nucleus

# to do:
# -> add index tracking to the motif function ✅
# -> add a way to check for index / window overlap ✅
# -> add implementation of blosum function in the make edges
# -> combine blosum functions into one ✅
# -> create function to make like a nice summary file
# -> implement fisher's exact test for cluster significance
# -> edit scoring functions to either consistently use vector of dictionary or tuple

# =============================================================================================== #
# Reading
# =============================================================================================== #

function is_t_gene(::Any)

    false

end

function is_t_gene(st::AbstractString)

    startswith(st, 'T')

end

function is_cdr3(st)

    st[1] == 'C' && st[end] == 'F'

end

function load_cdr3s(csvpath)
    df = CSV.read(csvpath, DataFrame)
    chain_keep_a = df.chain .∈ Ref(["TRG", "TRD", "TRA", "TRB"])
    chain_keep_b = .!(df.chain .∈ Ref(["None", "none", "NA", "na", " ", "Multi"]))
    cdr3_keep  = .!(df.cdr3 .∈ Ref(["None", "none", "NA", "na", " "]))
    c_check = startswith.(df.cdr3, "C")
    f_check = endswith.(df.cdr3, "F")

    mask = c_check .& chain_keep_a .& chain_keep_b .& cdr3_keep .& f_check

    filtered_df = df[mask, :]

    filtered_df_2 = transform(groupby(filtered_df, :cdr3), nrow => :duplicate_count)
    return filtered_df_2
end

# =============================================================================================== #
# Writing
# =============================================================================================== #

function make_trace(an_)

    um = 1000000

    Dict(
        "name" => "Naive",
        "type" => "histogram",
        "histnorm" => "probability",
        "x" => lastindex(an_) <= um ? an_ : rand(an_, um),
    )

end

function writ(ht, st, a1_, a2_)

    Nucleus.Plotly.writ(
        ht,
        (make_trace(a1_), make_trace(a2_)),
        Dict(
            "yaxis" => Dict("title" => Dict("text" => "Probability")),
            "xaxis" => Dict("title" => Dict("text" => st)),
        ),
    )

end

# =============================================================================================== #
# MISC
# =============================================================================================== #

function make_dictionary(an_)

    Dict(an_[nd] => nd for nd in eachindex(an_))

end

# =============================================================================================== #
# CDR3
# =============================================================================================== #

function make_hamming_distance(s1, s2)

    sum(s1[nd] != s2[nd] for nd in eachindex(s1))

end

const blosum62 = Dict{Tuple{Char,Char},Int}(
    ('A','A')=>4,  ('A','R')=>-1, ('A','N')=>-2, ('A','D')=>-2, ('A','C')=>0,
    ('A','Q')=>-1, ('A','E')=>-1, ('A','G')=>0,  ('A','H')=>-2, ('A','I')=>-1,
    ('A','L')=>-1, ('A','K')=>-1, ('A','M')=>-1, ('A','F')=>-2, ('A','P')=>-1,
    ('A','S')=>1,  ('A','T')=>0,  ('A','W')=>-3, ('A','Y')=>-2, ('A','V')=>0,

    ('R','R')=>5,  ('R','N')=>0,  ('R','D')=>-2, ('R','C')=>-3, ('R','Q')=>1,
    ('R','E')=>0,  ('R','G')=>-2, ('R','H')=>0,  ('R','I')=>-3, ('R','L')=>-2,
    ('R','K')=>2,  ('R','M')=>-1, ('R','F')=>-3, ('R','P')=>-2, ('R','S')=>-1,
    ('R','T')=>-1, ('R','W')=>-3, ('R','Y')=>-2, ('R','V')=>-3,

    ('N','N')=>6,  ('N','D')=>1,  ('N','C')=>-3, ('N','Q')=>0,  ('N','E')=>0,
    ('N','G')=>0,  ('N','H')=>1,  ('N','I')=>-3, ('N','L')=>-3, ('N','K')=>0,
    ('N','M')=>-2, ('N','F')=>-3, ('N','P')=>-2, ('N','S')=>1,  ('N','T')=>0,
    ('N','W')=>-4, ('N','Y')=>-2, ('N','V')=>-3,

    ('D','D')=>6,  ('D','C')=>-3, ('D','Q')=>0,  ('D','E')=>2,  ('D','G')=>-1,
    ('D','H')=>-1, ('D','I')=>-3, ('D','L')=>-4, ('D','K')=>-1, ('D','M')=>-3,
    ('D','F')=>-3, ('D','P')=>-1, ('D','S')=>0,  ('D','T')=>-1, ('D','W')=>-4,
    ('D','Y')=>-3, ('D','V')=>-3,

    ('C','C')=>9,  ('C','Q')=>-3, ('C','E')=>-4, ('C','G')=>-3, ('C','H')=>-3,
    ('C','I')=>-1, ('C','L')=>-1, ('C','K')=>-3, ('C','M')=>-1, ('C','F')=>-2,
    ('C','P')=>-3, ('C','S')=>-1, ('C','T')=>-1, ('C','W')=>-2, ('C','Y')=>-2,
    ('C','V')=>-1,

    ('Q','Q')=>5,  ('Q','E')=>2,  ('Q','G')=>-2, ('Q','H')=>0,  ('Q','I')=>-3,
    ('Q','L')=>-2, ('Q','K')=>1,  ('Q','M')=>0,  ('Q','F')=>-3, ('Q','P')=>-1,
    ('Q','S')=>0,  ('Q','T')=>-1, ('Q','W')=>-2, ('Q','Y')=>-1, ('Q','V')=>-2,

    ('E','E')=>5,  ('E','G')=>-2, ('E','H')=>0,  ('E','I')=>-3, ('E','L')=>-3,
    ('E','K')=>1,  ('E','M')=>-2, ('E','F')=>-3, ('E','P')=>-1, ('E','S')=>0,
    ('E','T')=>-1, ('E','W')=>-3, ('E','Y')=>-2, ('E','V')=>-2,

    ('G','G')=>6,  ('G','H')=>-2, ('G','I')=>-4, ('G','L')=>-4, ('G','K')=>-2,
    ('G','M')=>-3, ('G','F')=>-3, ('G','P')=>-2, ('G','S')=>0,  ('G','T')=>-2,
    ('G','W')=>-2, ('G','Y')=>-3, ('G','V')=>-3,

    ('H','H')=>8,  ('H','I')=>-3, ('H','L')=>-3, ('H','K')=>-1, ('H','M')=>-2,
    ('H','F')=>-1, ('H','P')=>-2, ('H','S')=>-1, ('H','T')=>-2, ('H','W')=>-2,
    ('H','Y')=>2,  ('H','V')=>-3,

    ('I','I')=>4,  ('I','L')=>2,  ('I','K')=>-3, ('I','M')=>1,  ('I','F')=>0,
    ('I','P')=>-3, ('I','S')=>-2, ('I','T')=>-1, ('I','W')=>-3, ('I','Y')=>-1,
    ('I','V')=>3,

    ('L','L')=>4,  ('L','K')=>-2, ('L','M')=>2,  ('L','F')=>0,  ('L','P')=>-3,
    ('L','S')=>-2, ('L','T')=>-1, ('L','W')=>-2, ('L','Y')=>-1, ('L','V')=>1,

    ('K','K')=>5,  ('K','M')=>-1, ('K','F')=>-3, ('K','P')=>-1, ('K','S')=>0,
    ('K','T')=>-1, ('K','W')=>-3, ('K','Y')=>-2, ('K','V')=>-2,

    ('M','M')=>5,  ('M','F')=>0,  ('M','P')=>-2, ('M','S')=>-1, ('M','T')=>-1,
    ('M','W')=>-1, ('M','Y')=>-1, ('M','V')=>1,

    ('F','F')=>6,  ('F','P')=>-4, ('F','S')=>-2, ('F','T')=>-2, ('F','W')=>1,
    ('F','Y')=>3,  ('F','V')=>-1,

    ('P','P')=>7,  ('P','S')=>-1, ('P','T')=>-1, ('P','W')=>-4, ('P','Y')=>-3,
    ('P','V')=>-2,

    ('S','S')=>4,  ('S','T')=>1,  ('S','W')=>-3, ('S','Y')=>-2, ('S','V')=>-2,

    ('T','T')=>5,  ('T','W')=>-2, ('T','Y')=>-2, ('T','V')=>0,

    ('W','W')=>11, ('W','Y')=>2,  ('W','V')=>-3,

    ('Y','Y')=>7,  ('Y','V')=>-1,

    ('V','V')=>4
)

function make_blosum_score(s1, s2)

    @assert length(s1) == length(s2)

    mismatch_indices = Int[]

    for nd in eachindex(s1)
        if s1[nd] != s2[nd]
            push!(mismatch_indices, nd)
        end
    end

    if isempty(mismatch_indices) || length(mismatch_indices) > 1
            return length(mismatch_indices)
    end

    nd = mismatch_indices[1]

    a1 = s1[nd]
    a2 = s2[nd]

    scores = blosum62[(a1,a2)], blosum62[(a2,a1)]

    return min(scores)

end

function make_distance(st_)

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

        in__[i1+=1] = i2, i3

        po_[i1] = make_hamming_distance(s1, s2)

    end

    resize!(in__, i1), resize!(po_, i1)

end

# =============================================================================================== #
# Motif
# =============================================================================================== #

# count whether motif is present at least once per cdr3 for all cdr3s

function get_motif_counts(
    motif::AbstractVector{<:AbstractString},
    cdrs::AbstractVector{<:AbstractString},
)

    motif_counts = Dict{String,Int}()

    for m in motif

        motif_counts[m] = count(cdr -> occursin(m, cdr), cdrs) # count occurances per cdr3 for each motif

    end

    return motif_counts

end

# makes list of motifs in cdr3s, counts them, filters based on a cutoff, then provides set count
function get_motifs(st_::AbstractVector{<:AbstractString}, min::Int, max::Int)

    mo_ = Dict{String,Int}() # make dictionary

    @showprogress desc = "Finding and counting motifs..." for s1 in st_

        lastindex(s1) < 7 && continue # only keep cdr3 7 or longer

        s2 = s1[4:(end-3)] # remove first and last 3 aa
        for um in min:max

            um > lastindex(s2) && continue # make sure cdr3 is longer than the motif size

            i1 = 0
            i2 = i1 + um - 1

            while i2 < lastindex(s2)

                m = s2[(i1+=1):(i2+=1)] # get the motif
                mo_[m] = get!(mo_, m, 0) + 1 # count total motif occurances

            end

        end

    end

    mo_ = Dict(m => num for (m, num) in mo_ if num >= 3) # cutoff of 3
    mo2_ = get_motif_counts(collect(keys(mo_)), st_)
    return mo2_

end

function find_significant_motifs(motifs, cdrs1, cdrs2, nsim, ove_cutoff)

    Random.seed!(1)

    motifs_list = collect(keys(motifs))
    L = length(motifs_list)
    counts_orig = [motifs[m] for m in motifs_list]

    counts_sim = Array{Float64}(undef, nsim, length(motifs_list))

    significant_motifs = Dict{String,Float64}()

    @showprogress desc = "Calculating significant motifs..." for i in 1:nsim

        random_cdrs = sample(cdrs2, length(cdrs1); replace=true, ordered=false)
        random_counts = get_motif_counts(motifs_list, random_cdrs)
        for j in 1:L
            counts_sim[i, j] = get(random_counts, motifs_list[j], 0)
        end

    end

    # everything above this line works

    for (index, m) in enumerate(motifs_list)

        if ove_cutoff == true
            ove = counts_orig[index] / mean(counts_sim[:, index])

            if counts_orig[index] < 2
                continue

            elseif (counts_orig[index] == 2 && ove >= 1000) ||
                (counts_orig[index] == 3 && ove >= 100) ||
                (counts_orig[index] >= 4 && ove >= 10)

                wins = count(x -> x >= counts_orig[index], counts_sim[:, index])
                p_val = (wins + 1) / (nsim + 1)
                #print(p_val)
                if p_val <= 0.05
                    significant_motifs[m] = get!(significant_motifs, m, 0.0) + p_val
                end

            end
        end

        if ove_cutoff != true
            wins = count(x -> x >= counts_orig[index], counts_sim[:, index])
            p_val = (wins) / (nsim)
            print(p_val)
            if p_val <= 0.05
                significant_motifs[m] = get!(significant_motifs, m, 0.0) + p_val
            end
        end

    end

    println(string("Number of significant motifs identified:", length(significant_motifs)))
    return significant_motifs

end


#_________#

# takes list of cdrs and checks for a motif, making pairs
function make_motif_pairs(st_::AbstractVector{<:AbstractString}, motif)
    u1 = lastindex(st_)
    u2 = div(u1 * (u1 - 1), 2)

    in__ = Vector{Tuple{Int,Int}}(undef, u2)
    po_ = Vector{Int}(undef, u2)

    i1 = 0

    for i2 in 1:u1, i3 in (i2+1):u1
        s1, s2 = st_[i2], st_[i3]

        has1 = occursin(motif, s1)
        has2 = occursin(motif, s2)

        if has1 && has2
            i1 += 1
            in__[i1] = (i2, i3)
            po_[i1] = 1
        end
    end

    resize!(in__, i1)
    resize!(po_, i1)

    return in__, po_
end

function make_motif_pairs_new(st_::AbstractVector{<:AbstractString}, motif)
    u1 = lastindex(st_)
    u2 = div(u1 * (u1 - 1), 2)
    l = length(motif)

    in__ = Vector{Tuple{Int,Int}}(undef, u2)
    po_ = Vector{Int}(undef, u2)

    i1 = 0

    for i2 in 1:u1, i3 in (i2+1):u1
        s1, s2 = st_[i2], st_[i3]

        has1 = occursin(motif, s1)
        has2 = occursin(motif, s2)

        if has1 && has2

            nd1 = findfirst(motif, s1)
            nd2 = findfirst(motif, s2)

            a1, a2 = start(nd1), stop(nd1)
            b1, b2 = start(nd2), stop(nd2)

            if max(a1, b1) ≤ min(a2, b2) + 3

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

# =============================================================================================== #
# Graphing
# =============================================================================================== #

function make_vertices!(g, cdrs)
    isblank(x) = x === nothing || x === missing || x in ("None", "none", "NA", "na", "Na", " ", "  ", "")

    for row in 1:nrow(cdrs)

        label = String(cdrs.cdr3[row])
        add_vertex!(g, label, Dict(:index => row)) # add vertex

        # check for vgenes, etc. and store if they exist
        if !isblank(cdrs.v_gene[row])
            g[label][:vgene] = cdrs.v_gene[row]
        end

        if !isblank(cdrs.j_gene[row])
            g[label][:jgene] = cdrs.j_gene[row]
        end

        if !isblank(cdrs.d_gene[row])
            g[label][:dgene] = cdrs.d_gene[row]
        end

        if !isblank(cdrs.donor[row])
            g[label][:donor] = cdrs.donor[row]
        end

    end

    return g # return graph

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

    cdr3_vec = String.(cdrs.cdr3)

    g = MetaGraph(
        Graph();
        label_type=String,
        vertex_data_type=Dict{Symbol,Any},  # test storing multiple labels?
        edge_data_type=Dict{Symbol,Any},
    )

    g = make_vertices!(g, cdrs)

    # start w/ global distances
    # # changed so label is cdr3
    if isglobal == true
        pairs, dists = make_distance(cdr3_vec)

            @showprogress desc = "Making edges..." for (index, (i, j)) in enumerate(pairs)
                d = dists[index]
                if d <= 1
                    # i, j are indices into cdr3_vec
                    u = cdr3_vec[i]
                    v = cdr3_vec[j]

                    if haskey(g, u, v)
                        g[u, v][:distance] = d
                    else
                        add_global_edge!(g, u, v, d)
                    end
                end
            end
    end

    # TODO: next do local edges
    if islocal == true

       @showprogress desc = "Making edges..." for (motif, pval) in motifs

            mask  = occursin.(Ref(motif), cdr3_vec)
            matches  = findall(mask)
            nmatches = length(matches)

            # nothing to do if fewer than 2 sequences have this motif
            if nmatches <= 1
                continue
            end

            # all unique pairs of those indices
            for i in 1:(nmatches-1), j in (i+1):nmatches
                u = cdr3_vec[matches[i]]
                v = cdr3_vec[matches[j]]

                if haskey(g, u, v)
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


# =============================================================================================== #
# Clusters / Significance
#
# TO DO: multiple comparisons correction?
# =============================================================================================== #

# length score
# for each group, resample 1000 groups of group size n and measure cdr3 length distribution.

using StatsBase: mode

function find_length_pvals(g, cdrs2, sim_depth)
    clusters = connected_components(g)

    modes = Vector{Int}()
    props = Vector{Float64}()

    # get most common length for each cluster
    for cluster in clusters

        labels = label_for.(Ref(g), cluster) # should be fixed now!
        cluster_lengths = [length(s) for s in labels]
        m = mode(cluster_lengths)
        push!(modes, m)
        push!(props, count(==(m), cluster_lengths) / length(cluster))

    end

    # randomly pick n sequences from the null distribution and compare to
    cluster_sizes = [length(cluster) for cluster in clusters]
    p_vals = Float64[]

    for (index, size) in enumerate(cluster_sizes)

        sim_props = Float64[]

        for i in 1:sim_depth
            random_cdrs = sample(cdrs2, size; replace=true, ordered=false)
            cluster_lengths = [length(s) for s in random_cdrs]
            push!(sim_props, (count(==(modes[index]), cluster_lengths)) / size)
        end

        # count number of times sim beats data?
        p_val = ((count(>=(props[index]), sim_props)) + 1) / (sim_depth + 1)
        push!(p_vals, p_val)

    end

    return p_vals

end

function score_lengths(g, cdrs2)

    clusters = connected_components(g)

    p_vals = find_length_pvals(g, cdrs2)

    for (index, cluster) in enumerate(clusters)

        labels = label_for.(Ref(g), cluster)

        for label in labels
            g[label][:length_pval] = p_vals[index] # adding cluster pval to each vertex
        end

    end

    return p_vals

end

#_______________#

# v gene score
# hyper geomtric or parameteric when over 200 sequences in one group

using Distributions

function score_vgene(g, sim_depth)
    clusters = connected_components(g)

    # make vector w/ dictionary of vgenes
    counts = Vector{Dict{String,Int}}(undef, length(clusters))
    totals = Dict{String,Int}()
    sizes = Vector{Int}(undef, length(clusters))

    for (index, cluster) in enumerate(clusters)
        di = Dict{String,Int}()
        sizes[index] = length(cluster)
        for vertex in cluster
            vgene = g[label_for(g, vertex)][:vgene] # should be fixed now!
            di[vgene] = get!(di, vgene, 0) + 1
            totals[vgene] = get!(totals, vgene, 0) + 1
        end
        counts[index] = di
    end

    # actually do p-vals after counting
    v_all = [g[cdr][:vgene] for cdr in labels(g) if haskey(g[cdr], :vgene)]
    cluster_pvals = Vector{Dict{String,Float64}}(undef, length(clusters))

    for (index, di) in enumerate(counts)

        p_vals = Dict{String,Float64}()

        if sizes[index] <= 200 # for less than 200 in a cluster

            # TO DO = verify if below is correct

            for (vgene, count_) in di
                distribution = Hypergeometric(totals[vgene], sum(sizes) - totals[vgene], sizes[index])
                p = ccdf(distribution, count_ - 1)
                p_vals[vgene] = p
            end

        else # for when > 200 in a cluster

            sim_wins = Dict{String,Int}()

            for i in 1:sim_depth # do sim_depth random pulls from the samples and count vgenes

                # sizes[index] has cluster size - stored in earlier loop
                random_vgenes = sample(v_all, sizes[index]; replace=true, ordered=false)

                # go through every vgene
                for (vgene, orig_count) in di

                    count_ = count(x -> x == vgene, random_vgenes)

                    if count_ >= orig_count
                        sim_wins[vgene] = get!(sim_wins, vgene, 0) + 1
                    else
                        sim_wins[vgene] = get!(sim_wins, vgene, 0) + 0
                    end

                end

            end

            # go through sim wins + gen new dictionary

            for (vgene, wins) in sim_wins
                p = (wins + 1) / (sim_depth + 1)
                p_vals[vgene] = p
            end

        end

        val, key = findmin(p_vals)   # val = pval (Float64), key = vgene (String)
        cluster_pvals[index] = Dict(key => (val * length(p_vals)))

    end

    # add cluster vgene pvals to graph
    for (index, cluster) in enumerate(clusters)

        labels = label_for.(Ref(g), cluster)

        for label in labels
            g[label][:vgene_pval] = cluster_pvals[index] # adding cluster pval to each vertex
        end

    end

    return cluster_pvals

end

#_______________#

# hla score

# for each CDR, get donor
# check hla file for shared hlas
# do enrichment exactly as done for v-genes

function allele_group(allele::String)
    m = match(r"^(HLA-)?([A-Z]{1,4}\d?)", allele)
    return m === nothing ? allele : m.captures[end]
end

function score_hla(g)
    clusters = connected_components(g, hla_df)

    # make dictionaries and vectors to store stuff in ✅
    counts = Vector{Dict{String,Int}}(undef, length(clusters))

    totals = Dict{String,Int}()
    sizes = Vector{Int}(undef, length(clusters))

    for (index, cluster) in enumerate(clusters) # count things ✅

        di = Dict{String,Int}()
        sizes[index] = length(cluster)

        for vertex in cluster

            donor = g[label_for(g, vertex)][:donor]

            row_index = findfirst(==(donor), hla_df.donor)
            row_index === nothing && continue # skip typing for donors without HLA type
            row = hla_df[row_index, 2:end]

            for allele in row

                if ismissing(allele)

                    continue

                end

                a = String(allele) # needs to be string

                totals[allele] = get(totals, allele, 0) + 1 # store total counts for distribution to compare cluster to
                di[allele] = get!(di, allele, 0) + 1

            end


        end

        counts[index] = di # store individual counts for each cluster

    end

    cluster_pvals = Vector{Dict{String,Float64}}(undef, length(clusters))

    for (index, di) in enumerate(counts)

        counts[index]
        p_vals = Dict{String,Float64}()

        for (hla_, count_) in di

            distribution = Hypergeometric(totals[index][hla_], sum(sizes) - totals[index][hla_], sizes[index])
            p = ccdf(distribution, count_ - 1)
            p_vals[hla_] = p

        end

        cluster_pvals[index] = p_vals

    end

    # generate table

    pval_df = DataFrame(cluster=Int[], allele=String[], allele_group=String[], pval=Float64[])

    for (cluster, d) in pairs(cluster_pvals)
        for (allele, p) in d
            group = allele_group(allele)
            push!(pval_df, (
                cluster = cluster,
                allele  = allele,
                group   = group,
                pval    = p
            ))
        end
    end

    return pval_df
    # heatmap_df = unstack(df, :cluster, :allele, :pval)

end


# motif scoring

using Graphs

function collect_cluster_motifs(g)
    clusters = connected_components(g)
    cluster_motifs = Vector{Vector{Tuple{String,Float64}}}(undef, length(clusters))

    for (index, cluster) in enumerate(clusters)
        in_cluster = falses(nv(g))
        in_cluster[cluster] .= true # find cluster members

        pairs = Tuple{String,Float64}[]

        for u in cluster

            for v in neighbors(g, u)

                if in_cluster[v] && u < v   # avoid double counting

                    if haskey(g[u, v], :motifs) && haskey(g[u, v], :motif_pvals)

                        motifs = g[u, v][:motifs]
                        pvals  = g[u, v][:motif_pvals]

                        @assert length(motifs) == length(pvals)

                        for ind in eachindex(motifs) # store paired motif and pval

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




# =============================================================================================== #
# Summarising
# =============================================================================================== #
#

# create function to report for each cluster:
# -> number of members
# -> number of unique donors and unique clones (unique alpha and beta?) - less of a priority
# -> length of cdr3 p-val
# -> motif p-val
# -> v-gene pval
# -> hla p-vals
# -> list of probable hla types

# length -> vector of p-vals
# indexed dictionary with correct p-val of most significant v-gene

function summarize_data(g, vgene_pvals, length_pvals)

    # get cluster sizes

    clusters = connected_components(g)

    # things to count:
    # unique members, unique cdr3s, (later - unique clones)
    # find most significant motif p-val

    # make vector w/ dictionary of vgenes
    cluster_sizes = Vector{Dict{String,Int}}(undef, length(clusters))
    totals = Dict{String,Int}()
    cluster_sizes = Vector{Int}(undef, length(clusters))
    donor_sizes = Vector{Int}(undef, length(clusters))
    length_pval = Vector{Float64}(undef, length(clusters))

    # get cluster sizes and no. donors
    for (index, cluster) in enumerate(clusters)
        donors = Vector{String}()
        cluster_sizes[index] = length(cluster)

        for vertex in cluster
            donor = g[label_for(g, vertex)][:donor] # should be fixed now!
            donors[vertex] = donor
        end

        donor_sizes[index] = length(unique(donors))

    end

    vgenes = [first(keys(i)) for i in vgene_pvals]
    v_pvals  = [first(values(i)) for i in vgene_pvals]

    df = DataFrame(
        group = Int64[],
        no_member = cluster_sizes,
        no_donor = donor_sizes,
        length_pval = length_pvals,
        vgene_pval = v_pvals,
        vgene = vgenes,
        motif_pval = Int64[],
        #motif = String[]
    )


    # create DataFrame
    # fill in iteratively with everything
    #



end
