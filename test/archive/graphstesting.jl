# testing for graphs

g = MetaGraph(
    Graph();
    label_type=Symbol,
    vertex_data_type=Vector{String},  # test storing multiple labels?
    edge_data_type=Vector{String},
)

# plans?
# labels will be cdr3s
# vertices will have vector data storing v-gene, j-gene, d-gene (None if no gene), patient,
# edges will store label (global -> distance, local -> motif) and p-val
# maybe store vector of dictionary with String, Any so I can store label followed by thing?
# vertices will store

#vertex_data_type = Dict{Symbol,Any}

add_vertex!(g, :Owo, ["hey", "hello", "hi"])
add_vertex!(g, :Uwu, ["goodbye", "seeya", "bye"])

add_edge!(g, :Uwu, :Owo, ["tidings", "bonjour", "nihao"])

edge_data = g[:Uwu, :Owo]
push!(edge_data, "shalom")
g[:Uwu, :Owo] = edge_data

g[:Uwu, :Owo]


# now trying with dictionary
# Dict{Symbol,Any}

g = MetaGraph(
    Graph();
    label_type=Symbol,
    vertex_data_type=Dict{Symbol,Any},  # test storing multiple labels?
    edge_data_type=Dict{Symbol,Any},
)

add_vertex!(g, :Owo, Dict(:label => "ACGATAG", :vgene => "V13B")) # cdr3 and gene is fake for testing
add_vertex!(g, :Uwu, Dict(:label => "CCGGTAT", :vgene => "V12A"))

add_edge!(g, :Uwu, :Owo, Dict(:type => "local", :pval => 0.02))
add_edge!(g, :Uwu, :Owo, Dict(:type => "global", :distance => 2)) # replaces prior edge

# solution? add a vector of dicionaries to edges and have a global distance, global p-val, motif and motif p-val entry.
# to determine if theres a global or local edge use bool to check for each edge?


# attempt 3

g = MetaGraph(
    Graph();
    label_type=Symbol,
    vertex_data_type=Dict{Symbol,Any},  # test storing multiple labels?
    edge_data_type=Dict{Symbol,Any},
)

add_vertex!(g, :Owo, Dict(:label => "ACGATAG", :vgene => "V13B", :jgene => "J7")) # cdr3 and gene is fake for testing
add_vertex!(g, :Uwu, Dict(:label => "CCGGTAT", :vgene => "V12A", :jgene => "J5"))

add_edge!(g, :Uwu, :Owo, Dict(:distance => 2, :dpval => 0.02, :motif => "TA", :mpval => "0.03"))

g[:Uwu, :Owo][:distance]
g[:Uwu, :Owo][:dpval]

g[:Uwu, :Owo][:new] = "testing" # omg it works!


# _______________________________________ #


function make_vertices!(g, cdrs)
    isblank(x) = x === nothing || x === missing || x in ("None", "none", "NA", "na", "Na", " ", "  ", "")

    for row in nrow(df)

        add_vertex!(g, Symbol(row), Dict(:cdr => cdrs.cdr3[row])) # add vertex

        # check for vgenes, etc. and store if they exist
        if !isblank(cdrs.v_gene[row])
            g[Symbol(row)][:vgene] = cdrs.v_gene[row]
        end

        if !isblank(cdrs.j_gene[row])
            g[Symbol(row)][:jgene] = cdrs.j_gene[row]
        end

        if !isblank(cdrs.d_gene[row])
            g[Symbol(row)][:dgene] = cdrs.d_gene[row]
        end

    end

    return g # return graph

end

# make a global edge w/ annotation of hamming disrance
function add_global_edge!(g, u, v, distance::Int64)
    add_edge!(g, u, v, Dict(:distance => 2))
end

# add a local edge w/ annotation of motif and the pval (from get_significant_motifs)
function add_local_edge!(g, u, v, motif::String, pval::Int64)
    add_edge!(g, u, v, Dict(:motif => motif, :mpval => pval))
end

# go through cdrs
function make_edges(cdrs, motifs)
    g = MetaGraph(Graph(length(cdrs)))
    g = add_vertices(g, cdrs)

    # start w/ global distances
    pairs, dists = make_distance(cdrs.cdr3)
    for (index, (cdr1, cdr2)) in enumerate(pairs)
        d = dists[index]
        if dists[index] <= 1
            add_global_edge!(g, cdr1, cdr2, dists[index])
        end
    end

    # TODO: next do local edges
    for (m, pval) in motifs
        pairs2, hasmotif = make_motif_pairs(cdrs.cdr3, m)

        for (index, (cdr1, cdr2)) in enumerate(pairs2)
            if hasmotif[index] == 1
                add_local_edge!(g, cdr1, cdr2, m, pval)
            end

        end

    end

    return g

end

# only make global edges
function make_global_edges(cdrs, motifs)
    g = MetaGraph(Graph(length(cdrs)))
    g = add_vertices(g, cdrs)

    # start w/ global distances
    pairs, dists = make_distance(cdrs.cdr3)
    for (index, (cdr1, cdr2)) in enumerate(pairs)
        d = dists[index]
        if dists[index] <= 1
            add_global_edge!(g, cdr1, cdr2, dists[index])
        end
    end

    return g

end

# only make local edges
function make_local_edges(cdrs, motifs)
    g = MetaGraph(Graph(length(cdrs)))
    g = add_vertices(g, cdrs)

    # local edge making
    for (m, pval) in motifs
        pairs2, hasmotif = make_motif_pairs(cdrs.cdr3, m)

        for (index, (cdr1, cdr2)) in enumerate(pairs2)
            if hasmotif[index] == 1
                add_local_edge!(g, cdr1, cdr2, m, pval)
            end

        end

    end

    return g

end
