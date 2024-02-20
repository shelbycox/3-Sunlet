include("sampling.jl")
include("lambda_cyclic.jl")

## get the following dictionary: keys = 0-1 vectors, values = (neighbors, rank)
function gen_graph_poset(group::FiniteCyclicGroup, N::Int64)
    H = generateSunletArr(group)
    S = sample(group, H, N)
    regions = collect(keys(S))
    return Dict(k => poset_data(k, S[k], regions, group) for k in regions)
end

function poset_data(I, sample_point, regions, group)
    neighbors = regions[map(regions -> areNeighbors(regions, I), regions)]

    zeroRegion = [Bool(0) for i=1:length(I)]
    distToZero = distToRegion(I, zeroRegion, regions)

    r = getRank(sample_point, group)
    
    return (neighbors, distToZero, r)
end

NODE_FORMAT = "node [shape=point];"
COLORS = ["blue", "red", "yellow", "purple", "green", "black", "goldenrod1"]
SHAPES = ["rectangle", "triangle", "circle", "egg", "pentagon", "hexagon", "pentagon"]

function toDot(group::FiniteCyclicGroup, N::Int64)
    P = gen_graph_poset(group, N)
    regions = collect(keys(P))
    ranks = sort(collect(Set(P[I][3] for I in regions)))
    dists = sort(collect(Set(P[I][2] for I in regions)))

    output_file = Base.string("graphs-DOT/", string(group), ".dot")
    
    open(output_file, "w+") do to_write
        write(to_write, "strict graph G {\n")
        # write(to_write, Base.string("\t", NODE_FORMAT, "\n")) ## set node shape to dot

        for d in dists
            ddist = join([i for i=eachindex(regions) if P[regions[i]][2] == d], ' ')
            for i=eachindex(ranks)
                r = ranks[i]
                drrank = [i for i=eachindex(regions) if (P[regions[i]][3] == r) & (P[regions[i]][2] == d)]
                if length(drrank) > 0
                    drrank_string = join(drrank, ' ')
                    write(to_write, "\tnode [shape=point, color=$(COLORS[i]), shape=$(SHAPES[i])];\n\t{ $drrank_string };\n")
                    # write(to_write, "\t{group=G$d$r $drrank_string};\n") ## doesn't seem to do anything :(
                end
            end
            write(to_write, "\t{ rank=same $ddist };\n")
        end

        for i=eachindex(regions)
            R = regions[i]
            targets = join([j for j=eachindex(regions) if areNeighbors(R, regions[j])], ' ')
            write(to_write, "\t$i -- {$targets};\n")
        end
        
        ## rank nodes by dist to zero (vertical) and rank (horizontal)
    
        write(to_write, "}")
    end
end