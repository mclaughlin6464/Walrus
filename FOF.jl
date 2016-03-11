# Friends of Friends clustering.  (Tom Abel, 1/2016)


#using FixedLengthVectors: FixedLengthVector

"""
Friends of Friends Algorithm: for sets of points x, using linking length l and returning groups with more than minlength members in an array of arrays of indeces inot the point data provided.
    using Gadfly
    x = rand((2,100))
    gps = FOF.groups(x,.1, 3)
    n = 1
    plot(layer(x=x[1,gps[n]], y=x[2,gps[n]], Geom.point,Theme(default_color=colorant"red")), layer(x=x[1,:], y=x[2,:], Geom.point, Theme(default_color=colorant"blue"))  )

"""
#TODO verbose option
function groups(x, l, minlength, v = nothing, l_v = nothing)
    if (v == nothing && l_v != nothing) || (v != nothing && l_v == nothing)
        throw(ArgumentError("Both v and l_v must be specified"))
    end

    Npart = size(x)[2]
    tree = KDTree(x) # Build tree with point data provided: KDtree is twice as fast as a balltree
    if v != nothing
        tree_v = KDTree(v)
    end
    println("FOF: built tree")
    ds = IntDisjointSets(Npart)
    for i in 1:Npart
        idxs = IntSet(inrange(tree, x[:,i], l, false)) # within search radius
        #TODO no velcoties or positions are linked?
        if v != nothing
            #println(idxs)
            idxs_v = IntSet(inrange(tree_v, v[:,i], l_v, false))
            #println(idxs_v)
            #println(size(idxs_v))
            #changes idxs to an IntSet. Probably fine.
            intersect!(idxs, idxs_v)

            #println(idxs)
            #println("\n")
        end
        for idx in idxs #
            union!(ds,i,idx)
        end
        if (num_groups(ds) == 1) # just in case people use too large a linking length don't waste time
            println("FOF: All points were linked. Exiting." )
            break
        end  # in case everything has been joined already exit
    end
    println("FOF: finished grouping")

    idxs = find(ds.ranks) # all non-zero ranks are parent particles in groups
    groupid = [ds.parents[idxs[i]] => i  for i in eachindex(idxs)]
    grouplen = Dict{Int,Int}()
    for i in 1:Npart
        if get(groupid, ds.parents[i], 0) > 0
            grouplen[ds.parents[i]] = get(grouplen, ds.parents[i], 0) + 1
        end
    end

    # now we collect the actual particles in the groups of the length we are interested in
    for (k,v) in grouplen
        if (v < minlength)
            delete!(grouplen, k)
        end
    end

    Ngroups = length(grouplen)
    # and provide them in reverse order with the biggest group first
    sid = sort(collect(grouplen), by = tuple -> last(tuple),rev=true)
    grouplo = [sid[i].first => i for i in 1:length(sid)]
    gps = [Int[] for i in 1:Ngroups]
    for i in 1:Npart
        if get(grouplen, ds.parents[i], 0) > 0
            push!(gps[grouplo[ds.parents[i]]], i)
        end
    end

    println("FOF: Found ", Ngroups, " with ", minlength, " or more points")
    gps # An array of different  sized arrays containing the particle ids in the groups is returned
end

#TODO ranks, parents. I want this to work with the above api. 
type ID_Set#datatype that is similar to IntDisjointSets but more general
    #doesn't require sequential integers
    main_dict::Dictionary
    ngroups::Int

    function ID_Set(I::Array{Int, 1})
        main_dict = Dict{Int,Set}(i => Set(i) for i in I)
        ngroups = size(I,1)
        new(main_dict, ngroups)
    end
end

function union!(ds::ID_Set, i1::Int, i2::Int)
    if i2 in ds.main_dict[i1]
        return nothing
    end

    push!(ds.main_dict[i1], i2)
    ds.main_dict[i2] = ds.main_dict[i1]
    ngroups-=1
    nothing
end

function num_groups(ds::ID_Set) = ID_Set.ngroups



#NOTE should I have it do all dims? Need that for 1 box but not several
function link_boundaries(P_box1::Particles, P_box2::Particles, l::Float64, ds::IntDisjointSets)
    for dim in 1:3
        box_correction = zeros(3)
        box_correction[dim] = BoxSize
        for p1 in P_box1[P_box1.x[:,dim]< l]
            for p2 in P_box2[P_box2.x[:, dim]> BoxSize - l]
                if in_same_set(ds, p1, p2)#NOTE this will fail.
                    #set needs to track id's; not implemented!
                    continue
                elseif euclidian(p1.x, p2.x+box_correction) < l
                    union!(ds, p1, p2)#NOTE Not implemented!
                end
            end
        end
    end
end
