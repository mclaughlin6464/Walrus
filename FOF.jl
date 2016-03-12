# Friends of Friends clustering.  (Tom Abel, 1/2016)

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

        #add boundary conditions here
        boundary_idxs = link_boundaries(tree, x[:,i], l)
        union!(idxs, boundary_idxs)

        if v != nothing
            idxs_v = IntSet(inrange(tree_v, v[:,i], l_v, false))
            #changes idxs to an IntSet. Probably fine.
            intersect!(idxs, idxs_v)
        end
        for idx in idxs #
            union!(ds,i,idx)
        end
        if (num_groups(ds) == 1) # just in case people use too large a linking length don't waste time
            println("FOF: All points were linked. Exiting." )
            break
        end  # in case everything has been joined already exit
    end
    """
    println("FOF: Doing boundaries")
    #take care of periodic boundary conditions
    if num_groups(ds) != 1
        for dim in 1:size(x,1)
            if v != nothing
                @time link_boundaries!(ds, dim, x, x,l, v, v, l_v)
            else
                link_boundaries!(ds, dim, x, x, l)
            end
        end
    end
    println("FOF: finished grouping")
    """
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

function link_boundaries(tree::KDTree, dim::Int, x::Array{Float64,1}, l::Real)
    if x[dim] > l
        return IntSet()
    end

    BoxCorrection = zeros(size(x,1))
    BoxCorrection[dim] = BoxSize
    idxs = IntSet(inrange(tree, x+BoxCorrection, l, false)) # within search radius
    return idxs
end

function link_boundaries(tree::KDTree, x::Array{Float64,1}, l::Real)
    Ndim = size(x,1)
    output_set = IntSet()

    for dim in 1:Ndim
        idxs = link_boundaries(tree, dim, x, l)
        union!(output_set, idxs)
    end
    return output_set
end

#TODO combine these 2 methods!
#TODO not efficient
function link_boundaries!(ds::IntDisjointSets, dim::Int, x1::Array{Float64,2}, x2::Array{Float64, 2},l::Float64)
    if dim < 0 || dim > 4
        throw(ArgumentError("dim must be between 1 and 3 inclusive."))
    end
    NP1 = size(x1,2)
    NP2 = size(x2,2)
    for i in 1:NP1
        if x1[dim, i] > l
            continue
        end
        for j in 1:NP2
            if in_same_set(ds, i, j)
                continue
            elseif x2[dim, j] < BoxSize-l
                continue
            elseif periodic_euclidean(x1[:,i], x2[:,j], BoxSize)[1] < l
                union!(ds, i,j)
            end
        end
    end
    nothing
end

link_boundaries!(ds::IntDisjointSets, dim::Int, x1::Array{Float64,2}, x2::Array{Float64, 2},l::Int) = link_boundaries!(ds, dim,x1, x2, convert(Float64, l))

function link_boundaries!(ds::IntDisjointSets, dim::Int, x1::Array{Float64,2}, x2::Array{Float64, 2},l::Float64, v1::Array{Float64,2}, v2::Array{Float64,2}, l_v::Float64)
    if dim < 0 || dim > 4
        throw(ArgumentError("dim must be between 1 and 3 inclusive."))
    end
    NP1 = size(x1,2)
    NP2 = size(x2,2)
    for i in 1:NP1
        if x1[dim, i] < l
            continue
        end
        for j in 1:NP2
            if in_same_set(ds, i, j)
                continue
            elseif x2[dim, j] > BoxSize-l
                continue
            elseif euclidean(v1[:,i],v2[:,j]) > l_v
                continue
            elseif periodic_euclidean(x1[:,i], x2[:,j], BoxSize)[1] < l
                union!(ds, i,j)
            end
        end
    end
    nothing
end

link_boundaries!(ds::IntDisjointSets, dim::Int, x1::Array{Float64,2}, x2::Array{Float64, 2},l::Int, v1::Array{Float64,2}, v2::Array{Float64,2}, l_v::Int) = link_boundaries!(ds, dim,x1, x2, convert(Float64, l), v1, v2, convert(Float64, l_v))


#TODO find a better place to put htis.
function periodic_euclidean(a::AbstractArray, b::AbstractArray)
    ld = abs(b .- a)
    res = zeros(size(a,2))
    for j in 1:size(a,2)
        d = 0.
        for i in 1:size(a,1)
            @inbounds c = (ld[i,j] > .5) ? 1-ld[i,j] : ld[i,j]
            d += c*c
        end
        res[j] = sqrt(d)
    end
    res
end

periodic_euclidean(a::AbstractArray, b::AbstractArray, D::Real) = D*periodic_euclidean(a,b)

#TODO ranks, parents. I want this to work with the above api.
#type ID_Set#datatype that is similar to IntDisjointSets but more general
    #doesn't require sequential integers
#    main_dict::Dictionary
#    ngroups::Int

#    function ID_Set(I::Array{Int, 1})
#        main_dict = Dict{Int,Set}(i => Set(i) for i in I)
#        ngroups = size(I,1)
#        new(main_dict, ngroups)
#    end
#end

#function union!(ds::ID_Set, i1::Int, i2::Int)
#    if i2 in ds.main_dict[i1]
#        return nothing
#    end

#    push!(ds.main_dict[i1], i2)
#    ds.main_dict[i2] = ds.main_dict[i1]
#    ngroups-=1
#    nothing
#end

#function num_groups(ds::ID_Set) = ID_Set.ngroups
