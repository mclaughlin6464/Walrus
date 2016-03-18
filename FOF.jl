# Friends of Friends clustering.  (Tom Abel, 1/2016)

"""
Friends of Friends Algorithm: for sets of points x, using linking length l and returning groups with more than minlength members in an array of arrays of indeces inot the point data provided.
    using Gadfly
    x = rand((2,100))
    gps = groups(x,.1, 3)
    n = 1
    plot(layer(x=x[1,gps[n]], y=x[2,gps[n]], Geom.point,Theme(default_color=colorant"red")), layer(x=x[1,:], y=x[2,:], Geom.point, Theme(default_color=colorant"blue"))  )

"""
#TODO verbose option
function groups(x, l, minlength, v = nothing, l_v = nothing)
    #make sure the 6D option is propertly specified
    if (v == nothing && l_v != nothing) || (v != nothing && l_v == nothing)
        throw(ArgumentError("Both v and l_v must be specified"))
    end

    Npart = size(x,2)
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
        #6D FOF
        if v != nothing
            #if l_v large, this line is very slow!
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

<<<<<<< HEAD
    #make halo idxs from the sets

    groupDict = Dict{Int, Array{Int,1}}()
=======
    groupDict = Dict{Int, Array{Int,1}}()

>>>>>>> pFOF
    for i in 1:Npart
        p = find_root(ds,i)
        if !haskey(groupDict, p)
            groupDict[p] = Int[]
        end
        push!(groupDict[p], i)
    end
    
    gps = sort(collect(values(groupDict)), by = x-> length(x), rev = true)
    Ngroups = -1

<<<<<<< HEAD
=======
    gps = sort(collect(values(groupDict)), by = x-> length(x), rev = true)
    Ngroups = -1
>>>>>>> pFOF
    for i in eachindex(gps)
        if length(gps[i])< minlength
            Ngroups = i
            break
        end
    end

    println("FOF: Found ", Ngroups, " with ", minlength, " or more points")

    return gps[1:Ngroups]# An array of different  sized arrays containing the particle ids in the groups is returned
end

function link_boundaries(tree::KDTree, dim::Int, x::Array{Float64,1}, l::Real)
    #link across periodic boundary conditions in one dim
    if x[dim] > l
        return IntSet()
    end

    BoxCorrection = zeros(size(x,1))
    BoxCorrection[dim] = BoxSize
    idxs = IntSet(inrange(tree, x+BoxCorrection, l, false)) # within search radius
    return idxs
end

#link periodic boundary conditions across all 3 dimensions
function link_boundaries(tree::KDTree, x::Array{Float64,1}, l::Real)
    Ndim = size(x,1)
    output_set = IntSet()

    for dim in 1:Ndim
        idxs = link_boundaries(tree, dim, x, l)
        union!(output_set, idxs)
    end
    return output_set
end
