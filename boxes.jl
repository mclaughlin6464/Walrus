#@Author Sean McLaughlin 3/2016
#this file contains definitiosn relating to the box object
#box object allows for partitioning of the particle space

type Box
    dim::Array{Tuple{Float64, Float64}, 1}#bounadries in dims 1:3
    NPart::Int
    orig_idxs::Array{Int, 1}#original indexs of the passed in particles
    x::Array{Float64, 2} # positions
    #v::Array{Float64, 2} # velocities
    tree::KDTree
    ds::IntDisjointSets
    #Boxes to the North, eastst and up.
    #Will be set after constructing.
    N_Box::Box
    U_Box::Box
    E_Box::Box

    function Box(dim::Array{Tuple{Float64, Float64}, 1}, orig_idxs::Array{Int, 1}, x::Array{Float64, 2})
        NPart = size(orig_idxs, 1)
        b = new(dim,NPart, orig_idxs, x, KDTree(x), IntDisjointSets(NPart))
        b.N_Box = b
        b.U_Box = b
        b.E_Box = b
    end
end

function find_groups!(b::Box, l::Real)

    for i in 1:b.NPart
        idxs = IntSet(inrange(b.tree, b.x[:,i], l, false)) # within search radius

        for idx in idxs #
            union!(b.ds,i,idx)
        end
        if (num_groups(b.ds) == 1) # just in case people use too large a linking length don't waste time
            #println("FOF: All points were linked. Exiting." )
            break
        end  # in case everything has been joined already exit
    end
end

function find_groups(k::Tuple, b::Box, l::Real)

    find_groups!(b, l)
    return k, b
end

function make_boxes(BoxSize::Real, NBoxes::Int, x::Array{Float64,2})
    NDim = size(x,1)
    NBoxesDim = NBoxes
    if NDim ==2
        NBoxesDim = convert(Int, floor(sqrt(NBoxes) ))
    elseif NDim == 3
        NBoxesDim = convert(Int, floor(cbrt(NBoxes)))
    elseif NDim != 1
        throw(ArgumentError("Dimensions other than 1:3 unsupported"))
    end

    SubBoxSize = BoxSize/NBoxesDim
    IdxDict = Dict{Tuple, Array{Int,1}}()
    BoxDict = Dict{Tuple, Box}()
    #TODO make this more efficient with less copy pasting
    #NOTE 1 based indexing. #Julia
    if NDim == 1
        for i in 1:size(x,2)
            xBoxIdx = convert(Int, floor(x[1, i]/SubBoxSize))+1
            key = (xBoxIdx,)
            if !haskey(IdxDict, key)
                IdxDict[key] = Int[]
            end
            push!(IdxDict[key], i)
        end
    elseif NDim == 2
        for i in 1:size(x,2)
            xBoxIdx = convert(Int, floor(x[1, i]/SubBoxSize))+1
            yBoxIdx = convert(Int, floor(x[2, i]/SubBoxSize))+1
            key = (xBoxIdx,yBoxIdx)
            if !haskey(IdxDict, key)
                IdxDict[key] = Int[]
            end
            push!(IdxDict[key], i)
        end
    else
        for i in 1:size(x,2)
            #NOTE could combine this into a comprehension; why bother tho?
            xBoxIdx = convert(Int, floor(x[1, i]/SubBoxSize))+1
            yBoxIdx = convert(Int, floor(x[2, i]/SubBoxSize))+1
            zBoxIdx = convert(Int, floor(x[3, i]/SubBoxSize))+1
            key = (xBoxIdx, yBoxIdx, zBoxIdx)
            if !haskey(IdxDict, key)
                IdxDict[key] = Int[]
            end
            push!(IdxDict[key], i)
        end
    end

    #problems if there aren't particles in each cell.
    if NDim == 1
        for i in 1:NBoxesDim
            dim = [((i-1)*SubBoxSize, i*SubBoxSize)]

            BoxDict[(i,)] = Box(dim, IdxDict[(i,)], x[:, IdxDict[(i,)]])
        end
    elseif NDim == 2
        for i in 1:NBoxesDim, j in 1:NBoxesDim
            dim = [((i-1)*SubBoxSize, i*SubBoxSize),((j-1)*SubBoxSize, j*SubBoxSize)]

            BoxDict[(i,j)] = Box(dim, IdxDict[(i,j)], x[:, IdxDict[(i,j)]])
        end
    else
        for i in 1:NBoxesDim, j in 1:NBoxesDim, k in 1:NBoxesDim
            dim = [((i-1)*SubBoxSize, i*SubBoxSize),((j-1)*SubBoxSize, j*SubBoxSize),((k-1)*SubBoxSize, k*SubBoxSize)]

            BoxDict[(i,j,k)] = Box(dim, IdxDict[(i,j,k)], x[:, IdxDict[(i,j,k)]])
        end
    end
    #link up neighbors
    for (tup, b) in BoxDict
        if NDim == 1
            if NBoxesDim == 1
                key = (1,)
            else
                key = tup[1] == NBoxesDim-1 ? (NBoxesDim, ) : (rem((tup[1]+1),NBoxesDim), )
            end
            b.N_Box = BoxDict[key]
        elseif NDim == 2
            if NBoxesDim == 1
                key = (1,1)
                b.N_Box = BoxDict[key]
                b.E_Box = BoxDict[key]

            else
                key = tup[1] == NBoxesDim-1 ? (NBoxesDim, tup[2]) : (rem((tup[1]+1),NBoxesDim),tup[2] )
                b.N_Box = BoxDict[key]

                key = tup[2] == NBoxesDim-1 ? (tup[1], NBoxesDim) : (tup[1], rem((tup[2]+1),NBoxesDim))
                b.E_Box = BoxDict[key]
            end
        else
            if NBoxesDim == 1
                key = (1,1,1)
                b.N_Box = BoxDict[key]
                b.E_Box = BoxDict[key]
                b.U_Box = BoxDict[key]
            else
                key = tup[1] == NBoxesDim-1 ? (NBoxesDim, tup[2], tup[3]) : (rem((tup[1]+1), NBoxesDim),tup[2], tup[3] )
                b.N_Box = BoxDict[key]

                key = tup[2] == NBoxesDim-1 ? (tup[1], NBoxesDim, tup[3]) : (tup[1], rem((tup[2]+1),NBoxesDim), tup[3])
                b.E_Box = BoxDict[key]

                key = tup[3] == NBoxesDim-1 ? (tup[1],tup[2], NBoxesDim) : (tup[1],tup[2], rem((tup[3]+1),NBoxesDim))
                b.U_Box = BoxDict[key]
            end
        end

    end
    return BoxDict
end

#NPart, BoxSize from where?
function link_boundaries(BoxDict::Dict,BoxSize::Real, l::Real, NPart::Int )#Not sure about NPart
    gds = IntDisjointSets(NPart)# a global disjoint set
    gds.ngroups = 0
    #make the global disjoint set.
    for (key, box) in BoxDict
        gds.parents[box.orig_idxs] = box.orig_idxs[box.ds.parents]
        gds.ranks[box.orig_idxs] = box.ds.ranks
        gds.ngroups+=box.ds.ngroups
    end

    for (key, box) in BoxDict
        link_boundaries!(gds, box,BoxSize, l)
    end
    return gds
end

#TODO get global boxsize here?
function link_boundaries!(gds::IntDisjointSets, b::Box, BoxSize::Real, l::Real)
    NDims = size(b.dim,1)
    other_boxes = [b.N_Box]
    if NDims > 1
        push!(other_boxes, b.E_Box)
    end
    if NDims>2
        push!(other_boxes, b.U_Box)
    end

    for (dim, ob) in zip(1:3, other_boxes)
        BoxCorrection = zeros(size(ob.x,1))
        if b.dim[dim][2] != ob.dim[dim][1] #we need to do periodic boundaries!
            BoxCorrection[dim] = BoxSize
        end
        for i in 1:ob.NPart
            if ob.x[dim, i] > ob.dim[dim][1]+l#we know it's out of range
                continue
            end

            idxs = IntSet(inrange(b.tree, ob.x[:,i]+BoxCorrection, l, false)) # within search radius

            for idx in idxs #
                union!(gds,ob.orig_idxs[i],b.orig_idxs[idx])
            end
        end
    end
end

function link_boundaries(BoxDict::Dict,BoxSize::Real, l::Real, NPart::Int )#Not sure about NPart
    gds = IntDisjointSets(NPart)# a global disjoint set
    gds.ngroups = 0
    #make the global disjoint set.
    for (key, box) in BoxDict
        gds.parents[box.orig_idxs] = box.orig_idxs[box.ds.parents]
        gds.ranks[box.orig_idxs] = box.ds.ranks
        gds.ngroups+=box.ds.ngroups
    end

    for (key, box) in BoxDict
        link_boundaries!(gds, box,BoxSize, l)
    end
    return gds
end

function get_halo_idxs(ds::IntDisjointSets, minlength::Int)#return the idxs that define halos
    #make halo idxs from the sets
    groupDict = Dict{Int, Array{Int,1}}()
    Npart = size(ds.parents, 1)
    #TODO slow, not sure how to cleverly parallelize
    for i in 1:Npart
        p = find_root(ds,i)
        if !haskey(groupDict, p)
            groupDict[p] = Int[]
        end
        push!(groupDict[p], i)
    end

    gps = sort(collect(values(groupDict)), by = x-> length(x), rev = true)
    Ngroups = -1
    for i in eachindex(gps)
        if length(gps[i])< minlength
            Ngroups = i
            break
        end
    end

    println("FOF: Found ", Ngroups, " with ", minlength, " or more points")

    return gps[1:Ngroups] # An array of different  sized arrays containing the particle ids in the groups is returned
end
