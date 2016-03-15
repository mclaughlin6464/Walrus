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
        NPart = size(orig_idxs, 2)
        b = new(dim,NPart, orig_idxs, x, KDTree(x), IntDisjointSets(NPart))
        b.N_Box = b
        b.U_Box = b
        b.E_Box = b
    end
end

function find_groups!(b::Box, l::Real)
    for i in 1:b.NPart
        idxs = IntSet(inrange(b.tree, b.x[:,i], l, false)) # within search radius

        println(b.ds)
        println(idxs)
        println(i)
        for idx in idxs #
            union!(b.ds,i,idx)
        end
        if (num_groups(b.ds) == 1) # just in case people use too large a linking length don't waste time
            #println("FOF: All points were linked. Exiting." )
            break
        end  # in case everything has been joined already exit
    end
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
            key = tup[1] == NBoxesDim-1 ? (NBoxesDim, ) : (rem((tup[1]+1),NBoxesDim), )
            b.N_Box = BoxDict[key]
        elseif NDim == 2
            key = tup[1] == NBoxesDim-1 ? (NBoxesDim, tup[2]) : (rem((tup[1]+1),NBoxesDim),tup[2] )
            b.N_Box = BoxDict[key]

            key = tup[2] == NBoxesDim-1 ? (tup[1], NBoxesDim) : (tup[1], rem((tup[2]+1),NBoxesDim))
            b.E_Box = BoxDict[key]
        else
            key = tup[1] == NBoxesDim-1 ? (NBoxesDim, tup[2], tup[3]) : (rem((tup[1]+1), NBoxesDim),tup[2], tup[3] )
            b.N_Box = BoxDict[key]

            key = tup[2] == NBoxesDim-1 ? (tup[1], NBoxesDim, tup[3]) : (tup[1], rem((tup[2]+1),NBoxesDim), tup[3])
            b.E_Box = BoxDict[key]

            key = tup[3] == NBoxesDim-1 ? (tup[1],tup[2], NBoxesDim) : (tup[1],tup[2], rem((tup[3]+1),NBoxesDim))
            b.U_Box = BoxDict[key]
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
        link_boundaries!(gds, box, l)
    end
    return gds
end

#TODO get global boxsize here?
function link_boundaries!(gds::IntDisjointSets, b::Box, BoxSize::Real, l::Real)
    other_boxes = [b.N_Box]
    if ! is(b, b.E_Box)
        push!(other_boxes, b.E_box)
    end
    if ! is(b, b.U_Box)
        push!(other_boxes, b.U_Box)
    end

    for (dim, ob) in zip(1:3, other_boxes)
        BoxCorrection = zero(size(ob.x,1))
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

function get_halo_idxs(ds::IntDisjointSets, minlength::Int)#return the idxs that define halos
    #make halo idxs from the sets
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
    #println("FOF: Found ", Ngroups, " with ", minlength, " or more points")
    return gps # An array of different  sized arrays containing the particle ids in the groups is returned
end
