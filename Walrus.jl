#Main module for the Walrus FOF clustering code
#Author: Sean McLaughlin, March 2016

module Walrus

export get_halos,get_halos_p, halo_output
export read_gadget_data
export Particles, Halo, Box

#may have to be changed to require's
using NearestNeighbors
using DataStructures
using Distances

include("gadget_load.jl")
#include("FOF.jl")
include("boxes.jl")
include("halos.jl")

const G = 6.67E-11

function get_halos(particles::Particles, BoxSize::Real,H::Real, dx=2, NBoxes = 8)
    Npart = size(particles, 2)

    BoxDict = make_boxes(BoxSize, NBoxes, particles.x)
    #Use pmap here
    for box in values(BoxDict)
        find_groups!(box, dx)
    end

    gds = link_boundaries(BoxDict, BoxSize, dx, Npart)

    gps = get_halo_idxs(gds, 10)

    if size(gps,1) == 0
        println("No halos found; exiting.")
        return Halo(gps[1])
    end

    halos = Array{Halo}(size(gps,1))
    halo_ids = collect(1:size(gps,1))

    for (halo_id, halo_parts) in zip(halo_ids, gps)
        halos[halo_id] = Halo(halo_id, particles[halo_parts], H, BoxSize)
    end
    return halos
end

function get_halos_p(particles::Particles, BoxSize::Real,H::Real, dx=2, NBoxes = 8)
    Npart = size(particles, 2)

    BoxDict = make_boxes(BoxSize, NBoxes, particles.x)
    #Use pmap here
    #for box in values(BoxDict)
    #    find_groups!(box, dx)
    #end

    boxes = pmap(d -> find_groups(d[1], d[2], dx), BoxDict)

    for (k,b) in boxes
        BoxDict[k] = b
    end

    gds = link_boundaries(BoxDict, BoxSize, dx, Npart)

    gps = get_halo_idxs(gds, 10)

    if size(gps,1) == 0
        println("No halos found; exiting.")
        return Halo(gps[1])
    end

    halos = Array{Halo}(size(gps,1))
    halo_ids = collect(1:size(gps,1))

    for (halo_id, halo_parts) in zip(halo_ids, gps)
        halos[halo_id] = Halo(halo_id, particles[halo_parts], H, BoxSize)
    end
    return halos
end

end
