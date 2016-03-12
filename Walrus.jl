#Main module for the Walrus FOF clustering code
#Author: Sean McLaughlin, March 2016

module Walrus

using ArgParse
using NearestNeighbors
using DataStructures
using Distances

include("gadget_load.jl")
include("FOF.jl")
include("halos.jl")

const G = 6.67E-11

#things to add:
#linking length
#minlength
function read_filenames()
    s = ArgParseSettings()

    @add_arg_table s begin
        "gadget_fname"
            help = "the filename of the gadget output file"
            required = true

        "output_fname"
            help = "the filename to write the output catalog to"
            required = true
    end

    d = parse_args(s)
    return d["gadget_fname"], d["output_fname"]
end

input_fname, output_fname = read_filenames()

particles, header = read_gadget_data(input_fname, false)

const BoxSize = header.BoxSize
const H = header.HubbleParam

#masses are in units of 1e10 m_sun/h, and distances are Mpc/h
#h=0.633657

Npart = size(particles, 1)
Ndim = size(particles, 2)
Npart_1D = Npart
if Ndim == 2
    Npart_1D = sqrt(Npart)
elseif Ndim == 3
    Npart_1D = cbrt(Npart)
end

q = .5
dx = q/Npart_1D #heard this is a good guess for linking length; not convinced.
dx = 2

#minlength = floor(Npart/1000) #?

gps = groups(particles.x, dx, 10, particles.v, 1000)

if size(gps,1) == 0
    println("No halos found; exiting.")
    quit()
end

halos = Array{Halo}(size(gps,1))
halo_ids = collect(1:size(gps,1))
h = Halo(1, particles[gps[1]])

for tup in zip(halo_ids, gps)
    halo_id = tup[1]
    halo_parts = tup[2]
    halos[halo_id] = Halo(halo_id, particles[halo_parts])
end

#TODO split up the vectors into their dims
open(output_fname, "w") do io
    write(io, "id, N, total_mass, com, R_200, velocity, vel_disp, J, E, spin_param\n")
    for h in halos
        write(io, halo_output(h) )
        write(io, "\n")
    end
end

end
