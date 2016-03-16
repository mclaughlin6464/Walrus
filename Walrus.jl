#Main module for the Walrus FOF clustering code
#Author: Sean McLaughlin, March 2016

module Walrus

using ArgParse
using NearestNeighbors
using DataStructures
using Distances
using Glob

include("gadget_load.jl")
include("FOF.jl")
include("boxes.jl")
include("halos.jl")

const G = 6.67E-11

#things to add:
#linking length
#minlength
function read_filenames()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--one"
            help = "Run on the file passed in, don't treat as root."
            nargs = 0
        "gadget_fname_root"
            help = "the root filename of the gadget output file(s)"
            required = true

        "output_fname"
            help = "the filename to write the output catalog to"
            required = true
    end

    d = parse_args(s)
    return d["one"], d["gadget_fname_root"], d["output_fname"]
end

use_file, input_fname, output_fname = read_filenames()

#use the file as passed in
if use_file
    particles, header = read_gadget_data(input_fname, false)
else #use the file as a root to load all files
    dir_idx = rsearchindex(input_fname, "/")#split into pattern and dir
    fnames = glob(input_fname[dir_idx+1:end]* "*", input_fname[1:dir_idx])
    particles, header = read_gadget_data(fnames, false)
end

const BoxSize = header.BoxSize
const H = header.HubbleParam

#masses are in units of 1e10 m_sun/h, and distances are Mpc/h
#h=0.633657

Npart = size(particles, 2)
particles = particles[1:10:Npart]
Npart = size(particles, 2)
Ndim = size(particles, 1)
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

#gps = groups(particles.x, dx, 10, particles.v, 1000)
BoxDict = make_boxes(BoxSize, 8, particles.x)
@time( begin
for box in values(BoxDict)
    find_groups!(box, dx)
end
end )

@time gds = link_boundaries(BoxDict, BoxSize, dx, Npart)

@time gps = get_halo_idxs(gds, 10)

if size(gps,1) == 0
    println("No halos found; exiting.")
    quit()
end

halos = Array{Halo}(size(gps,1))
halo_ids = collect(1:size(gps,1))

for (halo_id, halo_parts) in zip(halo_ids, gps)
    halos[halo_id] = Halo(halo_id, particles[halo_parts])
end

sort!(halos, by = x-> x.com[1])#sort by x value

#TODO split up the vectors into their dims
open(output_fname, "w") do io
    write(io, "id, N, total_mass, com, R_200, velocity, vel_disp, J, E, spin_param\n")
    for h in halos
        write(io, halo_output(h) )
        write(io, "\n")
    end
end

end
