#Main module for the Walrus FOF clustering code
#Author: Sean McLaughlin, March 2016

module Walrus

using ArgParse
using gadget_load
using FOF

#things to add:
#linking length
#minlength

const G = 6.67E-11

type halo
    id::Int
    N::Int
    total_mass::Float64
    com::Array{Float64,2}
    R_200_h::Float64
    velocity::Float64
    vel_disp::Float64
    J::Float64
    E::Float64
    spin_param::Float64

    function halo(id::Int, particles::Array{particle}, halo_parts::Array{Int})
        total_mass = sum(particles.m[halo])
        com = sum(particles.x[halo].*particles.m[halo], 2)/total_mass
        R_200_h= cbrt(G*total_mass)/100#R200 times h^2/3. Not sure how to simplify more
        N = size(halo)
        velocity = mean(particles.v[halo], 2)
        vel_disp = std(particles.v[halo], 2)
        J = sum(cross(particles.m[halo]*(particles.x[halo]-com), (particles.v[halo]-velocity)) )
        E = sum(particles.pot[halo],2) #TODO check sign
        spin_param = J*sqrt(E)/(G*sqrt(total_mass)^5)

        new(id, N, total_mass, com, R_200_h, N, velocity, vel_disp, J, E, spin_param)

    end
end


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

particles = read_gadget_data(input_fname)

Npart= size(particles.x,2)
Ndim = size(particles.x,1)
Npart_1D = Npart
if Ndim == 2
    Npart_1D = sqrt(Npart)
elseif Ndim == 3
    Npart_1D = cbrt(Npart)
q = .05
dx = q/Npart_1D #heard this is a good guess for linking length; not convinced.

minlength = Npart/100 #?

gps = groups(particles.x, dx, minlength)

halos = Array{halo}(size(gps))
halo_ids = collect(1:size(gps))
for halo_id, halo_parts in zip(halo_ids, gps):
    halos[halo_id] = halo(halo_id, particles, halo_parts)

end
