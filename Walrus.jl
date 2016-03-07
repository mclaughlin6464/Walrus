#Main module for the Walrus FOF clustering code
#Author: Sean McLaughlin, March 2016

module Walrus

using ArgParse
using gadget_load
using FOF
using Distances

#things to add:
#linking length
#minlength

const G = 6.67E-11

type Halo
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

    function Halo(id::Int, particles::particle, halo_parts::Array{Int, 1})
        total_mass = sum(particles.m[halo_parts])
        com = sum(particles.x[halo_parts].*particles.m[halo_parts], 1)/total_mass
        R_200_h= cbrt(G*total_mass)/100#R200 times h^2/3. Not sure how to simplify more
        N = size(halo_parts, 1)
        velocity = mean(particles.v[halo_parts], 1)
        vel_disp = std(particles.v[halo_parts], 1)
        for i in halo_parts
            println("$(particles.x[i])", "$(com) " )
        end
        #println("$(particles.m[halo_parts].*particles.x[halo_parts] )" )
        #println("$(particles.m[halo_parts].*(particles.x[halo_parts]-com) )" )
        J = sum(cross(particles.m[halo_parts].*(particles.x[halo_parts]-com),
            (particles.v[halo_parts]-velocity)),1 )
        E = abs(sum(particles.pot[halo_parts],1)) #TODO pot of halo or all parts?
        if E == 0 #gadget didn't calculate the potentials!
            #TODO better flag than this!
            calc_potential!(particles, halo_parts)
            E = abs(sum(particles.pot[halo_parts],1))
        end
        spin_param = J*sqrt(E)/(G*sqrt(total_mass)^5)

        new(id, N, total_mass, com, R_200_h, velocity, vel_disp, J, E, spin_param)

    end
end

function halo_output(halo::Halo) #ret"urns csv row for a given halo
    output = Array{AbstractString}(10)
    output[1] = "$(halo.id)"
    output[2] = "$(halo.N)"
    output[3] = "$(halo.total_mass)"
    output[4] = "$(halo.com)"
    output[5] = "$(halo.R_200_h)"
    output[6] = "$(halo.velocity)"
    output[7] = "$(halo.vel_disp)"
    output[8] = "$(halo.J)"
    output[9] = "$(halo.E)"
    output[10] = "$(halo.spin_param)"

    return join(output, ", ")
end

#TODO doesn't do periodic boundaries!
function calc_potential!(particles::particle)
    #sometiems gadget doesn't do the work for us
    Npart = size(particles.x,2)
    for i in 1:Npart
        pot_holder = 0
        for j in 1:Npart
            if i == j
                continue
            end
            r = euclidean(particles.x[i], particles.x[j])
            pot_holder+=- G*particles.m[i]*particles.m[j]/r
        end
        particles.pot[i] = pot_holder
    end
    nothing
end

function calc_potential!(particles::particle,halo_parts::Array{Int} )
    #sometiems gadget doesn't do the work for us
    #this method takes a subset of particles, not all of them.
    Npart = size(particles.x,2)
    halo_part_set = IntSet(halo_parts)
    for i in 1:Npart
        if !(i in halo_part_set)
            continue
        end
        pot_holder = 0
        for j in 1:Npart
            if i == j || !(j in halo_part_set)
                continue
            end
            r = euclidean(particles.x[i], particles.x[j])
            pot_holder+=- G*particles.m[i]*particles.m[j]/r
        end
        particles.pot[i] = pot_holder
    end
    nothing
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

particles = read_gadget_data(input_fname, true)

Npart= size(particles.x,2)
Ndim = size(particles.x,1)
Npart_1D = Npart
if Ndim == 2
    Npart_1D = sqrt(Npart)
elseif Ndim == 3
    Npart_1D = cbrt(Npart)
end

q = .5
dx = q/Npart_1D #heard this is a good guess for linking length; not convinced.
dx = 10
println("$(dx)")
minlength = floor(Npart/1000) #?

gps = groups(particles.x, dx, minlength)

halos = Array{Halo}(size(gps))
halo_ids = collect(1:size(gps,1))

h = Halo(1, particles, gps[1])

for tup in zip(halo_ids, gps)
    halo_id = tup[1]
    halo_parts = tup[2]
    halos[halo_id] = Halo(halo_id, particles, halo_parts)
end

open(output_fname, "w") do io
    write(io, "id, N, total_mass, com, R_200_h, velocity, vel_disp, J, E, spin_param\n")
    for h in halos
        write(io, halo_output(h) )
        write(io, "\n")
    end
end

end
