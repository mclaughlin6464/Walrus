#This module defines the objects used in the halo finder
#Author: Sean McLaughlin

#Particles

type Particles
    x::Array{Float64} # positions
    v::Array{Float64} # velocities
    id::Array{Int,1}  # ids
    m::Array{Float64,1} # mass
    pot::Array{Float64,1} # gravitational potential

    function Particles(x, v, id, m , pot)
        #assert dimensionality
        @assert size(x,1) == size(v,1)
        @assert size(x,2) == size(v,2)
        @assert size(x,1) == size(id,1) || size(x,2) === size(id, 1)
        @assert size(x,1) == size(m,1)|| size(x,2) === size(m, 1)
        @assert size(x,1) == size(pot, 1)|| size(x,2) === size(pot, 1)

        new(x,v,id,m,pot)
    end
end
Base.show(P::Particles) = "$(size(P.id)) Particles"
Base.size(P::Particles) = size(P.id)
Base.size(P::Particles, i::Int) = size(P.x, i)
Base.endof(P::Particles) = endof(P.id)

#Base.linearindexing(::Particles) = Base.LinearFast()
Base.getindex(P::Particles, i::Int) = Particles(Array(P.x[:, i] ), Array(P.v[:, i]),
                                        [P.id[i]], [P.m[i]], [P.pot[i]] )
Base.getindex(P::Particles, I) = Particles(P.x[:, I], P.v[:, I],
                                                    P.id[I], P.m[I], P.pot[I])
Base.start(P::Particles) = P[1]
Base.next(P::Particles, state) = P[state+1]
Base.done(P::Particles, state) = state == size(P,1)

<<<<<<< HEAD

=======
>>>>>>> pFOF
function join!(P::Particles, P2::Particles)
    P.x = cat(2, P.x, P2.x)
    P.v = cat(2,P.v, P2.v)
    P.id = cat(1,P.id, P2.id)
    P.m= cat(1, P.m, P2.m)
    P.pot = cat(1, P.pot, P2.pot)
end
#didn't define setIndex!
type Halo
    id::Int
    N::Int
    total_mass::Float64
    com::Array{Float64,1}
    R_200::Float64
    velocity::Array{Float64,1}
    vel_disp::Array{Float64,1}
    J::Float64
    E::Float64
    spin_param::Float64

    #TODO units?
    function Halo(id::Int, P::Particles, H::Real, BoxSize::Real)

        total_mass = sum(P.m)
        com = reshape(transpose(sum(P.x.*transpose(P.m), 2)./total_mass), 3)
        R_200= cbrt(G*total_mass)/(100*cbrt(H^2))
        N = size(P, 2)
        velocity = reshape(transpose(mean(P.v, 2)), 3)
        vel_disp = reshape(transpose(std(P.v, 2)), 3)
        cproducts = [cross(P.m[i]*(P.x[:, i]-com),(P.v[:,i]-velocity) ) for i in 1:N ]
        J = sum( [sqrt(dot(cproducts[i],cproducts[i])) for i in 1:N])
        E = abs(sum(P.pot)) #TODO pot of halo or all parts?
        if E == 0 #gadget didn't calculate the potentials!
            #TODO better flag than this!
            calc_potential!(P, BoxSize)
            E = abs(sum(P.pot))
        end
        spin_param = J*sqrt(E)/(G*sqrt(total_mass)^5)

        new(id, N, total_mass, com, R_200, velocity, vel_disp, J, E, spin_param)

    end
end

function halo_output(halo::Halo) #ret"urns csv row for a given halo
    output = Array{AbstractString}(10)
    output[1] = "$(halo.id)"
    output[2] = "$(halo.N)"
    output[3] = "$(halo.total_mass)"
    output[4] = "$(halo.com)"
    output[5] = "$(halo.R_200)"
    output[6] = "$(halo.velocity)"
    output[7] = "$(halo.vel_disp)"
    output[8] = "$(halo.J)"
    output[9] = "$(halo.E)"
    output[10] = "$(halo.spin_param)"

    return join(output, ", ")
end

function calc_potential!(P::Particles, BoxSize::Real)
    #sometiems gadget doesn't do the work for us
    Npart = size(P.x,2)
    for i in 1:Npart
        pot_holder = 0
        for j in 1:Npart
            if i == j
                continue
            end
            r = periodic_euclidean(P.x[:, i], P.x[:, j], BoxSize)
            pot_holder+=- G*P.m[i]*P.m[j]/r[1]
        end
        P.pot[i] = pot_holder
    end
    nothing
end

#TODO is wrong!
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
