#handles I/O of gadget files
#Written by Tom Abel 1/2016
#modified by Sean McLaughlin 3/2016

type gadget_header
    npart::Array{Int32,1}
    mass::Array{Float64,1}
    time::Float64
    redshift::Float64
    flag_sfr::Int32
    flag_feedback::Int32
    npartTotal::Array{UInt32,1}
    flag_cooling::Int32
    num_files::Int32
    BoxSize::Float64
    Omega0::Float64
    OmegaLambda::Float64
    HubbleParam::Float64
    flag_stellarage::Int32
    flag_metals::Int32
    npartTotalHighWord::Array{Int32,1}
    flag_entropy_instead_u::Int32
    fill::AbstractString
end

function read_gadget_header(filename, verbose = false)
    blk = UInt32(0)
    head = gadget_header(zeros(Int32,6), zeros(Float64,6), Float64(0.), Float64(0.),
                         Int32(0), Int32(0), zeros(UInt32,6), Int32(0), Int32(0), Float64(0.),
                         Float64(0.), Float64(0.), Float64(0.), Int32(0), Int32(0), zeros(Int32,6), Int32(0), lpad(' ', 60, ' '))


    #    fn = "/Users/tabel/Research/TetraPart/data/23eV/Gadget/NP128/snapshot_039"
    istream = open(filename, "r")
    if verbose
        println("am at position ", position(istream), " in this file")
    end

    blk = read(istream, Int32, (1))
    if verbose
        println(blk)
    end

    head.npart = read(istream, Int32, (6))
    head.mass = read(istream, Float64, (6))
    head.time = read(istream, Float64)
    head.redshift = read(istream, Float64)
    head.flag_sfr = read(istream, Int32)
    head.flag_feedback = read(istream, Int32)
    head.npartTotal = read(istream, UInt32, (6))
    head.flag_cooling = read(istream, Int32)
    head.num_files = read(istream, Int32)
    head.BoxSize = read(istream, Float64)
    head.Omega0 = read(istream, Float64)
    head.OmegaLambda = read(istream, Float64)
    head.HubbleParam = read(istream, Float64)
    head.flag_stellarage = read(istream, Int32)
    head.flag_metals = read(istream, Int32)
    head.npartTotalHighWord = read(istream, UInt32, (6))
    head.flag_entropy_instead_u = read(istream, Int32)
    dummy = read(istream, Char,60)
    blk = read(istream, Int32)
    if verbose
        println(blk)

        println("read header from file:", filename)
        println("redshift:", head.redshift, "\t scale factor:", head.time)
        println("npart: " , head.npart)
        println("masses:", head.mass)
        println("am at position ", position(istream), " in this file")
    end

    close(istream)
    return head
end

#TODO Quiet option?
function read_gadget_data(filename::AbstractString, verbose = false)

    head = read_gadget_header(filename, verbose)

    istream = open(filename, "r")
    seek(istream, 264)# jump past header
    if verbose
        println("am at position ", position(istream), " in this file")
    end

    Npart = sum(head.npart)

    blk = read(istream, Int32)
    if verbose
        println(blk)
    end
    #p = particle(zeros(Float32, Npart, 3), zeros(Float32, Npart, 3), zeros(Float32, Npart))
    p = Particles(zeros(Float32, 1, 3), zeros(Float32, 1, 3), zeros(Int32, 1), zeros(Float32, 1), zeros(Float32, 1))
    p.x = read(istream, Float32, (3, Npart))
    blk = read(istream, Int32)
    if verbose
        println(blk)
    end

    blk = read(istream, Int32)
    if verbose
        println(blk)
    end
    p.v = read(istream, Float32, (3, Npart))  # read velocities
    blk = read(istream, Int32)
    if verbose
        println(blk)
    end

    blk = read(istream, Int32)
    if verbose
        println(blk)
    end
    p.id = read(istream, Int32, (Npart))        # read ids
    blk = read(istream, Int32)
    if verbose
        println(blk)
    end

    NwithMass = 0
    for i in 1:6
       if head.mass[i] > Float32(0.) && head.npart[i] > 0
           NwithMass += head.npart[i]
       end
    end
    if verbose
        println("# of particles with masses stored in file: ", NwithMass)
    end

    if NwithMass > 0
       blk = read(istream, Int32)
       if verbose
           println(blk)
       end
       masses = read(istream, Float32, (NwithMass)) # read masses
       blk = read(istream, Int32)
       if verbose
           println(blk)
       end
    end
    # define the mass array even in the case where we have all equal masses.
    # could do this better  ....
    p.m = zeros(Float32,Npart)
    count = 1
    wcoutn = 1
    for i in 1:6
       if head.npart[i] > 0
           ncount = (count+head.npart[i])
           if head.mass[i] > Float32(0.)
               if verbose
                   println(i," ", count, " ",ncount)
               end
               p.m[count:(ncount-1)] = head.mass[i]
               count = ncount
           elseif head.npart[i] > 0
               nwcount = wcount + head.npart[i]
               p.m[count:ncount-1] = masses[wcount:nwcount-1]
               count  = ncount
               wcount = nwcount
           end
       end
    end

    blk = read(istream, Int32)
    if verbose
        println(blk)
    end
    try
        p.pot = read(istream, Float32, (Npart)) # read potential
        blk = read(istream, Int32)
        if verbose
            println(blk)
        end
    catch LoadError#no potentials
        p.pot = zeros(Float32, Npart)
    end

    close(istream)
    #TODO doesn't return redshift and other relevant info
    return p, head
end

function read_gadget_data(filenames, verbose = false)
    p, head = read_gadget_data(filenames[1], verbose)
    for fname in filenames[2:end]
        newp, newhead = read_gadget_data(fname, verbose)
        join!(p, newp)
    end
    return p, head
end
