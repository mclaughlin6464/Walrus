#Main module for the Walrus FOF clustering code
#Author: Sean McLaughlin, March 2016

module Walrus

using ArgParse
using gadget_load
using FOF

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

end
