# exports
export readdata, writedata

using DelimitedFiles

function readdata(filename)
    data = readdlm(filename, skipstart=2)
    x = 1*transpose(data[:,1:3])
    v = 1*transpose(data[:,4:6])
    a_ids = 1*transpose(data[:,7])
    m_ids = 1*transpose(data[:,8])
    (x, v, a_ids, m_ids)
end


function writedata(filename, vals)
    N = size(vals[1])[1]
    open(filename; write=true) do f
        write(f, "$N\n\n")
        writedlm(f, hcat(vals...))
    end
end
