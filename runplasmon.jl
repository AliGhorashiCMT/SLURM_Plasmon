frac = parse(Float64, ARGS[1])
eV=1/27.2
using PyCall
np=pyimport("numpy")

function hwannier(wannier_file::AbstractString, cell_map_file::AbstractString, nbands::Integer) 
    cell_map_numlines = countlines(cell_map_file)
    Hwannier = permutedims(reshape(np.loadtxt(wannier_file), (cell_map_numlines, nbands, nbands)), [1, 3, 2])
    return Hwannier
end
function wannier_bands(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, k::Vector{<:Real}, nbands) 
    phase = np.exp(2im*np.pi*cell_map*k); H = np.tensordot(phase, Hwannier, axes=1); E, _=np.linalg.eigh(H);
    return E ./eV 
end
function wannier_bands(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, k::Vector{<:Real}) 
    phase = np.exp(2im*np.pi*cell_map*k); H = np.tensordot(phase, Hwannier, axes=1); E, _=np.linalg.eigh(H);
    return E[1] ./eV 
end

HWannier=hwannier("wannierUp.txt", "wannierUp.map.txt", 1)
cellmap=np.loadtxt("wannierUp.map.txt")
function im_polarization(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, 
    q::Vector{<:Real}, mu::Real; spin::Integer=1, mesh::Integer=100, histogram_width::Integer=100) 
    Polarization_Array=zeros(histogram_width*100)
    V=86.65976167500001
    qnormalized=q
    for (i, j) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        kvector=[i/mesh, j/mesh, 0]
        E1 = wannier_bands(HWannier, cell_map, kvector)
        E2 = wannier_bands(HWannier, cell_map, kvector+qnormalized)
        f1 = np.heaviside(mu-E1, 0.5)
        f2 = np.heaviside(mu-E2, 0.5)
        DeltaE = E2-E1
        DeltaE >0 || continue
        Polarization_Array[round(Int, histogram_width*DeltaE+1)] += π*(f2-f1)/V*(1/mesh)^2*histogram_width*spin
    end
    return Polarization_Array
end
qnorm=[2/3, -1/3, 0]*frac
pols=im_polarization(HWannier, cellmap, qnorm, -2.3, spin=1, mesh=600, histogram_width=1000);
for pol in pols
	println(pol)
end

