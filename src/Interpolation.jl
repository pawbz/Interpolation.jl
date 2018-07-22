__precompile__()

"""
## TODO:
* add dimension checks to interp_spray!
Reference: https://www.ibiblio.org/e-notes/Splines/bezier.html
"""
module Interpolation


using LinearMaps


abstract type Kernel end

mutable struct Kernel_2D{T<:Real} <: Kernel
	F::SparseMatrixCSC{T,Int64}
	x::Vector{T} # just allocating input
	y::Vector{T} # just allocating output
end

mutable struct Kernel_1D{T<:Real} <: Kernel
	F::SparseMatrixCSC{T,Int64}
end


function Kernel(x::Array{Vector{Float64}}, 
		  xi::Array{Vector{Float64}}, Battrib::Symbol=:B1)


	pa=P_core(x,xi,Battrib)

	xin=zeros(flipdim(pa.nx,1)...)
	xout=zeros(flipdim(pa.nxi,1)...)
	# constructor for y=Ax
	function interp!(y, x)
		for i in eachindex(y)
			y[i]=0.0
		end
		copy!(xin,x)
		Interpolation.interp_spray!(xin, xout, pa, :interp)
		copy!(y,xout)
	end

	F=LinearMap(interp!, nothing, prod(pa.nxi), prod(pa.nx), ismutating=true)

	F=SparseArrays.sparse(F)

	if(length(pa.nx)==1)
		return Kernel_1D(F)
	elseif(length(pa.nx)==2)
		x=zeros(prod(pa.nx))
		y=zeros(prod(pa.nxi))
		return Kernel_2D(F,x,y)
	else
		error("dimension should be <=2")

	end

end

function interp_spray!(y, yi, pa::Kernel_1D, attrib)
	if(attrib == :interp)
		A_mul_B!(yi, pa.F, y)
	elseif(attrib == :spray)
		Ac_mul_B!(y, pa.F, yi)
	end
end
function interp_spray!(y, yi, pa::Kernel_2D, attrib)
	if(attrib == :interp)
		copy!(pa.x,y)
		A_mul_B!(pa.y, pa.F, pa.x)
		copy!(yi,pa.y)
	elseif(attrib == :spray)
		copy!(pa.y,yi)
		Ac_mul_B!(pa.x, pa.F, pa.y)
		copy!(y,pa.x)
	end
end



include("Misc.jl")

include("Core.jl")
include("Weights.jl")

end # module
