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
	F1::SparseMatrixCSC{T,Int64}
	F2::SparseMatrixCSC{T,Int64}
	xx::Matrix{Float64}
end

mutable struct Kernel_1D{T<:Real} <: Kernel
	F::SparseMatrixCSC{T,Int64}
end


function Kernel(x::Array{Vector{Float64}}, 
		  xi::Array{Vector{Float64}}, Battrib::Symbol=:B1)



	# constructor for y=Ax
	function interp!(y, x, pa)
		for i in eachindex(y)
			y[i]=0.0
		end
		Interpolation.interp_spray!(x, y, pa, :interp)
	end

	if(length(x)==1)
		pa=P_core(x,xi,Battrib)
		F1=LinearMap((y,x)->interp!(y,x,pa), 
		      nothing, length(xi[1]), length(x[1]), ismutating=true)
		F1=SparseArrays.sparse(F1)
		return Kernel_1D(F1)
	elseif(length(x)==2)
		pa1=P_core(x[2:2],xi[2:2],Battrib)
		F1=LinearMap((y,x)->interp!(y,x,pa1), 
		      nothing, length(xi[2]), length(x[2]), ismutating=true)
		F1=SparseArrays.sparse(F1)
		pa2=P_core(x[1:1],xi[1:1],Battrib)
		F2=LinearMap((y,x)->interp!(y,x,pa2), 
		      nothing, length(xi[1]), length(x[1]), ismutating=true)
		F2=SparseArrays.sparse(F2)
		xx=zeros(size(F1,1), length(x[1]))
		return Kernel_2D(F1,F2,xx)
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
		A_mul_B!(pa.xx, pa.F1, y)
		A_mul_Bc!(yi, pa.xx, pa.F2)
	elseif(attrib == :spray)
		A_mul_B!(pa.xx, yi, pa.F2)
		Ac_mul_B!(y, pa.F1, pa.xx)
	end
end



include("Misc.jl")

include("Core.jl")
include("Weights.jl")

end # module
