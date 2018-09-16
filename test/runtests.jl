using Interpolation
using BenchmarkTools
using Test

@testset "no extrapolations" begin
	nx=Array(range(21,100,length=80));
	mx=Array(range(1,100,length=100));
	ny=ones(length(nx));
	c=randn()
	my=c.*ones(length(mx))

	pa=Interpolation.P_core([nx], [mx], :B1)

	Interpolation.interp_spray!(ny, my, pa, :interp)
	@test all(my[1:20] .== c)

	nx=Array(range(3,7,length=64));
	nz=Array(range(3,7,length=95));
	mx=Array(range(1,10,length=10));
	mz=Array(range(1,10,length=10));
	ny=ones(length(nz), length(nx));
	c=randn()
	my=c.*ones(length(mz), length(mx))

	pa=Interpolation.P_core([nx, nz], [mx, mz], :B1)

	Interpolation.interp_spray!(ny, my, pa, :interp)

	@test all(my[my .≠1] .== c)

end


# also testing behaviour when out of bounds!!
nx=Array(range(1,20,length=100));
mx=Array(range(3,25,length=200));
nz=Array(range(1,20,length=30));
mz=Array(range(5,27,length=40));

println("=====================================")
## 1D
ny=randn(length(nx));
# interpolate (my=f(ny))
my=zeros(length(mx));

for Battrib in [:B1, :B2]

	pa=Interpolation.Kernel([nx], [mx], Battrib)
	println("=====================================")
	## 1D
	ny=randn(length(nx));
	# interpolate (my=f(ny))
	my=zeros(length(mx));
	@time Interpolation.interp_spray!(ny, my, pa, :interp)

	myp = randn(length(my));

	# spray (nyp=fˣ(myp))
	nyp=zeros(length(nx));
	@time Interpolation.interp_spray!(nyp, myp, pa, :spray)

	# dot product test
	@test dot(my, myp) ≈ dot(ny, nyp)
end


@testset "Interpolation on the same grid should not change for B1" begin
	Battrib=:B1
	pa=Interpolation.Kernel([nx, nz], [nx, nz], Battrib)
	println("=====================================")
	ny=randn(length(nz), length(nx));
	my=zeros(length(nz), length(nx));
	@time Interpolation.interp_spray!(ny, my, pa, :interp)
	@test ny≈my


	nyvec=vcat(vec(ny),vec(ny))
	myvec=zeros(length(my)*2)
	@time Interpolation.interp_spray!(nyvec, myvec, pa, :interp, 2)
	@test nyvec≈myvec



	myp = randn(size(my));
	nyp=zeros(length(nz), length(nx));
	@time Interpolation.interp_spray!(nyp, myp, pa, :spray)
	@test dot(my, myp) ≈ dot(ny, nyp)
end


for Battrib in [:B1, :B2]
	pa=Interpolation.Kernel([nx, nz], [mx, mz], Battrib)
	## 2D
	ny=randn(length(nz), length(nx));
	# interpolate
	my=zeros(length(mz), length(mx));
	@time Interpolation.interp_spray!(ny, my, pa, :interp)


	myp = randn(size(my));

	# spray
	nyp=zeros(length(nz), length(nx));
	@time Interpolation.interp_spray!(nyp, myp, pa, :spray)

	# dot product test
	@test dot(my, myp) ≈ dot(ny, nyp)


	nmod=3
	ny=randn(length(nz)*length(nx)*nmod);
	# interpolate
	my=zeros(length(mz)*length(mx)*nmod);
	@time Interpolation.interp_spray!(ny, my, pa, :interp,nmod)


	myp = randn(size(my));

	# spray
	nyp=zeros(length(nz)*length(nx)*nmod);
	@time Interpolation.interp_spray!(nyp, myp, pa, :spray,nmod)

	# dot product test
	@test dot(my, myp) ≈ dot(ny, nyp)




end

