module PureEOS
export greet

using Roots
using Dierckx

export LJCritical
export LJTriple
export LJSinglephase
export LJ2Singlephase

mutable struct singlephase
	T::Float64
	p::Float64
	ro::Float64
	z::Float64
	s::Float64
	s0::Float64
	ds::Float64
	g::Float64
	g0::Float64
	dg::Float64
	f::Float64
	f0::Float64
	df::Float64
	u::Float64
	du::Float64
	h::Float64
	dh::Float64
	w::Float64
	cp::Float64
	cv::Float64
	dpdt::Float64
end

mutable struct critical
	T::Float64
	p::Float64
	ro::Float64
end

mutable struct triple
	T::Float64
	p::Float64
	vRo::Float64
	lRo::Float64
end

function newSingle()::singlephase
	return singlephase(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
end

function newCritical():critical
	return critical(0,0,0)
end

function newTriple():triple
	return triple(0,0,0,0)
end

include("LJ.jl")
include("LJ2.jl")

end # module PureEOS
