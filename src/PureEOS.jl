module PureEOS

# package code goes here
using Roots

export waterSaturate, waterSinglephase
export waterCritical, waterTriple
export waterVisc, waterViscMax
export waterTherm
export LJSinglephase

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

mutable struct saturate
	T::Float64
	p::Float64
	vRo::Float64
	lRo::Float64
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
	lro::Float64
end

mutable struct viscmax
	p::Float64
	ro::Float64
end

function newSingle()::singlephase
	return singlephase(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
end

function newSaturate()::saturate
	return saturate(0,0,0,0)
end

function newCritical():critical
	return critical(0,0,0)
end

function newTriple():triple
	return triple(0,0,0,0)
end

function newViscMax():viscmax
	return viscmax(0,0)
end

include("water.jl")
include("LJ.jl")
end # module
