module PureEOS

# package code goes here
using Roots
using Dierckx

#water
export waterSaturate, waterSinglephase, waterMax
export waterIceBoundary
export waterCritical, waterTriple
export waterVisc, waterViscMax
export waterTherm

#Lennard Jones
export LJSinglephase

#ethanol
export ethanolSaturate, ethanolCritical
export ethanolTriple, ethanolSinglephase
export ethanolTherm, ethanolVisc

#methanol
export methanolSaturate, methanolCritical
export methanolTriple, methanolSinglephase

#ethane
export ethaneSinglephase,ethaneCritical
export ethaneSaturate,ethaneTriple

#stocmayer
export StECH

#mathane
export methaneSaturate, methaneSinglephase
export methaneCritical, methaneTriple

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
	lRo::Float64
end

mutable struct viscmax
	p::Float64
	ro::Float64
end

mutable struct bound
	pUp::Float64
	roUp::Float64
	pDown::Float64
	roDown::Float64
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

function newBound()::bound
	return bound(0.0, 0.0, 0.0, 0.0)
end

include("water.jl")
include("LJ.jl")
include("ethanol.jl")
include("ethane.jl")
include("Stockmayer.jl")
include("methanol.jl")
include("methane.jl")
end # module
