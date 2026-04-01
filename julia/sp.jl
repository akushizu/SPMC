module Sps

export Sp
export Spn

mutable struct Sp

	number::Int
	id::Int

	m::Float64
	x::Float64
	xcell::Int
	p::Float64
	pcell::Int
	weight::Int

	time::Float64

	status::Bool

	Sp(number) = new(number)
end

mutable struct Spn

	dimension::Int

	xn::Vector{Float64}
	xncell::Vector{Int}
	pn::Vector{Float64}
	pncell::Vector{Int}

	number::Int
	id::Int

	m::Float64
	x::Float64
	xcell::Int
	p::Float64
	pcell::Int
	weight::Int

	time::Float64

	status::Bool

	Spn(number) = new(number)
end

end