using Random
using Dates
include("wmc_muscato.jl")
include("save.jl")
using .Wmc
using .Wmc.Annihilation
using .Wmc.Annihilation.Init
using .Wmc.Annihilation.Init.Sps
using .Save

Nmax = 2000000

hbar = 6.62607015 / 1.602176634 / 2 / pi
m = 0.32
e = 1

Nx = 200
Np = 400
Nt = 500

xrange = [-30, 30]
prange = [-10 * hbar, 10 * hbar]

initialParticles = 160000

#int fanni = 10;

dx = (xrange[2] - xrange[1]) / Nx
dp = (prange[2] - prange[1]) / Np
dt = 0.05

x0 = -15
s0 = 2.852
p0 = 0.7 * hbar
xc = 0
a = 0.3
sigma = 1

function wFunction(x, p)::Float64

	return 1. / pi / hbar * exp(-0.5 * ((x - x0) / s0)^2) * exp(-2. * (s0 * (p - p0) / hbar)^2)
end

function potential(x)::Float64

	return e * a * exp(-0.5 * ((x - xc) / sigma)^2)
end

function wKernel(x, p)::Float64

	return 2. * a * sigma * sqrt(2 * pi) / pi / hbar^2 * exp(-2. * (sigma * p / hbar)^2) * sin(2. * p * x / hbar)
end

cutoff = 8 * hbar
vwh = 2 * a * sigma * sqrt(2 * pi) / pi / hbar^2
gamma = vwh * cutoff

function @main()(args)
	
	Random.seed!(12345)
	
	sp = Sp[]

	x = range(xrange[1] + dx / 2, xrange[2]; step=dx)
	p = range(prange[1] + dp / 2, prange[2]; step=dp)
	t = range(0, Nt * dt; step=dt)

	println("Maximum number of particles allowed = $Nmax\n")

	t0 = time_ns()

	numberOfParticles = init(sp, Nmax, m, x, p, wFunction, initialParticles)

	save(sp, t, 0, x, p)

	t1 = t0
	t2 = time_ns()
	println("Elapsed time is $((t2-t1)/1e8) seconds.")
	println("Total elapsed time is $(canonicalize(Nanosecond(t2-t0)))")
	flush(stdout)

	particleTime = 0

	for i = 1:Nt

		particleTime += dt

		println("\nStep $i of $Nt --- time = $particleTime\n")

		numberOfParticles = wmc(sp, Nmax, x, p, dt, wKernel, cutoff, vwh, gamma)

		if numberOfParticles > 480000
			
			sp = annihilation_a(sp, x, p)
		end
		
		save(sp, t, i, x, p)

		t1 = t2
		t2 = time_ns()
		println("Elapsed time is $((t2-t1)/1e8) seconds.")
		println("Total elapsed time is $(canonicalize(Nanosecond(t2-t0)))")
		flush(stdout)

	end
	
	println("\nOutput files saved\n")
	
end

@main
