module Wmc

using Random

include("annihilation.jl")
using .Annihilation.Init.Sps

export wmc

function propagate(sp, x0, dx, dt)

	sp.x += sp.p / sp.m * dt
	sp.xcell = ceil((sp.x - (x0 - dx / 2)) / dx)
	sp.time += dt

	return
end

function range_check(sp, Nx, Np)

	if (sp.xcell > 0 && sp.xcell <= Nx) && (sp.pcell > 0 && sp.pcell <= Np)
		
		return true
	else

		return false
	end
end

function create_particles(sp, n, currentParticles, newParticles, Nmax, dt, p0, dp, Vw, cutoff, vwh, Gamma)

	if rand() < (Gamma * dt)
	
		ps = (2 * rand() - 1) * cutoff

		if rand() < abs(Vw(sp[n].x, ps)) / vwh
		
			currentParticles += 2
			newParticles += 2

			if currentParticles > Nmax
			
				print("The number of particles has reached the limit.")
				exit(0)
			end

			push!(sp, Sp(currentParticles-1))
			sp[currentParticles-1].m = sp[n].m
			sp[currentParticles-1].x = sp[n].x
			sp[currentParticles-1].xcell = sp[n].xcell
			sp[currentParticles-1].p = sp[n].p + ps
			sp[currentParticles-1].pcell = ceil((sp[currentParticles-1].p - (p0 - dp / 2)) / dp)
			sp[currentParticles-1].weight = sp[n].weight * sign(Vw(sp[n].x, ps))
			sp[currentParticles-1].time = sp[n].time
			sp[currentParticles-1].status = true

			push!(sp, Sp(currentParticles))
			sp[currentParticles].m = sp[n].m
			sp[currentParticles].x = sp[n].x
			sp[currentParticles].xcell = sp[n].xcell
			sp[currentParticles].p = sp[n].p - ps
			sp[currentParticles].pcell = ceil((sp[currentParticles].p - (p0 - dp / 2)) / dp)
			sp[currentParticles].weight = -sp[n].weight * sign(Vw(sp[n].x, ps))
			sp[currentParticles].time = sp[n].time
			sp[currentParticles].status = true

		end
	end

	return (currentParticles, newParticles)
end

function wmc(sp, Nmax, x, p, dt, Vw, cutoff, vwh, Gamma)

	println("Evolving particles...")
	
	Nx, Np, x0, dx, p0, dp = (length(x), length(p), x[1], x[2] - x[1], p[1], p[2] - p[1])

	originalParticles = length(sp)
	
	currentParticles, newParticles, outsideParticles = originalParticles, 0, 0

	for n = 1:originalParticles

		if sp[n].status == false

			outsideParticles += 1
			continue
		end
		
		propagate(sp[n], x0, dx, dt)
		
		if range_check(sp[n], Nx, Np) == false

			sp[n].status = false
			outsideParticles += 1
			continue
		end

		currentParticles, newParticles = create_particles(sp, n, currentParticles, newParticles, Nmax, dt, p0, dp, Vw, cutoff, vwh, Gamma)
	end

	println("Number of created particles = $newParticles")
	println("Total number of particles = $currentParticles")

	for n = originalParticles+1:currentParticles

		if range_check(sp[n], Nx, Np) == false
		
			sp[n].status = false
			outsideParticles += 1
		end
	end

	println("Number of outside particles = $outsideParticles")

	return currentParticles
end

end