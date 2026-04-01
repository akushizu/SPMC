module Init

include("sp.jl")
using .Sps

export init

function init(sp, Nmax, m, x, p, fw, number)::Int
		
	println("Initializing fw...")

	Nx, Np, dx, dp = (length(x), length(p), x[2] - x[1], p[2] - p[1])
	
	actualNumber = 0
	count = 0

	for i = 1:Nx
		for j = 1:Np

			Fw = round(Int, number * fw(x[i], p[j]) * dx * dp)

			actualNumber += abs(Fw)

			if number > Nmax

				println("The number of particles has reached the limit.")
				exit(0)

			end

			for k = 1:abs(Fw)

				count += 1

				push!(sp, Sp(count))

				sp[count].m = m
				sp[count].xcell = i
				sp[count].x = x[i] + (rand() - 0.5) * dx
				sp[count].pcell = j
				sp[count].p = p[j] + (rand() - 0.5) * dp
				sp[count].weight = sign(Fw)
				sp[count].time = 0

				sp[count].status = true

			end
		end
	end

	if count == 0

		println("There are no particles in the system.")
		exit(0)

	end

	println("Initial number of particles = $count")

	return count

end

end