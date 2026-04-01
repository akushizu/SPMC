module Annihilation

using Random

include("init.jl")
using .Init.Sps

export annihilation, annihilation_a

function density(sp, Nx, Np, Fw)

	for n in eachindex(sp)

		i = sp[n].xcell
		j = sp[n].pcell

		if (i > 0 && i <= Nx) && (j > 0 && j <= Np)

			Fw[i, j] += sp[n].weight
		end
	end

	return
end

function comp(sp1, sp2)

	if sp1.status == true && sp2.status == false
		return true
	
	elseif sp1.status == false && sp2.status == false
		return false
	
	elseif sp1.status == true && sp2.status == true

		if sp1.xcell < sp2.xcell
			return true
		
		elseif sp1.xcell > sp2.xcell
			return false
		
		elseif sp1.pcell < sp2.pcell
			return true
		
		elseif sp1.pcell > sp2.pcell
			return false
		
		else
			return false
		end
		
	else
		return false
	end
end

function check_range(sp, Nx, Np)

	for n in eachindex(sp)

		if sp[n].status == false
			continue

		else 
			if !(sp[n].xcell > 0 && sp[n].xcell <= Nx) || !(sp[n].pcell > 0 && sp[n].pcell <= Np)

				sp[n].status = false
			end
		end
	end

	return
end

function sparse_matrix(sp, sparse)

	if sp[1].status != false

		push!(sparse, ([sp[1].xcell, sp[1].pcell, sp[1].weight]))
	end
	
	for i in eachindex(sp)[2:end]

		if sp[i].status == false
			continue
		
		elseif sp[i].xcell != sp[i - 1].xcell

			if sparse[end][3] == 0

				pop!(sparse)
			end

			push!(sparse, ([sp[i].xcell, sp[i].pcell, sp[i].weight]))
		
		elseif sp[i].pcell != sp[i - 1].pcell

			if sparse[end][3] == 0

				pop!(sparse)
			end

			push!(sparse, ([sp[i].xcell, sp[i].pcell, sp[i].weight]))
		
		else

			sparse[end][3] += sp[i].weight
		end
	end
end

function annihilation(sp, x, p)
	#grid, O(n)

	println("\nAnnihilating particles...")
	
	Nx, Np, dx, dp = (length(x), length(p), x[2] - x[1], p[2] - p[1])

	println("# of particles before annihilation = $(length(sp))")
	
	Fw = zeros(Int, length(x), length(p))
	
	density(sp, Nx, Np, Fw)

	count = 0

	for i = 1:Nx
		for j = 1:Np
			for k = 1:abs(Fw[i, j])

				count += 1
				
				sp[count].number = count

				sp[count].xcell = i
				sp[count].x = x[i] + (rand() - 0.5) * dx
				sp[count].pcell = j
				sp[count].p = p[j] + (rand() - 0.5) * dp
				sp[count].weight = sign(Fw[i, j])

				sp[count].status = true
			end
		end
	end

	for i = eachindex(sp)[count+1:end]

		pop!(sp)
	end
					
	println("# of particles after annihilation  = $count\n")

	return sp
end

function annihilation_a(sp, x, p)
	#sort, O(n log n)

	println("\nAnnihilating particles...")

	Nx, Np, dx, dp = (length(x), length(p), x[2] - x[1], p[2] - p[1])

	println("# of particles before annihilation = $(length(sp))")

	check_range(sp, Nx, Np)

	sort!(sp, lt=comp)

	sparse = []

	sparse_matrix(sp, sparse)

	count = 0

	m = sp[1].m
	t = sp[1].time

	empty!(sp)

	for i in eachindex(sparse)
		for n = 1:abs(sparse[i][3])

			count += 1

			push!(sp, Sp(count))

			sp[count].m = m
			sp[count].xcell = sparse[i][1]
			sp[count].x = x[sparse[i][1]] + (rand() - 0.5) * dx
			sp[count].pcell = sparse[i][2]
			sp[count].p = p[sparse[i][2]] + (rand() - 0.5) * dp
			sp[count].weight = sign(sparse[i][3])
			sp[count].time = t

			sp[count].status = true
		end
	end

	println("# of particles after annihilation  = $(length(sp))\n")

	return sp
end

end