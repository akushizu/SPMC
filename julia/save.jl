module Save

using Printf

export save

function saveRho(sp, rho, t, ti, x, Np)

	Nx, Nt = (length(x), length(t))

	for n in eachindex(sp)

		i = sp[n].xcell
		j = sp[n].pcell

		if (i > 0 && i <= Nx) && (j > 0 && j <= Np)

			rho[i] += sp[n].weight
		end
	end

	if ti == 0

		open("rho.dat", "w") do den

			println(den, "x\trho")
		end
	end

	open("rho.dat", "a") do den

		for i = 1:20

			print(den, "-")
		end

		println(den, "\n# = $ti of $(Nt-1)")
		println(den, "t = $(t[ti+1])")
		println(den, "")

		for i = 1:Nx

			@printf(den, "%.3f\t%d\n", x[i], rho[i])
		end

		println(den, "")
	end

	return
end

function saveFw(sp, Fw, t, ti, x, p)

	Nx, Np, Nt = (length(x), length(p), length(t))

	for n in eachindex(sp)

		i = sp[n].xcell
		j = sp[n].pcell

		if (i > 0 && i <= Nx) && (j > 0 && j <= Np)

			Fw[i, j] += sp[n].weight
		end
	end

	if ti == 0 

		open("Fw.dat", "w") do fw

			println(fw, "x")

			for i = 1:Nx

				@printf(fw, "%.3f\t", x[i])
			end

			println(fw, "\np")

			for i = 1:Np

				@printf(fw, "%.3f\t", p[i])
			end

			println(fw, "\nFw")
		end
	end

	open("Fw.dat", "a") do fw

		for i = 1:20

			print(fw, "-")
		end

		println(fw, "\n# = $ti of $(Nt-1)")
		println(fw, "t = $(t[ti+1])")
		println(fw, "")

		for j = Np:-1:1

			for i = 1:Nx

				print(fw, "$(Fw[i, j])\t")
			end

			println(fw, "")
		end

		println(fw, "\n")
	end

	return
end

function save(sp, t, ti, x, p)

	Fw = zeros(Int, length(x), length(p))
	rho = zeros(Int, length(x))

	saveRho(sp, rho, t, ti, x, length(p))
	saveFw(sp, Fw, t, ti, x, p)
	
	return
end

end