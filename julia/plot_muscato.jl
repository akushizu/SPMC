# %%
using Pkg
Pkg.activate(".")
Pkg.add("JLD2")
Pkg.add("GLMakie")

using GLMakie
using JLD2

# %%

adata = load("adata.jld2", "adata")

fig = Figure()

record(fig, "muscato_rho.mp4"; framerate = 10) do v
    open("rho.dat", "r") do f
        
        readline(f)
        readline(f)

        N = 0.

        while !eof(f)

            _, _, i, _, n = split(readline(f))
            i = parse(Int, i)

            _, _, t = split(readline(f))
            t = parse(Float64, t)

            readline(f)

            line = readline(f)
            C = []

            while !isempty(line)

                C = isempty(C) ? parse.(Float64, split(line)) : [C parse.(Float64, split(line))]

                line = readline(f)

            end

            readline(f)

            x = C[1, :]
            rho = C[2, :]

            if i == 0
                
                N = sum(rho)
                
            end

            rho = rho / N / (x[2] - x[1])

            xa  = adata[i + 1].x
            amp = adata[i + 1].Amp

            empty!(fig)

            ax = Axis(fig[1, 1],
                title = "t = $t fs", xlabel = "x (nm)", ylabel = "density")
            xlims!(ax, -30, 30)

            p1 = lines!(ax, x, rho; linewidth = 2)
            lines!(ax, xa, 0.05 .* exp.(-(xa .- 0).^2 / 2 / 1^2))
            p2 = lines!(ax, xa, abs2.(amp))
            axislegend(ax, [p1, p2], ["Monte Carlo", "backward difference"])
            
            recordframe!(v)

        end

    end
end

save("rho.png", fig)

# %%

fig = Figure()

record(fig, "muscato_fw.mp4"; framerate = 10) do v
    open("Fw.dat", "r") do f

        readline(f)

        x = parse.(Float64, split(readline(f)))

        readline(f)

        p = reverse(parse.(Float64, split(readline(f))))

        readline(f)
        readline(f)

        N = 0.

        while !eof(f)

            _, _, i, _, n = split(readline(f))
            i = parse(Int, i)

            _, _, t = split(readline(f))
            t = parse(Float64, t)

            readline(f)
            
            line = readline(f)
            fw = []

            while !isempty(line)

                fw = isempty(fw) ? parse.(Float64, split(line)) : [fw parse.(Float64, split(line))]

                line = readline(f)
                
            end

            readline(f)
            readline(f)

            if i == 0
            
                N = sum(fw)
        
            end

            fw = fw ./ N ./ (x[2] - x[1]) ./ (p[1] - p[2])

            empty!(fig)

            ax = Axis(fig[1, 1], title = "t = $t fs", xlabel = "x (nm)", ylabel = "p",
                        limits = ((-30, 30), (p[end] - 0.5 * (p[1] - p[2]) , p[1] + 0.5 * (p[1] - p[2]))))
            contourf!(ax, x, p, fw, colormap = :bwr, levels = range(-0.3, 0.3, length=11),
                        extendhigh = :auto, extendlow = :auto)
            contour!(ax, x, p, fw, color = :black, levels = range(-0.3, 0.3, length=11))
            Colorbar(fig[1, 2], colormap = :bwr, colorrange = (-0.3, 0.3),
                        highclip = :red, lowclip = :blue)

            recordframe!(v)

        end

    end
end

save("Fw.png", fig)
