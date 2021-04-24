### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ 38304fe1-5acb-4785-b7c3-bb08fa5481a6
begin
	import PlutoUI as pl
	using FITSIO, Glob, Statistics
	using CairoMakie
    # using WGLMakie, JSServe
    # Page(exportable = true, offline = true)
end

# ╔═╡ eb500902-d124-4d17-bbba-b60855035f2a
cube = Array{UInt16}(undef, 1088, 2112, 8);

# ╔═╡ 86b4ace1-a891-41e5-a29f-f7eee5f8fb17
fpaths_glob = sort(glob("data/raw_frames/ut131219/ift0026c*.fits"));

# ╔═╡ 35c9020e-0a57-4006-95b3-43d48510ddf1
chips = ["c1", "c6", "c2", "c5", "c3", "c8", "c4", "c7"]

# ╔═╡ d4fed004-9747-4765-b0b0-995d8d54c30d
fpaths = [
	filter(s -> occursin(c, s), fpaths_glob)[1]
	for c in chips
]

# ╔═╡ b37c402b-738e-42fa-aecb-9dafe706f832
for (i, fpath) in enumerate(fpaths)
	cube[:, :, i] .= FITS(fpath, "r") do f
		read(f[1])
	end
end

# ╔═╡ 3a6ab0c0-ba08-4151-9646-c19d45749b9f
let
	fig = Figure(resolution = (1_000, 500))
	
	chip_idx = 1
	for j in 1:4, i in 1:2
		if "$(chips[chip_idx])" ∈ ["c5", "c6", "c7", "c8"]
			f = x -> x
			g = reverse
		else
			f = reverse #x -> x
			g = x -> x
		end
		heatmap(
			fig[i, j],
			cube[g(begin:16:end), f(begin:16:end), chip_idx],
			colormap = :magma,
			axis = (title="$(chips[chip_idx])",),
			colorrange = (0, 4_000)
		)
		chip_idx += 1
	end
	
	axs = reshape(copy(fig.content), 2, 4)
	
# 	for ax in axs
# 		ax.colorrange = (1, 3)
# 	end
	
	scene = current_axis()
	linkaxes!(axs...)
	hidedecorations!.(axs)
	
# 	fig[1:2, 0] = Label(fig, "Counts", rotation=π/2)
# 	fig[end+1, 2:3] = Label(fig, "Index")
	
	fig
end

# ╔═╡ Cell order:
# ╠═eb500902-d124-4d17-bbba-b60855035f2a
# ╠═86b4ace1-a891-41e5-a29f-f7eee5f8fb17
# ╠═d4fed004-9747-4765-b0b0-995d8d54c30d
# ╠═35c9020e-0a57-4006-95b3-43d48510ddf1
# ╠═b37c402b-738e-42fa-aecb-9dafe706f832
# ╠═3a6ab0c0-ba08-4151-9646-c19d45749b9f
# ╠═38304fe1-5acb-4785-b7c3-bb08fa5481a6
