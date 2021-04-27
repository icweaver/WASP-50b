### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ 38304fe1-5acb-4785-b7c3-bb08fa5481a6
begin
	import PlutoUI as pl
	using FITSIO, Glob, Statistics, CCDReduction
	using CairoMakie
end

# ╔═╡ fb39c593-86bd-4d4c-b9ec-e5e212a4de98
md"""
# Raw data frames

$(pl.TableOfContents())

In this notebook we will take a quick look at a sample bias subtracted science frame using `CCDReductions.jl`. This is done to verify the mask slit placement, as well as the positions of WASP-50 and its comparison stars on each CCD chip.
"""

# ╔═╡ 8d370ec7-c633-4eb6-9216-f3b88a814b38
md"""
## Load data 📁

For each of the 8 fits files (1 for each chip), we extract just the portion on each chip specified by `DATASEC`, and then subtract out the median bias level measured in the overscan region (defined by `BIASSEC`). We then store the subtracted data in the ``2048 \times 1024 \times 8`` array `cube`:
"""

# ╔═╡ 86b4ace1-a891-41e5-a29f-f7eee5f8fb17
begin
	fpaths_glob = sort(glob("data/raw_frames/ut131219/ift0026c*.fits"))
	
	# Pre-sort chips into IMACS order
	chips = ["c1", "c6", "c2", "c5", "c3", "c8", "c4", "c7"]
	fpaths = [
		filter(s -> occursin(c, s), fpaths_glob)[1]
		for c in chips
	]
	
	# Load science frames and perform bias subtraction
	cube = Array{Float64}(undef, 2048, 1024, 8);	
	for (i, fpath) in enumerate(fpaths)
		ccd = CCDData(fpath)
		dxsec, dysec = Meta.parse(ccd.hdr["DATASEC"]) |> eval
		oxsec, oysec = Meta.parse(ccd.hdr["BIASSEC"]) |> eval
		cube[:, :, i] .= ccd[dysec, dxsec] .- median(ccd[oysec, oxsec])
	end
end

# ╔═╡ 2f1b6036-9ed8-429d-92a9-0f80222c0d68
md"""
## Plot
"""

# ╔═╡ 3a6ab0c0-ba08-4151-9646-c19d45749b9f
let
	fig = Figure(resolution = (800, 600))
	step = 32
	chip_idx = 1
	for j in 1:4, i in 1:2
		if "$(chips[chip_idx])" in ["c5", "c6", "c7", "c8"]
			x_flip = x -> x
			y_flip = reverse
		else
			x_flip = reverse
			y_flip = x -> x
		end
		d = cube[x_flip(begin:step:end), y_flip(begin:step:end), chip_idx]'
		nrows, ncols = size(d)
		global ax, hm = heatmap(
			fig[i, j],
			d,
			colormap = :magma,
			#axis = (title="$(chips[chip_idx])",),
			colorrange = (0, 2_000),
			#interpolate = true,
		)
		annotations!(
			["$(chips[chip_idx])"],
			[Point2f0(2, 0.92*ncols)],
			textsize =0.1*nrows,
			color = :white,
		)
		chip_idx += 1
	end
	axs = reshape(copy(fig.content), 2, 4)
	
	Colorbar(fig[:, end+1], hm, width=20, labelcolor=:black, label="Counts",)
	
	scene = current_axis()
	linkaxes!(axs...)
	hidedecorations!.(axs)
	
	fig
end

# ╔═╡ Cell order:
# ╟─fb39c593-86bd-4d4c-b9ec-e5e212a4de98
# ╟─8d370ec7-c633-4eb6-9216-f3b88a814b38
# ╠═86b4ace1-a891-41e5-a29f-f7eee5f8fb17
# ╟─2f1b6036-9ed8-429d-92a9-0f80222c0d68
# ╠═3a6ab0c0-ba08-4151-9646-c19d45749b9f
# ╠═38304fe1-5acb-4785-b7c3-bb08fa5481a6
