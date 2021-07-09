### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 38304fe1-5acb-4785-b7c3-bb08fa5481a6
begin
	import Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
	
	using AlgebraOfGraphics
	using CSV
	using CairoMakie
	using Colors
	using DataFrames
	using DataFramesMeta
	using Dates
	using DelimitedFiles
	using FITSIO
	using Glob
	using ImageFiltering
	using KernelDensity
	using Latexify
	using Measurements
	using NaturalSort
	using OrderedCollections
	using Printf
	using PyCall
	using Statistics
	using PlutoUI: TableOfContents, Select, Slider, as_svg, with_terminal
end

# ╔═╡ fb39c593-86bd-4d4c-b9ec-e5e212a4de98
md"""
# Raw Frames

$(TableOfContents())

In this notebook we will take a quick look at a sample bias subtracted science frame. We do this to verify the mask slit placement, as well as the positions of WASP-50 and its comparison stars on each CCD chip.

For each of the 8 fits files (1 for each chip), we extract just the portion on each chip specified by `DATASEC`, and then subtract out the median bias level measured in the overscan region (defined by `BIASSEC`). We then store the subtracted data in the ``2048 \times 1024 \times 8`` array `cube`:
"""

# ╔═╡ 86b4ace1-a891-41e5-a29f-f7eee5f8fb17
begin
	fpaths_glob = sort(glob("data/raw/ut131219/ift0026c*.fits"))
	
	# Pre-sort chips into IMACS order
	chips = ["c1", "c6", "c2", "c5", "c3", "c8", "c4", "c7"]
	fpaths = [
		filter(s -> occursin(c, s), fpaths_glob)[1]
		for c in chips
	]
	
	# Load science frames and perform bias subtraction
	cube = Array{Float64}(undef, 2048, 1024, 8);	
	for (i, fpath) in enumerate(fpaths)
		FITS(fpath, "r") do f
			data = read(f[1]) |> permutedims
			dxsec, dysec = read_key(f[1], "DATASEC")[1] |> Meta.parse |> eval
			oxsec, oysec = read_key(f[1], "BIASSEC")[1] |> Meta.parse |> eval
			cube[:, :, i] .= data[dysec, dxsec] .- median(data[oysec, oxsec])
		end
	end
end

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
			colorrange = (0, 2_000),
		)
		
		text!(
			"$(chips[chip_idx])",
			position = (2, 2),
			color = :white,
		)
		chip_idx += 1
	end
	axs = reshape(copy(fig.content), 2, 4)
	
	Colorbar(fig[:, end+1], hm, width=20, labelcolor=:black, label="Counts",)
	
	linkaxes!(axs...)
	hidedecorations!.(axs)
	
	fig #|> as_svg
end

# ╔═╡ f017cddc-5793-4273-a884-05df392149a8
md"""
## Plot configs
"""

# ╔═╡ efe48a46-9ba6-4e26-92a1-0dfe3a7b133c
begin
	const FIG_TALL = (900, 1_200)
	const FIG_WIDE = (1_350, 800)
	#const COLORS = to_colormap(:seaborn_colorblind6, 8)[[8, 6, 4, 1]]
	const COLORS = parse.(Colorant,
		[
			"#a6cee3",  # Cyan
			"#fdbf6f",  # Yellow
			"#ff7f00",  # Orange
			"#1f78b4",  # Blue
			# "plum",
			# "#956cb4",  # Purple
			# "mediumaquamarine",
			# "#029e73",  # Green
		]
	)
	
	set_aog_theme!()
	update_theme!(
		Theme(
			Axis = (xlabelsize=18, ylabelsize=18,),
			Label = (textsize=18,),
			Lines = (linewidth=3, cycle=Cycle([:color, :linestyle], covary=true)),
			Scatter = (linewidth=10,),
			palette = (color=COLORS, patchcolor=[(c, 0.35) for c in COLORS]),
			fontsize = 18,
			rowgap = 0,
			colgap = 0,
		)
	)
	
	COLORS
end

# ╔═╡ 4480ae72-3bb2-4e17-99be-28afc756332a
md"""
## Packages
"""

# ╔═╡ 6000db3d-0798-4f76-be31-617d43406b54
html"""
<style>
#launch_binder {
	display: none;
}
body.disable_ui main {
		max-width : 95%;
	}
@media screen and (min-width: 1081px) {
	body.disable_ui main {
		margin-left : 10px;
		max-width : 72%;
		align-self: flex-start;
	}
}
</style>
"""

# ╔═╡ Cell order:
# ╟─fb39c593-86bd-4d4c-b9ec-e5e212a4de98
# ╠═86b4ace1-a891-41e5-a29f-f7eee5f8fb17
# ╠═3a6ab0c0-ba08-4151-9646-c19d45749b9f
# ╟─f017cddc-5793-4273-a884-05df392149a8
# ╠═efe48a46-9ba6-4e26-92a1-0dfe3a7b133c
# ╟─4480ae72-3bb2-4e17-99be-28afc756332a
# ╠═38304fe1-5acb-4785-b7c3-bb08fa5481a6
# ╟─6000db3d-0798-4f76-be31-617d43406b54
