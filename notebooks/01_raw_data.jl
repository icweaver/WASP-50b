### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 3433ed02-c27c-4fe5-bfda-a5108a58407c
begin
	import Pkg
	Pkg.activate(Base.current_project())
	
	using AlgebraOfGraphics
	using CSV
	using CairoMakie
	using CCDReduction: fitscollection, getdata, CCDData
	using Colors
	using DataFrames
	using DataFrameMacros
	using Dates
	using DelimitedFiles
	using Glob
	using ImageFiltering
	import CairoMakie.Makie.KernelDensity: kde
	using Latexify
	using LaTeXStrings
	using Measurements
	using Measurements: value, uncertainty
	using NaturalSort
	using OrderedCollections
	using Printf
	using Statistics
	using PlutoUI: TableOfContents, Select, Slider, as_svg, with_terminal
	using Unitful
	
	# Python setup
	ENV["PYTHON"] = "/home/mango/miniconda3/envs/WASP-50b/bin/python"
	Pkg.build("PyCall")
	using PyCall
	
	##############
	# PLOT CONFIGS
	##############
	const FIG_TALL = (900, 1_200)
	const FIG_WIDE = (1_350, 800)
	
	set_aog_theme!()
	update_theme!(
		Theme(
			Axis = (xlabelsize=18, ylabelsize=18,),
			Label = (textsize=18,  padding=(0, 10, 0, 0)),
			Text = (; font=AlgebraOfGraphics.firasans("Medium")),
			Lines = (linewidth=3, cycle=Cycle([:color, :linestyle], covary=true)),
			Scatter = (linewidth=10,),
			fontsize = 18,
			rowgap = 0,
			colgap = 0,
		)
	)
	
	COLORS = Makie.wong_colors()
end

# ╔═╡ fb39c593-86bd-4d4c-b9ec-e5e212a4de98
md"""
# Raw Data

$(TableOfContents())

In this notebook we will summarize the raw data and view a sample frame from each instrument before moving on to the automated data reduction.

!!! note "rclone"

	```	
	rclone sync -P drive_cfa,shared_with_me:LDSS3/LDSS3_e150927 drive_ACCESS:papers/WASP-50b/data/raw/LDSS3/ut150927 --include "*606*"
	```
	
	```
	rclone sync -P drive_ACCESS:papers/WASP-50b/data/raw data/raw
	```
"""

# ╔═╡ 7e111d10-0aa8-47bd-9ea8-186c1aecc321
md"""
## Data overview
"""

# ╔═╡ 1e0db1d5-86b7-475c-897b-b0054575a5fa
begin
	DATA_DIR = "data/raw/IMACS/ut131219"
	
	df = fitscollection(
		DATA_DIR;
		recursive=false,
		exclude=r"^((?!ift[0-9]{4}c1.fits).)*$",
		exclude_key=("", "COMMENT"),
	)
	
	@. df[!, "TIMESTAMP"] = DateTime(
		parse(Date, df[!, "UT-DATE"]), parse(Time, df[!, "UT-TIME"])
	)
	
	@. df[!, "JD"] = datetime2julian(df.TIMESTAMP) - 2.456e6
			
	# Select header items to show
	cols = [
		"FILENAME",
		"TIMESTAMP",
		"JD",
		"OBJECT",
		"DISPERSR",
		"SLITMASK",
		"SPEED",
		"FILTER",
		"EXPTIME",
		"BINNING",
		"UT-DATE",
		"UT-TIME",
		"RA",
		"DEC",
		"EQUINOX",
		"EPOCH",
		"AIRMASS",
		"OBSERVER",
	]
	
	night_log = select(df, cols)
end

# ╔═╡ 9257eef5-c18a-4226-9ff3-fc0ea90a1262
combine(groupby(df, :OBJECT), nrow => :count)

# ╔═╡ c3def15d-d41a-4e3f-bb49-d57d46a6474e
obs_cols = [:TIMESTAMP, :JD, :EXPTIME, :AIRMASS, :FILTER, :SLITMASK, :DISPERSR]

# ╔═╡ fbb58e2e-6dd4-42f5-acb1-037b2cc231c6
df_sci = df[occursin.("sci", lowercase.(df[!, :OBJECT])), obs_cols];

# ╔═╡ 3ca865bb-078a-4bd5-a632-a4b4b9d6e1e9
combine(
	groupby(
		df_sci,
		[:EXPTIME, :FILTER, :SLITMASK, :DISPERSR]
	),
	nrow => :count,
)

# ╔═╡ c3ac2517-3bac-4686-9b1e-16a5c5b5865f
combine(
	groupby(
		df[df[!, :OBJECT] .== "Arc Lamps", :],
		[:FILTER, :SLITMASK, :DISPERSR]
	),
	nrow => :count,
)

# ╔═╡ 2fdfc049-11b1-4972-9c48-e8724765126f
let
	airmass = data(df_sci) * mapping(
		:JD=>"Time", :AIRMASS=>"Airmass",
	)
	draw(airmass)
end

# ╔═╡ 90845d70-35d9-402d-8936-74936b069577
md"""
## IMACS

Let's take a quick look at a sample bias subtracted science frame. We do this to verify the mask slit placement, as well as the positions of WASP-50 and its comparison stars on each CCD chip.

For each of the 8 fits files (1 for each chip), we extract just the portion on each chip specified by `DATASEC`, and then subtract out the median bias level measured in the overscan region (defined by `BIASSEC`). We then store the subtracted data in the ``2048 \times 1024 \times 8`` array `cube`:
"""

# ╔═╡ 86b4ace1-a891-41e5-a29f-f7eee5f8fb17
begin
	fpaths_glob = sort(glob("data/raw/IMACS/ut131219/ift0026c*.fits"))
	
	# Pre-sort chips into IMACS order
	chips = ["c1", "c6", "c2", "c5", "c3", "c8", "c4", "c7"]
	fpaths = [
		filter(s -> occursin(c, s), fpaths_glob)[1]
		for c ∈ chips
	]
	
	# Load science frames and perform bias subtraction
	cube = Array{Float64}(undef, 2048, 1024, 8);	
	for (i, fpath) ∈ enumerate(fpaths)
		ccd = CCDData(fpath)
		dxsec, dysec = ccd.hdr["DATASEC"] |> Meta.parse |> eval
		oxsec, oysec = ccd.hdr["BIASSEC"] |> Meta.parse |> eval
		cube[:, :, i] .= ccd[dysec, dxsec] .- median(ccd[oysec, oxsec])
	end
end

# ╔═╡ 26feb668-4e7e-4a9d-a2b6-a5dac81e3ab7
coords = CSV.read("data/raw/IMACS/ut131219/WASP50.coords", DataFrame;
	delim = ' ',
	ignorerepeated = true,
	header = ["target", "chip", "x", "y"],
)

# ╔═╡ 3a6ab0c0-ba08-4151-9646-c19d45749b9f
let
	fig = Figure(resolution = (800, 600))
	step = 1 # For quick testing
	chip_idx = 1
	hms = []
	axs = []
	for j ∈ 1:4, i ∈ 1:2
		# Set chip orientation based on IMACS conventions
		if "$(chips[chip_idx])" ∈ ["c5", "c6", "c7", "c8"]
			xflip = reverse
			yflip = reverse
			xreversed = false
			yreversed = true
		else
			xflip = x -> x
			yflip = x -> x
			xreversed = false
			yreversed = true
		end
		ax = Axis(fig[i, j], xreversed=xreversed, yreversed=yreversed)
		d = cube[yflip(begin:step:end), xflip(begin:step:end), chip_idx]'
		nrows, ncols = size(d)
		hm = heatmap!(
			ax,
			d,
			colormap = :magma,
			colorrange = (0, 2_000),
		)

		# Label objects on chip
		chip = chips[chip_idx]
		for obj in eachrow(@subset coords :chip == chip)
			coord = (obj.x, obj.y) ./ step
			text!(ax, obj.target;
				position = coord,
				font = AlgebraOfGraphics.firasans("Light"),
				textsize = 11,
				align = (:center, :baseline),
				color = :white,
			)
		end
		
		# Label chip
		text!(
			"$(chip)",
			position = (50, 2_000),
			color = :yellow,
		)
		chip_idx += 1
		push!(axs, ax)
		push!(hms, hm)
	end
	axs = reshape(axs, 2, 4)
	
	Colorbar(fig[:, end+1], hms[1], width=20, labelcolor=:black, label="Counts",)
	
	linkaxes!(axs...)
	hidedecorations!.(axs)
	
	path = "../../ACCESS_WASP-50b/figures/frames"
	mkpath(path)
	save("$(path)/sci_imacs.png", fig)
	
	fig #|> as_svg
end

# ╔═╡ a0d54171-3d3b-472f-baca-2883924ec75e
let
	chip_idx = 6
	step = 16
	d = cube[(begin:step:end), reverse(begin:step:end), chip_idx]
	
	fig = Figure(resolution = (400, 600))
	ax = Axis(fig[1, 1], xreversed=false, yreversed=true)

	heatmap!(ax, d', colorrange = (0, 2_000), colormap=:magma)

	chip = chips[chip_idx]
	for obj in eachrow(@subset coords :chip == chip)
		coord = (obj.x, obj.y) ./ step
		scatter!(ax, coord)
		text!(ax, obj.target, position=coord, color=:white,
		#font = AlgebraOfGraphics.firasans("Light")
		)
	end

	text!(ax, "c8", position = (4, 20), color=:yellow, space=:screen)

	fig
end

# ╔═╡ 06a834f0-8c90-4013-af34-725166970969
md"""
## LDSS3
"""

# ╔═╡ 5c6e6f7b-70e0-49a8-b064-60dcf1440223
begin
	DATA_DIR_LDSS3 = "data/raw/LDSS3/ut150927"
	
	df_LDSS3 = fitscollection(
		DATA_DIR_LDSS3;
		recursive=false,
		exclude=r"^((?!ccd[0-9]{4}c1.fits).)*$",
		exclude_key=("", "COMMENT"),
	)
	
# 	@. df_LDSS3[!, "TIMESTAMP"] = DateTime(
# 		parse(Date, df[!, "UT-DATE"]), parse(Time, df[!, "UT-TIME"])
# 	)
	
# 	@. df_LDSS3[!, "JD"] = datetime2julian(df.TIMESTAMP) - 2.456e6
			
# 	# Select header items to show
# 	cols = [
# 		"FILENAME",
# 		"TIMESTAMP",
# 		"JD",
# 		"OBJECT",
# 		"DISPERSR",
# 		"SLITMASK",
# 		"SPEED",
# 		"FILTER",
# 		"EXPTIME",
# 		"BINNING",
# 		"UT-DATE",
# 		"UT-TIME",
# 		"RA",
# 		"DEC",
# 		"EQUINOX",
# 		"EPOCH",
# 		"AIRMASS",
# 		"OBSERVER",
# 	]
	
# 	night_log = select(df, cols)
end

# ╔═╡ b5affee5-0322-43da-9f9a-05978fd90a21
function compute_cube(fpaths)	
	# Pre-sort chips into IMACS order
# 	chips = ["c1", "c6", "c2", "c5", "c3", "c8", "c4", "c7"]
# 	fpaths = [
# 		filter(s -> occursin(c, s), fpaths_glob)[1]
# 		for c in chips
# 	]
	
	# Load science frames and perform bias subtraction
	cube_list = []
	for (i, fpath) in enumerate(fpaths)
		ccd = CCDData(fpath)
		dxsec, dysec = ccd.hdr["DATASEC"] |> Meta.parse |> eval
		oxsec, oysec = ccd.hdr["BIASSEC"] |> Meta.parse |> eval
		# @show dxsec, dysec
		# @show oxsec, oysec
		#@show size(ccd[dysec, dxsec]), size(ccd[oysec, oxsec])
		push!(cube_list, ccd[dysec, dxsec] .- median(ccd[oysec, oxsec]))
	end
	
	return cube_list
end

# ╔═╡ c488270a-3126-4e38-a0c8-ee242115a3ea
#with_terminal() do
cube_LDSS3 = compute_cube(
	sort(glob("data/raw/LDSS3/ut150927/ccd0606c*.fits"))
)
#end

# ╔═╡ 71ba9181-90e4-4d12-97c0-462b3f1df077
let
	fig = Figure(resolution = (800, 600))
	step = 1
	chip_idx = 1
	hms = []
	axs = []
	for j in 1:2, i in 1
		ax = Axis(fig[i, j])
		if "$(chips[chip_idx])" in ["c5", "c6", "c7", "c8"]
			x_flip = x -> x
			y_flip = reverse
		else
			x_flip = reverse
			y_flip = x -> x
		end
		d = cube_LDSS3[chip_idx][x_flip(begin:step:end), y_flip(begin:step:end)]'
		nrows, ncols = size(d)
		hm = heatmap!(
			ax,
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
		push!(axs, ax)
		push!(hms, hm)
	end
	axs = reshape(axs, 1, 2)
	
	Colorbar(fig[:, end+1], hms[1], width=20, labelcolor=:black, label="Counts",)
	
	linkaxes!(axs...)
	hidedecorations!.(axs)
	
	path = "../../ACCESS_WASP-50b/figures/frames"
	mkpath(path)
	save("$(path)/sci_ldss3.png", fig)
	
	fig #|> as_svg
end

# ╔═╡ 4480ae72-3bb2-4e17-99be-28afc756332a
md"""
## Notebook setup
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
# ╟─7e111d10-0aa8-47bd-9ea8-186c1aecc321
# ╠═1e0db1d5-86b7-475c-897b-b0054575a5fa
# ╠═9257eef5-c18a-4226-9ff3-fc0ea90a1262
# ╠═c3def15d-d41a-4e3f-bb49-d57d46a6474e
# ╠═fbb58e2e-6dd4-42f5-acb1-037b2cc231c6
# ╠═3ca865bb-078a-4bd5-a632-a4b4b9d6e1e9
# ╠═c3ac2517-3bac-4686-9b1e-16a5c5b5865f
# ╠═2fdfc049-11b1-4972-9c48-e8724765126f
# ╟─90845d70-35d9-402d-8936-74936b069577
# ╠═86b4ace1-a891-41e5-a29f-f7eee5f8fb17
# ╠═3a6ab0c0-ba08-4151-9646-c19d45749b9f
# ╠═a0d54171-3d3b-472f-baca-2883924ec75e
# ╠═26feb668-4e7e-4a9d-a2b6-a5dac81e3ab7
# ╟─06a834f0-8c90-4013-af34-725166970969
# ╠═5c6e6f7b-70e0-49a8-b064-60dcf1440223
# ╠═b5affee5-0322-43da-9f9a-05978fd90a21
# ╠═c488270a-3126-4e38-a0c8-ee242115a3ea
# ╠═71ba9181-90e4-4d12-97c0-462b3f1df077
# ╟─4480ae72-3bb2-4e17-99be-28afc756332a
# ╠═3433ed02-c27c-4fe5-bfda-a5108a58407c
# ╟─6000db3d-0798-4f76-be31-617d43406b54
