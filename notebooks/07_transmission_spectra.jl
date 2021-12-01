### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ ef970c0c-d08a-4856-b10b-531bb5e7e53e
begin
	import Pkg
	Pkg.activate(Base.current_project())
	
	using AlgebraOfGraphics
	using CSV
	using CairoMakie
	using CCDReduction: fitscollection
	using Colors
	using DataFrames
	using DataFrameMacros
	using Dates
	using DelimitedFiles
	using Glob
	using ImageFiltering
	using Latexify
	using LaTeXStrings
	using Measurements
	using Measurements: value, uncertainty
	using NaturalSort
	using OrderedCollections
	using Printf
	using Statistics
	using StatsBase
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
			Lines = (linewidth=3, cycle=Cycle([:color, :linestyle], covary=true)),
			Scatter = (linewidth=10,),
			fontsize = 18,
			rowgap = 0,
			colgap = 0,
		)
	)
	
	# const COLORS = Makie.wong_colors()
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
end

# ╔═╡ 71672a5f-af6c-46f4-8e32-7bdd133ee039
using PhysicalConstants.CODATA2018: G

# ╔═╡ 7021cefd-f750-4422-b17b-c9abdc35dd2f
using UnitfulAstro

# ╔═╡ e8b8a0c9-0030-40f2-84e9-7fca3c5ef100
md"""
# Transmission Spectra

In this notebook we will load in the individual transmission spectra from each night, and combine them on a common wavelength basis.

$(TableOfContents(title="📖 Table of Contents"))
"""

# ╔═╡ 9413e640-22d9-4bfc-b4ea-f41c02a3bfde
md"""
## Load data
"""

# ╔═╡ 5100e6b4-03da-4e58-aad1-13376bcb4b59
md"""
First let's load up all of the data, including the white light transit depths from each night. This will be used for combining the transmission spectra later:
"""

# ╔═╡ c53be9cf-7722-4b43-928a-33e7b0463330
const DATA_DIR = "data/detrended/out_l/WASP50"

# ╔═╡ 1decb49e-a875-412c-938f-74b4fa0e2e85
maxmeasure(x, x_u, x_d) = x ± max(x_u, x_d)

# ╔═╡ 7b6d3a33-cb3b-4776-86f6-3af1663b9e49
dates_to_names = OrderedDict(
	"131219_IMACS" => "Transit 1 (IMACS)",
	"150927_IMACS" => "Transit 2 (IMACS)",
	"150927_LDSS3_flat" => "Transit 2 (LDSS3)",
	"150927_LDSS3_noflat" => "Transit 2 (LDSS3 noflat)",
	"161211_IMACS" => "Transit 3 (IMACS)",
 )

# ╔═╡ 1e8524c4-a732-4e2f-80a9-b5e7548ef2b2
function name(fpath, data_to_names)
	date_target = splitpath(split(glob(fpath)[1], "w50_")[2])[1]
	return dates_to_names[date_target]
end

# ╔═╡ 5c4fcb25-9a26-43f1-838b-338b33fb9ee6
begin
	cubes = OrderedDict{String, OrderedDict}()
	
	for dirpath in sort(glob("$(DATA_DIR)/w50_*"))
		# Read tspec file
		fpath = "$(dirpath)/transpec.csv"
		transit = name(fpath, dates_to_names)
		
		cubes[transit] = OrderedDict()
		
		df = cubes[transit]["tspec"] = CSV.read(fpath, DataFrame;
			normalizenames = true,
		)
		
		# Add wav bins for external instruments
		if occursin("LDSS3", dirpath)
			wbins = readdlm("$(dirpath)/wbins.dat", comments=true)
			cubes[transit]["tspec"][:, [:Wav_d, :Wav_u]] .= wbins
		end
		
		# For plotting later
		df.wav = mean([df.Wav_u, df.Wav_d])
		df.δ = maxmeasure.(df.Depth_ppm_, df.Depthup_ppm_, df.DepthDown_ppm_)
		
		fpath_WLC = "$(dirpath)/white-light/results.dat"
		df_WLC = CSV.File(
			fpath_WLC,
			comment = "#",
			normalizenames = true,
		) |> DataFrame
		
		symbol, p, p_u, p_d = eachcol(
			@subset(df_WLC, :Variable .== "p")
		)
	 	
		# Save WLC depths
		cubes[transit]["δ_WLC"] = maxmeasure(p[1], p_u[1], p_d[1])^2 * 1e6
	end
	
	cubes = sort(cubes)
end

# ╔═╡ e58ec082-d654-44e3-bcd4-906fc34171c8
md"""
## Combine spectra 🌈
"""

# ╔═╡ 11066667-9da2-4b36-b784-c3515c04a659
md"""
We start the combining process by saving the subset of the data sharing the same range of wavelength bins:
"""

# ╔═╡ cb1b277b-aa92-44de-91ce-88122bc34bb9
df_common = innerjoin(
	(cube["tspec"] for (transit, cube) in cubes)...,
	on = :wav,
	makeunique = true,
)

# ╔═╡ 029a9340-fd81-496c-a7e2-bbee60507710
cubes["Transit 1 (IMACS)"]["tspec"].wav

# ╔═╡ 461097e9-a687-4ef2-a5b4-8bf4d9e1c98f
antijoins = (
	"Transit 1 (IMACS)" => antijoin(
		cubes["Transit 1 (IMACS)"]["tspec"],
		cubes["Transit 2 (LDSS3)"]["tspec"],
		on = :wav,
	),
	
	"Transit 2 (IMACS)" => antijoin(
		cubes["Transit 2 (IMACS)"]["tspec"],
		cubes["Transit 2 (LDSS3)"]["tspec"],
		on = :wav,
	),
	
	"Transit 2 (LDSS3)" => antijoin(
		cubes["Transit 2 (LDSS3)"]["tspec"],
		cubes["Transit 1 (IMACS)"]["tspec"],
		on = :wav,
	),
	
	"Transit 3 (IMACS)" => antijoin(
		cubes["Transit 3 (IMACS)"]["tspec"],
		cubes["Transit 2 (LDSS3)"]["tspec"],
		on = :wav,
	),
)

# ╔═╡ 99757e41-0a88-4662-abef-0fad1bbbed1d
md"""
We next take these common transmission spectra and offset each, given by the difference between its corresponding white-light curve transit depth and the average WLC depth across nights:
"""

# ╔═╡ 84055852-1b9f-4221-95a7-ab48110bf78c
#depths_common = df_common[!, r"δ"] |>  x -> rename!(x, cubes.keys);

# ╔═╡ 2f377692-2abf-404e-99ea-a18c7af1a840
wlc_depths = [cube["δ_WLC"] for (transit, cube) in cubes]

# ╔═╡ 9141dba4-4c11-404d-b18a-b22f3466caba
Rₛ = 0.873u"Rsun"

# ╔═╡ 520d2cc3-00e0-46d8-83b2-5c740fd3bdd0
Mₚ = 1.78u"Mjup"

# ╔═╡ 45acc116-e585-4ddf-943d-128db7736921
function weightedmean2(m; corrected=true)
	if length(collect(m)) == 1
		return collect(m)[1] # Convert back to Measurement from skipmissing wrapper
	end
	x = value.(m)
	x_unc = uncertainty.(m)
	w = @. inv(x_unc^2)
	# Use ProbabilityWeights for bias correction
	a, b = mean_and_std(x, pweights(w), corrected=corrected)
	return a ± b
end

# ╔═╡ c405941d-bdcc-458f-b0bf-01abf02982e0
mean_wlc_depth = weightedmean2(wlc_depths)

# ╔═╡ 54c341d9-2065-48cf-89bd-11acf72bdf9d
Rₚ = √(mean_wlc_depth * 1e-6 * Rₛ^2)

# ╔═╡ cc3aec2c-6ca3-4817-9100-3e1c01df4651
Rₚ |> u"Rjup"

# ╔═╡ eaed62d7-5733-44b8-bd98-8b0fc4a18fe5
gₚ = G * Mₚ / Rₚ^2

# ╔═╡ 410644d5-e1e5-4107-aba7-e8a293bfff74
gₚ |> u"cm/s^2"

# ╔═╡ a915f236-8dae-4c91-8f96-fb9a805a0a7f
wlc_offsets = reshape(wlc_depths .- mean_wlc_depth, 1, :)

# ╔═╡ 4b9cfc02-5e18-422d-b18e-6301a659561a
begin
	depths_common = df_common[!, r"δ"] |>  x -> rename!(x, cubes.keys);
	depths_adj = depths_common #.- Measurements.value.(wlc_offsets)
	depths_adj.Combined = weightedmean2.(eachrow(Matrix(depths_adj)))
	insertcols!(depths_adj, 1,
		:Wav_d => df_common.Wav_d,
		:Wav_u => df_common.Wav_u,
		:Wav_cen => df_common.wav,
	)
end

# ╔═╡ 5d25caa3-916a-40b1-ba7c-ea1295afb775
md"""
Average precision per bin: $(round(Int, getproperty.(depths_adj[!, :Combined], :err) |> median)) ppm
"""

# ╔═╡ ed954843-34e5-49be-8643-e2671b659e06
depths_extra = let
	y = outerjoin((x[2] for x in antijoins)..., on=:wav, makeunique=true)
	z = y[!, r"δ"] |>  x -> rename!(x, cubes.keys)
	z.Combined = [weightedmean2(skipmissing(row)) for row ∈ eachrow(z)]
	insertcols!(z, 1,
		:Wav_cen => y.wav,
	)
	sort!(z, :Wav_cen)
end

# ╔═╡ 8c077881-fc5f-4fad-8497-1cb6106c6ed5
let
	fig = Figure(resolution=(1_000, 500))
		
	ax = Axis(
		fig[1, 1], xlabel="Wavelength (Å)", ylabel="Transit depth (ppm)",
		limits = (nothing, (15_500, 22_500)),
		grid = (linewidth=(0, 0),),
	)
	
	# Overplot lines
	vlines!(ax, [5892.9, 7682.0, 8189.0], color=:grey, linewidth=0.5)
	hlines!(ax, mean_wlc_depth.val, color=:grey, linewidth=3)
	hlines!.(ax, (mean_wlc_depth.val + mean_wlc_depth.err,
	mean_wlc_depth.val - mean_wlc_depth.err), linestyle=:dash, color=:grey)
	
	wav = depths_adj.Wav_cen
	ΔWLC_depth = mean_wlc_depth.err
	
	# Individual nights
	for (i, transit) in enumerate(keys(cubes))
		tspec = depths_adj[!, transit]
		depth = value.(tspec)
		depth_err = uncertainty.(tspec)

		CSV.write(
			"/home/mango/Desktop/tspec_$(CSV.normalizename(transit))w50.csv",
			DataFrame(
				:wav => wav,
				:depth => depth,
				:depth_err => depth_err,
			),
		)
		
		errorbars!(ax, wav, depth, depth_err;
			color = COLORS[i],
			linewidth = 3,
			whiskerwidth = 10,
		)
		scatter!(ax, wav, depth;
			color = COLORS[i],
			strokewidth = 0,
			markersize = 16,
			label = transit,
		)
	end
	
	# Combined
	tspec_combined = depths_adj.Combined
	depth_combined = value.(tspec_combined)
	depth_combined_err = uncertainty.(tspec_combined)

	CSV.write(
		"/home/mango/Desktop/tspec_combined_standard_w50.csv",
		DataFrame(
			:wav => wav,
			:depth => depth_combined,
			:depth_err => depth_combined_err,
		),
	)
	
	nudge = 25.0
	errorbars!(ax, wav .+ nudge, depth_combined, depth_combined_err;
		linewidth = 5,
		whiskerwidth = 10,
	)
	scatter!(ax, wav .+ nudge, depth_combined;
		color = :white,
		strokewidth = 3,
		markersize = 16,
		label = "Combined",
	)

	# Extra
	errorbars!(
		ax,
		depths_extra.Wav_cen .+ nudge,
		value.(depths_extra.Combined),
		uncertainty.(depths_extra.Combined);
		linewidth = 5,
		whiskerwidth = 10,
	)
	scatter!(ax, depths_extra.Wav_cen .+ nudge, value.(depths_extra.Combined);
		color = :grey,
		strokewidth = 3,
		markersize = 16,
		label = "Extra Combined",
	)
	
	# Write to file
	N = nrow(depths_adj)
	# CSV.write(
	# 	"/home/mango/Desktop/tspec_w50.csv",
	# 	DataFrame(
	# 	:Wlow => depths_adj.Wav_d,
	# 	:Wup => depths_adj.Wav_u,
	# 	:Depth => depth_combined,
	# 	:ErrUp => depth_combined_err,
	# 	:ErrLow => depth_combined_err,
	# 	:Instrument => fill("Magellan/IMACS", N),
	# 	:Offset => fill("NO", N)
	# 	),
	# )
	
	# Plot uncombined points
	for (i, (transit, df)) in enumerate(antijoins)
		wav = df.wav
		f = value.(df.δ)
		f_err = uncertainty.(df.δ)
		errorbars!(ax, wav, f, f_err;
			linewidth = 3,
			color = (COLORS[i], 0.5),
			whiskerwidth = 10,
		)
		scatter!(ax, wav, f;
			color = (COLORS[i], 0.5),
			markersize = 12,
		)
	end
	
	# All points
	# df_all = let
	
	# 	N = 25

	# 	df_blue = antijoins[1][2][1:4, :]

	# 	δ_5600 = weightedmean2(
	# 		(antijoins[1][2].δ[end], antijoins[2][2].δ[1], antijoins[4][2].δ[1])
	# 	)

	# 	df_middle = depths_adj

	# 	df_red = antijoins[3][2]

	# 	DataFrame(
	# 	"Wlow" => [df_blue.Wav_d..., 5600.0, df_middle.Wav_d..., df_red.Wav_d...],
	# 	"Wup" => [df_blue.Wav_u..., 5800.0, df_middle.Wav_u...,  df_red.Wav_u...],
	# 	"Depth" => [
	# 		value.(df_blue.δ)...,
	# 		δ_5600.val,
	# 		value.(df_middle.Combined)...,
	# 		value.(df_red.δ)...,
	# 	],
	# 	"ErrUp" => [
	# 		uncertainty.(df_blue.δ)...,
	# 		δ_5600.err,
	# 		uncertainty.(df_middle.Combined)...,
	# 		uncertainty.(df_red.δ)...,
	# 	],
	# 	"ErrLow" => [
	# 		uncertainty.(df_blue.δ)...,
	# 		δ_5600.err,
	# 		uncertainty.(df_middle.Combined)...,
	# 		uncertainty.(df_red.δ)...,
	# 	],
	# 	"Instrument" => fill("Magellan/IMACS", N),
	# 	"Offset?" => fill("NO", N)
	# )
	# end
	
	# CSV.write("/home/mango/Desktop/tspec_all.csv", df_all)
	
	#scatter!(ax, (df_all.Wlow .+ df_all.Wup) ./ 2, df_all.Depth, color=:red)

	axislegend(orientation=:horizontal, valign=:top, labelsize=16)
	
	path = "../../ACCESS_WASP-50b/figures/detrended"
	mkpath(path)
	save("$(path)/tspec.png", fig)
		
	fig #|> as_svg
end

# ╔═╡ 80b7b8a3-918f-4bc4-a59a-284cfbab5f8a
depths_extra

# ╔═╡ 07406bb6-7896-4cb9-83c0-73a92b5128c6
A = [
		3 ± 4 1 ± 2 2
		1     2     3
		4     5     6
		0     3     0
	]

# ╔═╡ 008c0cc0-ae58-4c05-a8b6-752cbe0d5eef
[weightedmean2(skipmissing(row)) for row in eachrow(A)]

# ╔═╡ 4dbc1714-f804-41d0-8cb6-0d15d0bc5cc5
mean_and_std(skipmissing([1, missing, 3]))

# ╔═╡ d85d48f6-7b40-4f94-b18d-8e6a2c350883
skipmissing([1, missing, 3])

# ╔═╡ 6e610bd5-8c4c-43a6-9a16-f3733d4a701a
# df_all = let
	
# 	N = 25
	
# 	df_blue = antijoins[1][2][1:4, :]
	
# 	δ_5600 = weightedmean2(
# 		(antijoins[1][2].δ[end], antijoins[2][2].δ[1], antijoins[4][2].δ[1])
# 	)
	
# 	df_middle = depths_adj
	
# 	df_red = antijoins[3][2]
	
# 	DataFrame(
# 	"Wlow" => [df_blue.Wav_d..., 5600.0, df_middle.Wav_d..., df_red.Wav_d...],
# 	"Wup" => [df_blue.Wav_u..., 5800.0, df_middle.Wav_u...,  df_red.Wav_u...],
# 	"Depth" => [
# 		value.(df_blue.δ)...,
# 		δ_5600.val,
# 		value.(df_middle.Combined)...,
# 		value.(df_red.δ)...,
# 	],
# 	"ErrUp" => [
# 		uncertainty.(df_blue.δ)...,
# 		δ_5600.err,
# 		uncertainty.(df_middle.Combined)...,
# 		uncertainty.(df_red.δ)...,
# 	],
# 	"ErrLow" => [
# 		uncertainty.(df_blue.δ)...,
# 		δ_5600.err,
# 		uncertainty.(df_middle.Combined)...,
# 		uncertainty.(df_red.δ)...,
# 	],
# 	"Instrument" => fill("Magellan/IMACS", N),
# 	"Offset?" => fill("NO", N)
# )
# end

# ╔═╡ f8a86915-f7d8-4462-980e-7b8124b13a3f
md"""
## Notebook setup
"""

# ╔═╡ 3510ead9-6e66-4fec-84ca-15c8a3ce4c3e
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
# ╟─e8b8a0c9-0030-40f2-84e9-7fca3c5ef100
# ╟─9413e640-22d9-4bfc-b4ea-f41c02a3bfde
# ╟─5100e6b4-03da-4e58-aad1-13376bcb4b59
# ╠═c53be9cf-7722-4b43-928a-33e7b0463330
# ╠═5c4fcb25-9a26-43f1-838b-338b33fb9ee6
# ╠═1decb49e-a875-412c-938f-74b4fa0e2e85
# ╠═1e8524c4-a732-4e2f-80a9-b5e7548ef2b2
# ╠═7b6d3a33-cb3b-4776-86f6-3af1663b9e49
# ╟─e58ec082-d654-44e3-bcd4-906fc34171c8
# ╟─11066667-9da2-4b36-b784-c3515c04a659
# ╠═cb1b277b-aa92-44de-91ce-88122bc34bb9
# ╠═029a9340-fd81-496c-a7e2-bbee60507710
# ╠═461097e9-a687-4ef2-a5b4-8bf4d9e1c98f
# ╟─99757e41-0a88-4662-abef-0fad1bbbed1d
# ╠═84055852-1b9f-4221-95a7-ab48110bf78c
# ╠═2f377692-2abf-404e-99ea-a18c7af1a840
# ╠═c405941d-bdcc-458f-b0bf-01abf02982e0
# ╠═9141dba4-4c11-404d-b18a-b22f3466caba
# ╠═54c341d9-2065-48cf-89bd-11acf72bdf9d
# ╠═cc3aec2c-6ca3-4817-9100-3e1c01df4651
# ╠═520d2cc3-00e0-46d8-83b2-5c740fd3bdd0
# ╠═eaed62d7-5733-44b8-bd98-8b0fc4a18fe5
# ╠═410644d5-e1e5-4107-aba7-e8a293bfff74
# ╠═71672a5f-af6c-46f4-8e32-7bdd133ee039
# ╠═7021cefd-f750-4422-b17b-c9abdc35dd2f
# ╠═a915f236-8dae-4c91-8f96-fb9a805a0a7f
# ╠═4b9cfc02-5e18-422d-b18e-6301a659561a
# ╠═45acc116-e585-4ddf-943d-128db7736921
# ╟─5d25caa3-916a-40b1-ba7c-ea1295afb775
# ╠═8c077881-fc5f-4fad-8497-1cb6106c6ed5
# ╠═ed954843-34e5-49be-8643-e2671b659e06
# ╠═80b7b8a3-918f-4bc4-a59a-284cfbab5f8a
# ╠═07406bb6-7896-4cb9-83c0-73a92b5128c6
# ╠═008c0cc0-ae58-4c05-a8b6-752cbe0d5eef
# ╠═4dbc1714-f804-41d0-8cb6-0d15d0bc5cc5
# ╠═d85d48f6-7b40-4f94-b18d-8e6a2c350883
# ╠═6e610bd5-8c4c-43a6-9a16-f3733d4a701a
# ╟─f8a86915-f7d8-4462-980e-7b8124b13a3f
# ╠═ef970c0c-d08a-4856-b10b-531bb5e7e53e
# ╠═3510ead9-6e66-4fec-84ca-15c8a3ce4c3e
