### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ ef970c0c-d08a-4856-b10b-531bb5e7e53e
begin
	import PlutoUI as Pl
	using CSV
	using CairoMakie
	using Colors
	using DataFrames
	using DataFramesMeta
	using DelimitedFiles
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
end

# ╔═╡ e8b8a0c9-0030-40f2-84e9-7fca3c5ef100
md"""
# Transmission spectra

In this notebook we will load in the individual transmission spectra from each night, and combine them on a common wavelength basis.

$(Pl.TableOfContents())
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
dates_to_names = Dict(
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
	cubes = Dict{String, Dict}()
	
	for dirpath in sort(glob("$(DATA_DIR)/w50_*"))
		# Read tspec file
		fpath = "$(dirpath)/transpec.csv"
		transit = name(fpath, dates_to_names)
		
		cubes[transit] = Dict()
		
		df = cubes[transit]["tspec"] = CSV.File(
			fpath,
			normalizenames = true,
		) |> DataFrame
		
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
			@where(df_WLC, :Variable .== "p")
		)
	 	
		cubes[transit]["δ_WLC"] = maxmeasure(p[1], p_u[1], p_d[1])^2 * 1e6
	end
	
	cubes = sort(cubes)
end

# ╔═╡ d6918a50-f75f-47f5-86c6-e251f7ef1e12
cubes.keys

# ╔═╡ e58ec082-d654-44e3-bcd4-906fc34171c8
md"""
## Combine spectra
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
depths_common = df_common[!, r"δ"] |>  x -> rename!(x, cubes.keys);

# ╔═╡ 2f377692-2abf-404e-99ea-a18c7af1a840
wlc_depths = [cube["δ_WLC"] for (transit, cube) in cubes]

# ╔═╡ c405941d-bdcc-458f-b0bf-01abf02982e0
mean_wlc_depth = mean(wlc_depths)

# ╔═╡ a915f236-8dae-4c91-8f96-fb9a805a0a7f
wlc_offsets = reshape(wlc_depths .- mean_wlc_depth, 1, :)

# ╔═╡ 4b9cfc02-5e18-422d-b18e-6301a659561a
begin
	depths_adj = depths_common .- Measurements.value.(wlc_offsets)
	#depths_adj = depths_common .- wlc_offsets
	depths_adj.Combined = weightedmean.(eachrow(depths_adj))
	insertcols!(depths_adj, 1,
		:Wav_d => df_common.Wav_d,
		:Wav_u => df_common.Wav_u,
		:Wav_cen => df_common.wav,
	)
end;

# ╔═╡ 5d25caa3-916a-40b1-ba7c-ea1295afb775
md"""
Average precision per bin: $(round(Int, getproperty.(depths_adj[!, :Combined], :err) |> median)) ppm
"""

# ╔═╡ bef0918c-c645-4557-a2e5-00b6c26573bc
begin
	const FIG_TALL = (900, 1_200)
	const FIG_WIDE = (1_350, 800)
	const COLORS = parse.(Colorant,
		[
			"#fdbf6f",  # Yellow
			"#a6cee3",  # Cyan
			"#1f78b4",  # Blue
			"#ff7f00",  # Orange
			"plum",
			"#956cb4",  # Purple
			"mediumaquamarine",
			"#029e73",  # Green
			"slategray",
		]
	)
end

# ╔═╡ 8c077881-fc5f-4fad-8497-1cb6106c6ed5
let
	fig = Figure(resolution=(1_000, 500))
		
	ax = Axis(
		fig[1, 1], xlabel="Wavelength (Å)", ylabel="Transit depth (ppm)",
		limits = (nothing, (mean_wlc_depth.val - 4_000, mean_wlc_depth.val + 4_000)),
		grid = (linewidth=(0, 0),),
	)
	
	vlines!(ax, [5892.9, 7682.0, 8189.0], color=:grey, linestyle=:dash)
	hlines!(ax, mean_wlc_depth.val, color=:grey, linestyle=:dash, linewidth=3)
	hlines!(ax, mean_wlc_depth.val + mean_wlc_depth.err, color=:grey, linewidth=3)
	hlines!(ax, mean_wlc_depth.val - mean_wlc_depth.err, color=:grey, linewidth=3)
	
	wav = depths_adj.Wav_cen
	ΔWLC_depth = mean_wlc_depth.err
	
	# Individual nights
	for (i, transit) in enumerate(keys(cubes))
		tspec = depths_adj[!, transit]
		depth = Measurements.value.(tspec)
		depth_err = Measurements.uncertainty.(tspec)
		
		errorbars!(ax, wav, depth, depth_err;
			color = COLORS[i],
			linewidth = 3,
		)
		scatter!(ax, wav, depth;
			color = COLORS[i],
			strokewidth = 0,
			markersize = 12,
			label = transit,
		)
	end
	
	# Combined
	tspec_combined = depths_adj.Combined
	depth_combined = Measurements.value.(tspec_combined)
	depth_combined_err = Measurements.uncertainty.(tspec_combined)

	errorbars!(ax, wav, depth_combined, depth_combined_err;
		linewidth = 5,
		whiskerwidth = 10,
	)
	
	scatter!(ax, wav, depth_combined;
		color = :white,
		strokewidth = 3,
		markersize = 12,
		label = "Combined",
	)
		
	# Plot uncombined points
	for (i, (transit, df)) in enumerate(antijoins)
		wav = df.wav
		f = Measurements.value.(df.δ)
		f_err = Measurements.uncertainty.(df.δ)
		errorbars!(ax, wav, f, f_err;
			linewidth = 3,
			color = (COLORS[i], 0.5),
		)

		scatter!(ax, wav, f;
			color = (COLORS[i], 0.5),
			markersize = 12,
		)
	end

	axislegend(orientation=:horizontal, valign=:bottom, labelsize=16)
		
	fig |> Pl.as_svg
end

# ╔═╡ Cell order:
# ╟─e8b8a0c9-0030-40f2-84e9-7fca3c5ef100
# ╟─9413e640-22d9-4bfc-b4ea-f41c02a3bfde
# ╟─5100e6b4-03da-4e58-aad1-13376bcb4b59
# ╠═c53be9cf-7722-4b43-928a-33e7b0463330
# ╠═d6918a50-f75f-47f5-86c6-e251f7ef1e12
# ╠═5c4fcb25-9a26-43f1-838b-338b33fb9ee6
# ╠═1decb49e-a875-412c-938f-74b4fa0e2e85
# ╠═1e8524c4-a732-4e2f-80a9-b5e7548ef2b2
# ╠═7b6d3a33-cb3b-4776-86f6-3af1663b9e49
# ╟─e58ec082-d654-44e3-bcd4-906fc34171c8
# ╟─11066667-9da2-4b36-b784-c3515c04a659
# ╠═cb1b277b-aa92-44de-91ce-88122bc34bb9
# ╠═461097e9-a687-4ef2-a5b4-8bf4d9e1c98f
# ╟─99757e41-0a88-4662-abef-0fad1bbbed1d
# ╠═84055852-1b9f-4221-95a7-ab48110bf78c
# ╠═2f377692-2abf-404e-99ea-a18c7af1a840
# ╠═c405941d-bdcc-458f-b0bf-01abf02982e0
# ╠═a915f236-8dae-4c91-8f96-fb9a805a0a7f
# ╠═4b9cfc02-5e18-422d-b18e-6301a659561a
# ╟─5d25caa3-916a-40b1-ba7c-ea1295afb775
# ╠═8c077881-fc5f-4fad-8497-1cb6106c6ed5
# ╠═bef0918c-c645-4557-a2e5-00b6c26573bc
# ╠═ef970c0c-d08a-4856-b10b-531bb5e7e53e
