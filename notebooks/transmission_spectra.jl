### A Pluto.jl notebook ###
# v0.14.6

using Markdown
using InteractiveUtils

# ╔═╡ 1dc0c992-a894-11eb-2eab-f1fccba239f0
using CairoMakie, Colors, CSV, DataFrames, DataFramesMeta, Glob, Measurements, Statistics, Latexify

# ╔═╡ c53be9cf-7722-4b43-928a-33e7b0463330
const DATA_DIR = "data/detrended/out_l_C/WASP50"

# ╔═╡ 5c4fcb25-9a26-43f1-838b-338b33fb9ee6
cube = Dict(
	"tspec" => Dict(
		"Transit $i" => CSV.File(fpath) |> DataFrame
		for (i, fpath) in enumerate(sort(
			glob("$(DATA_DIR)/w50_*/transpec.csv")
		))
	),

	"BMA_WLC" => Dict(
		"Transit $i" => CSV.File(
			fpath,
			comment = "#",
			normalizenames = true,
		) |> DataFrame
		for (i, fpath) in enumerate(sort(glob(
			"$(DATA_DIR)/w50_*/white-light/results.dat"
		)))	
	),
)

# ╔═╡ d6918a50-f75f-47f5-86c6-e251f7ef1e12
cube["tspec"] |> keys

# ╔═╡ 855095aa-c7e1-4799-9fcb-070a95bf7656
begin
	const N_NIGHTS = length(cube["tspec"])
	const N_BINS = length(cube["tspec"]["Transit 1"][!, "Depth (ppm)"])
	#wavs = [] # Can be different lengths based on different binning schemes
	#tspecs = []
	wavs = zeros(Float64,  N_BINS, N_NIGHTS)
	tspecs = zeros(Measurement,  N_BINS, N_NIGHTS)
	ps = Matrix{Measurement}(undef, 1, N_NIGHTS) # Rₚ/Rₛ
	const IDX_START = 5 # Index where first common row between nights starts
	
	for (j, (transit, df)) in enumerate(cube["tspec"])
		# Unpack tspec data
		wav, depth, depth_u, depth_d = eachcol(
			df[!, ["Wav_d", "Depth (ppm)", "Depthup (ppm)", "DepthDown (ppm)"]]
		)
		# Adopt max uncertainty
		max_depth_err = maximum(hcat(depth_u, depth_d), dims=2) |> vec
		#push!(wavs, wav)
		#push!(tspecs, depth .± max_depth_err)
		if transit == "Transit 1"
			wavs[:, j] .= wav
			tspecs[:, j] .= depth .± max_depth_err
		else
			wavs[IDX_START:end, j] .= wav
			tspecs[IDX_START:end, j] .= depth .± max_depth_err
		end
		# Extract BMA p ≡ Rₚ/Rₛ for each night
		symbol, p, p_u, p_d = eachcol(
			@where(cube["BMA_WLC"][transit], :Variable .== "p")
		)
		ps[1, j] = p[1] ± maximum((p_u[1], p_d[1]))
	end
	
	# Compute offsets from mean BMA WLC depth
	wlc_depths = ps.^2 * 1e6
	mean_wlc_depth = weightedmean(wlc_depths)
	wlc_offsets = (wlc_depths .- mean_wlc_depth)
	
	# Extract common subset of transmission spectra
	tspecs_common = tspecs[IDX_START:end, :]
	
	# Combine
	tspecs_wlc_adj = @. tspecs_common - wlc_offsets
	tspec_combined = mapslices(weightedmean, tspecs_wlc_adj, dims=2)
end;

# ╔═╡ 8373a47f-9596-43f2-b5c8-bb179de3eec8
begin
	wav_common = wavs[IDX_START:end, end]
	df_wav = DataFrame([wav_common], [:Wavelength])
	
	df_tspecs = DataFrame(
		hcat(tspecs_wlc_adj, tspec_combined),
		[keys(cube["tspec"])..., "Combined"]
	)
	
	hcat(df_wav, df_tspecs) |> latexify
end

# ╔═╡ 64fd6b14-028f-4396-8a6d-240c73896174
Measurements.uncertainty.(tspec_combined) |> mean

# ╔═╡ 817aac12-13b9-44fe-a7e5-aa5ec7b8939e
mean_wlc_depth

# ╔═╡ bef0918c-c645-4557-a2e5-00b6c26573bc
const COLORS =  parse.(Colorant,
	[
		"#5daed9",  # Cyan
		"plum",
		"#f7ad4d",  # Yellow
		"mediumaquamarine",
		"#126399",  # Blue
		"#956cb4",  # Purple
		"#ff7f00",  # Orange
		"#029e73",  # Green
		"slategray",
	]
);

# ╔═╡ 5c6d5b2d-3521-4626-b12d-d5685e0b09a5
begin
	fig = Figure()
	
	ax = Axis(fig[1, 1], xlabel="Wavelength (Å)", ylabel="Relative depth")
	
	hlines!(ax, 0, color=:grey, linestyle=:dash)
	
	for (i, tspec) in enumerate(eachcol(tspecs_wlc_adj))
		depth = Measurements.value.(tspec)
		depth_err = Measurements.uncertainty.(tspec)
		errorbars!(ax, wav_common, depth .- mean_wlc_depth.val, depth_err;
			color = COLORS[i],
			whiskerwidth = 0,
			linewidth = 3,
		)
		scatter!(ax, wav_common, depth .- mean_wlc_depth.val;
			color = COLORS[i],
			strokewidth = 0,
			markersize = 12,
			label = "Transit $i",
		)
	end

	f_combined = Measurements.value.(tspec_combined) |> vec
	f_combined_err = Measurements.uncertainty.(tspec_combined) |> vec

	errorbars!(ax, wav_common, f_combined .- mean_wlc_depth.val, f_combined_err;
		linewidth = 3,
	)
	scatter!(ax, wav_common, f_combined .- mean_wlc_depth.val;
		color = :white,
		strokewidth = 2,
		label = "Combined",
	)
	
	# Plot blue points from Transit 1
	wav_blue = wavs[begin:IDX_START-1, begin] |> vec
	f_blue = Measurements.value.(tspecs[begin:IDX_START-1, begin]) |> vec
	f_blue_err = Measurements.uncertainty.(tspecs[begin:IDX_START-1, begin]) |> vec
	errorbars!(ax, wav_blue, f_blue .- mean_wlc_depth.val, f_blue_err)
	scatter!(ax, wav_blue, f_blue .- mean_wlc_depth.val, strokewidth=0)

	axislegend(orientation=:horizontal, valign=:bottom, labelsize=16)
		
	fig
end

# ╔═╡ e68c4b84-09eb-41e5-a688-36a83fa03625
import PlutoUI as pl

# ╔═╡ e8b8a0c9-0030-40f2-84e9-7fca3c5ef100
md"""
# Transmission spectra

In this notebook we will

$(pl.TableOfContents())
"""

# ╔═╡ Cell order:
# ╟─e8b8a0c9-0030-40f2-84e9-7fca3c5ef100
# ╠═c53be9cf-7722-4b43-928a-33e7b0463330
# ╠═d6918a50-f75f-47f5-86c6-e251f7ef1e12
# ╠═5c4fcb25-9a26-43f1-838b-338b33fb9ee6
# ╠═855095aa-c7e1-4799-9fcb-070a95bf7656
# ╠═8373a47f-9596-43f2-b5c8-bb179de3eec8
# ╠═64fd6b14-028f-4396-8a6d-240c73896174
# ╠═817aac12-13b9-44fe-a7e5-aa5ec7b8939e
# ╠═5c6d5b2d-3521-4626-b12d-d5685e0b09a5
# ╠═bef0918c-c645-4557-a2e5-00b6c26573bc
# ╠═e68c4b84-09eb-41e5-a688-36a83fa03625
# ╠═1dc0c992-a894-11eb-2eab-f1fccba239f0
