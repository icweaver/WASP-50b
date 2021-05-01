### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 1dc0c992-a894-11eb-2eab-f1fccba239f0
using PlutoUI, CairoMakie, Colors, CSV, DataFrames, DataFramesMeta, Glob, Measurements, Statistics

# ╔═╡ 74cdebd4-4d28-42df-93a7-176d794a932e
using Latexify

# ╔═╡ c53be9cf-7722-4b43-928a-33e7b0463330
const DATA_DIR = "data/detrended/out_l/WASP50"

# ╔═╡ 5c4fcb25-9a26-43f1-838b-338b33fb9ee6
cube = Dict(
	"tspec" => Dict(
		"Transit $i" => CSV.File(fpath) |> DataFrame
		for (i, fpath) in enumerate(sort(glob("$(DATA_DIR)/w50_*/transpec.csv")))
	),

	"BMA_WLC" => Dict(
		"Transit $i" => CSV.File(
			fpath,
			comment = "#",
			delim = ' ',			
			ignorerepeated = true,
		) |> DataFrame
		for (i, fpath) in enumerate(sort(glob(
			"$(DATA_DIR)/w50_*/white-light/results.dat"
		)))
	),
)

# ╔═╡ 855095aa-c7e1-4799-9fcb-070a95bf7656
begin
	const N_NIGHTS = length(cube["tspec"])
	const N_BINS = length(cube["tspec"]["Transit 1"][!, "Depth (ppm)"])
	#wavs = [] # Can be different lengths based on different binning schemes
	#tspecs = []
	tspecs = Matrix{AbstractFloat}(undef, N_BINS, N_NIGHTS+1)
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
			tspecs[:, begin] .= wav
			tspecs[:, j+1] .= depth .± max_depth_err
		else
			tspecs[IDX_START:end, begin] .= wav
			tspecs[IDX_START:end, j+1] .= depth .± max_depth_err
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
	tspecs_common = tspecs[IDX_START:end, begin+1:end]
	
	# Combine
	tspecs_wlc_adj = @. tspecs_common - wlc_offsets
	tspec_combined = mapslices(weightedmean, tspecs_wlc_adj, dims=2)
end

# ╔═╡ 8373a47f-9596-43f2-b5c8-bb179de3eec8
begin
	wav_common = tspecs[IDX_START:end, begin]
	
	df_wav = DataFrame([wav_common], [:Wavelength])
	
	df_tspecs = DataFrame(
		hcat(tspecs_wlc_adj, tspec_combined),
		[keys(cube["tspec"])..., "Combined"]
	)
	
	hcat(df_wav, df_tspecs)
end

# ╔═╡ 5c6d5b2d-3521-4626-b12d-d5685e0b09a5
# begin
# 	fig = Figure()
	
# 	ax = Axis(fig[1, 1], xlabel="Wavelength (Å)", ylabel="Relative depth")
	
# 	hlines!(ax, 0, color=:grey, linestyle=:dash)
	
# 	for (i, tspec) in enumerate(eachcol(tspecs_wlc_adj))
# 		depth = Measurements.value.(tspec)
# 		depth_err = Measurements.uncertainty.(tspec)
# 		errorbars!(ax, wav_common, depth .- mean_wlc_depth.val, depth_err;
# 			color = COLORS[i],
# 			whiskerwidth = 0,
# 			linewidth = 3,
# 		)
# 		scatter!(ax, wav_common, depth .- mean_wlc_depth.val;
# 			color = COLORS[i],
# 			strokewidth = 0,
# 			markersize = 12,
# 			label = "Transit $i",
# 		)
# 	end

# 	tspec_combined = mapslices(weightedmean, tspecs_common .- offsets, dims=2)
# 	f_combined = Measurements.value.(tspec_combined) |> vec
# 	f_combined_err = Measurements.uncertainty.(tspec_combined) |> vec

# 	wavelengths = cube["tspec"]["Transit 3"][!, "Wav_d"]
# 	errorbars!(ax, wavelengths, f_combined .- mean_wlc_depth.val, f_combined_err;
# 		linewidth = 3,
# 	)
# 	scatter!(ax, wavelengths, f_combined .- mean_wlc_depth.val;
# 		color = :white,
# 		strokewidth = 2,
# 		label = "Combined",
# 	)

# 	axislegend(orientation=:horizontal, valign=:bottom, labelsize=16)
		
# 	fig
# end

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

# ╔═╡ Cell order:
# ╠═1dc0c992-a894-11eb-2eab-f1fccba239f0
# ╠═c53be9cf-7722-4b43-928a-33e7b0463330
# ╠═5c4fcb25-9a26-43f1-838b-338b33fb9ee6
# ╠═855095aa-c7e1-4799-9fcb-070a95bf7656
# ╠═8373a47f-9596-43f2-b5c8-bb179de3eec8
# ╠═74cdebd4-4d28-42df-93a7-176d794a932e
# ╠═5c6d5b2d-3521-4626-b12d-d5685e0b09a5
# ╠═bef0918c-c645-4557-a2e5-00b6c26573bc
