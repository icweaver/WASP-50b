### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ e917bd8d-7f4a-44e4-9eb9-84199dd061f5
begin
	import Pkg
	Pkg.activate(Base.current_project())

	using PlutoUI
	using AlgebraOfGraphics
	using CairoMakie
	import CairoMakie.Makie.KernelDensity: kde
	using CSV, DataFrames, DataFramesMeta
	using DelimitedFiles
	using Glob
	using ImageFiltering
	using Latexify
	using Measurements, StatsBase, Unitful, UnitfulAstro
	import Measurements: value, uncertainty
	
	using OrderedCollections
	using Printf
	using NaturalSort
end

# ╔═╡ 086555bd-a886-4025-a354-0388df77a591


# ╔═╡ e8b8a0c9-0030-40f2-84e9-7fca3c5ef100
md"""
# Transmission Spectra

In this notebook we will load in the individual transmission spectra from each night, and combine them on a common wavelength basis.

$(TableOfContents(title="📖 Table of Contents"))
"""

# ╔═╡ bd80f202-b4bb-4682-aa62-e8245867208c
const FIG_PATH = "figures/detrended"

# ╔═╡ 0c752bd5-5232-4e82-b519-5ca23fff8a52
md"""
## Load data

First let's load up all of the data, including the white-light transit depths from each night.

!!! note "Data download"
	```
	rsync -azRP $H:/pool/sao_access/iweaver/GPTransmissionSpectra/./"out_*" detrended/ --exclude={"*george*","*mnest*"}

	rclone sync -P drive_ACCESS:papers/WASP-50b/data/detrended data/detrended
	```
"""

# ╔═╡ c53be9cf-7722-4b43-928a-33e7b0463330
@bind DATA_DIR Select(glob("data/detrended/out_*/WASP50"))

# ╔═╡ 1decb49e-a875-412c-938f-74b4fa0e2e85
maxmeasure(x, x_u, x_d) = x ± max(x_u, x_d)

# ╔═╡ b8abb240-65d6-4358-bd95-955be047371a
#fpath = "data/detrended/out_sp/WASP50/w50_131219_sp_IMACS/transpec.csv"

# ╔═╡ 774c4ab2-ad34-4a0d-a523-7234ac3d98e5
#fpath = "data/detrended/out_l/WASP50/w50_161211_IMACS/transpec.csv"

# ╔═╡ 7b6d3a33-cb3b-4776-86f6-3af1663b9e49
dates_to_names = OrderedDict(
	"131219_IMACS" => "Transit 1 (IMACS)",
	"150927_IMACS" => "Transit 2 (IMACS)",
	"150927_LDSS3" => "Transit 2 (LDSS3)",
	"161211_IMACS" => "Transit 3 (IMACS)",
	
	"131219_sp_IMACS" => "Transit 1 (IMACS)",
	"150927_sp_IMACS" => "Transit 2 (IMACS)",
	"150927_sp_LDSS3" => "Transit 2 (LDSS3)",
	"161211_sp_IMACS" => "Transit 3 (IMACS)",
 )

# ╔═╡ 1e8524c4-a732-4e2f-80a9-b5e7548ef2b2
function name(fpath, data_to_names)
	date_target = splitpath(split(fpath, "w50_")[2])[1]
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
		
		# Compute transmission spectra specific values
		df.wav = mean([df.Wav_u, df.Wav_d])
		df.δ = maxmeasure.(df.Depth_ppm_, df.Depthup_ppm_, df.DepthDown_ppm_)
		
		# Extract WLC information
		df_WLC = CSV.read(
			"$(dirpath)/white-light/results.dat",
			DataFrame,
			comment = "#",
			normalizenames = true,
		)
		symbol, p, p_u, p_d = eachcol(
			@subset(df_WLC, :Variable .== "p")
		)
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
df_common_0 = innerjoin(
	(cube["tspec"] for (transit, cube) in cubes)...,
	on = :wav,
	makeunique = true,
);

# ╔═╡ 59cc9c0a-332a-4f32-9f22-824f4c8f81b3
df_common_0

# ╔═╡ acde40fd-8ed4-4175-9a52-13ed91dc5495
md"""
Conversely, we also store which points in the spectrum are not common between all nights. `Transit 2 (LDSS3)` encompasses the spectra from all other nights, so we `antijoin` relative to this dataset:
"""

# ╔═╡ 461097e9-a687-4ef2-a5b4-8bf4d9e1c98f
dfs_unique = (
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
);

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

# ╔═╡ 4b9cfc02-5e18-422d-b18e-6301a659561a
begin
	df_common = @chain df_common_0 begin
		select(_, :Wav_d=>:Wlow, :Wav_u=>:Wup, :wav=>:Wcen, r"δ")
		# @rtransform :δ = :δ - wlc_offsets[1]
		# @rtransform :δ_1= :δ_1 - wlc_offsets[2]
		# @rtransform :δ_2 = :δ_2 - wlc_offsets[3]
		# @rtransform :δ_3 = :δ_3 - wlc_offsets[4]
		@rtransform :Combined = weightedmean2([:δ, :δ_1, :δ_2, :δ_3])
		rename(_, names(_, r"δ") .=> cubes.keys)
	end
	
	insertcols!(df_common, 4,
		:ΔW => df_common.Wup .- df_common.Wlow
	)
	
	df_common
end

# ╔═╡ ed954843-34e5-49be-8643-e2671b659e06
df_extra = let
	y = outerjoin((x[2] for x in dfs_unique)...;
		on = [:Wav_d, :Wav_u, :wav],
		makeunique = true,
	)
	z = y[!, r"δ"] |>  x -> rename!(x, cubes.keys)
	z.Combined = [weightedmean2(skipmissing(row)) for row ∈ eachrow(z)]
	insertcols!(z, 1,
		:Wlow => y.Wav_d,
		:Wup => y.Wav_u,
		:Wcen => y.wav,
		:ΔW => y.Wav_u .- y.Wav_d
	)
	sort!(z, :Wcen)
end

# ╔═╡ b32273bc-1bb5-406a-acfe-57fd643ded51
df_tspecs = sort(vcat(df_common, df_extra), :Wcen)
#df_tspecs = df_common

# ╔═╡ 64f608b9-76df-402e-801c-006dc3096f94
@with_terminal begin
	latextabular(df_tspecs, latex=false) |> print
end

# ╔═╡ 4ce88d88-ff08-4ea8-b397-fc61c22adc4e
df_tspecs[3, 5].err

# ╔═╡ fe04e248-c47b-4913-8405-26365c6027f4
avg_prec = round(Int, getproperty.(df_tspecs[!, :Combined], :err) |> median)

# ╔═╡ 5d25caa3-916a-40b1-ba7c-ea1295afb775
md"""
Average precision per bin: $(avg_prec) ppm
"""

# ╔═╡ 2f377692-2abf-404e-99ea-a18c7af1a840
wlc_depths = [cube["δ_WLC"] for (transit, cube) in cubes]

# ╔═╡ c405941d-bdcc-458f-b0bf-01abf02982e0
mean_wlc_depth = weightedmean2(wlc_depths)

# ╔═╡ a915f236-8dae-4c91-8f96-fb9a805a0a7f
wlc_offsets = reshape(wlc_depths .- mean_wlc_depth, 1, :)

# ╔═╡ 618aa031-1b09-4122-975b-562506beaae4
basename(dirname(DATA_DIR))

# ╔═╡ 09887c41-022a-4109-8c5d-0ba033c50bcb
function plot_tspec!(ax, df, col;
	nudge = 0.0,
	kwargs_errorbars = Dict(),
	kwargs_scatter = Dict(),
	color = :black,
	label = "enter label",
)
	wav, m = eachcol(select(df, "Wcen", col) |> dropmissing)
	f, ferr = value.(m), uncertainty.(m)

	errorbars!(ax, wav .+ nudge, f, ferr; color=color, kwargs_errorbars...)
	
	scatter!(ax, wav .+ nudge, f; color=color, kwargs_scatter..., label=label)
end

# ╔═╡ f92ecf4d-aab8-413e-a32d-5f77a969d1ca
md"""
## Offset test

Here we combine only the IMACS data and plot it along with the LDSS3 data set.
"""

# ╔═╡ 1a6067ca-645a-448b-815f-6a2966548ca6
begin
	df_IMACS = @chain df_tspecs begin
		select(_, ["Wlow","Wup", "Wcen","Transit 1 (IMACS)", "Transit 2 (IMACS)", "Transit 3 (IMACS)"])
		#z.Combined = [weightedmean2(skipmissing(row)) for row ∈ eachrow(z)]
		rename(_, names(_, r"Tr") .=> [:x1, :x2, :x3])
		#@rtransform :Combined = sum(skipmissing([:x1, :x2, :x3]))
		dropmissing(_, :x1)
		@rtransform :Combined = weightedmean2(skipmissing([:x1, :x2, :x3]))
		rename(_, names(_, r"x") .=> ["Transit 1 (IMACS)", "Transit 2 (IMACS)", "Transit 3 (IMACS)"])
	end
	insertcols!(df_IMACS, 4,
		:ΔW => df_IMACS.Wup .- df_IMACS.Wlow
	)
end

# ╔═╡ b555e372-2292-4f0e-b511-88b92588ad14
begin
	df_LDSS3 = @chain df_tspecs begin
		select(_, ["Wlow","Wup", "Wcen", "Transit 2 (LDSS3)"])
		#z.Combined = [weightedmean2(skipmissing(row)) for row ∈ eachrow(z)]
		rename(_, names(_, r"Tr") .=> [:x1])
		#@rtransform :Combined = sum(skipmissing([:x1, :x2, :x3]))
		dropmissing(_)
		@transform :Combined = :x1
		#@rtransform :Combined = weightedmean2(skipmissing([:x1, :x2, :x3]))
		rename(_, names(_, r"x") .=> ["Transit 2 (LDSS3)"])
	end
	insertcols!(df_LDSS3, 4,
		:ΔW => df_LDSS3.Wup .- df_LDSS3.Wlow
	)
end

# ╔═╡ 940ad41b-910c-40a8-8752-e68e13ff4a1f
@with_terminal begin
	latextabular(df_IMACS, latex=false) |> print
end

# ╔═╡ 094bd22d-8e42-440f-a78c-3a2787f380ea
df_LDSS3

# ╔═╡ f37d9e45-575c-40d9-8f26-31bd6cc6d145
avg_prec_IMACS = round(Int, getproperty.(df_IMACS[!, :Combined], :err) |> median)

# ╔═╡ b6fa6c00-14cf-47af-9593-c70514373db5
avg_prec_LDSS3 = round(Int, getproperty.(df_LDSS3[!, :Combined], :err) |> median)

# ╔═╡ 3af0d3b0-c698-4845-a74e-c7186b03a721
#CSV.write("$(DATA_DIR)/tspec_w50_IMACS.csv",
#	create_df(df_IMACS, instrument="Magellan/IMACS")
#)

# ╔═╡ 0fced6dc-fe96-4216-93b1-f8493f2402e9
#CSV.write("$(DATA_DIR)/tspec_w50_LDSS3.csv",
#	create_df(df_LDSS3, instrument="Clay/LDSS3")
#)

# ╔═╡ 146a2be7-1c08-4d7c-802f-41f65aeae0d5
md"""
## Retrieval params

Finally, we save the final combined transmission spectrum to file for our retrieval analysis, along with planet/star parameters computed from the WLC fits:
"""

# ╔═╡ 5718672b-1bc6-4676-8703-5fc06b83f0f9
# CSV.write("$(DATA_DIR)/tspec_w50_all.csv", create_df(df_tspecs))

# ╔═╡ b27f5a0a-812d-44c8-9c84-c74b0c58c794
# CSV.write("$(DATA_DIR)/tspec_w50.csv", create_df(df_common))

# ╔═╡ 9141dba4-4c11-404d-b18a-b22f3466caba
Rₛ = 0.873u"Rsun"

# ╔═╡ 54c341d9-2065-48cf-89bd-11acf72bdf9d
Rₚ = √(mean_wlc_depth * 1e-6 * Rₛ^2)

# ╔═╡ cc3aec2c-6ca3-4817-9100-3e1c01df4651
Rₚ |> u"Rjup"

# ╔═╡ 520d2cc3-00e0-46d8-83b2-5c740fd3bdd0
Mₚ = 1.78u"Mjup"

# ╔═╡ cb02a053-d048-43d9-950a-de3106019520
function create_df(df; instrument="add instrument")
	@select df begin
		:Wlow
		:Wup
		:Depth = value.(:Combined)
		:Errup = uncertainty.(:Combined)
		:ErrLow = uncertainty.(:Combined)
		:Instrument = instrument
		:Offset = "NO"
	end
end

# ╔═╡ f8a86915-f7d8-4462-980e-7b8124b13a3f
md"""
## Notebook setup
"""

# ╔═╡ e9b22d93-1994-4c31-a951-1ab00dc4c102
function savefig(fig, fpath)
	mkpath(dirname(fpath))
    save(fpath, fig)
	@info "Saved to: $(fpath)"
end

# ╔═╡ ef970c0c-d08a-4856-b10b-531bb5e7e53e
begin
	##############
	# PLOT CONFIGS
	##############
	const FIG_TALL = (900, 1_200)
	const FIG_WIDE = (800, 600)
	const FIG_LARGE = (1_200, 1_000)
	const COLORS_SERIES = to_colormap(:seaborn_colorblind, 9)
	const COLORS = parse.(Makie.Colors.Colorant,
		[
			"#a6cee3",  # Cyan
			"#fdbf6f",  # Yellow
			"#ff7f00",  # Orange
			"#1f78b4",  # Blue
		]
	)
	
	set_aog_theme!()
	update_theme!(
		Theme(
			Axis = (
				xlabelsize = 18,
				ylabelsize = 18,
				topspinevisible = true,
				rightspinevisible = true,
				topspinecolor = :darkgrey,
				rightspinecolor = :darkgrey
			),
			Label = (
				textsize = 18,
				padding = (0, 10, 0, 0),
				font = AlgebraOfGraphics.firasans("Medium")
			),
			Lines = (linewidth=3, cycle=Cycle([:color, :linestyle], covary=true)),
			Scatter = (linewidth=10,),
			Text = (font = AlgebraOfGraphics.firasans("Medium"),),
			palette = (color=COLORS, patchcolor=[(c, 0.35) for c in COLORS]),
			fontsize = 18,
			rowgap = 5,
			colgap = 5,
		)
	)
	
	COLORS

end

# ╔═╡ 8c077881-fc5f-4fad-8497-1cb6106c6ed5
let
	fig = Figure(resolution=(1_200, 400))
		
	ax = Axis(
		fig[1, 1], xlabel="Wavelength (Å)", ylabel="Transit depth (ppm)",
		limits = (4_600, 9_800, 15_500, 22_500),
		grid = (linewidth=(0, 0),),
	)
	
	# Overplot lines
	vlines!(ax, [5892.9, 7682.0, 8189.0], color=:grey, linestyle=:dash, linewidth=0.5)
	hlines!(ax, mean_wlc_depth.val, color=:grey, linewidth=3)
	hlines!.(ax, (mean_wlc_depth.val + mean_wlc_depth.err,
	mean_wlc_depth.val - mean_wlc_depth.err), linestyle=:dash, color=:grey)
	
	# Individual nights
	kwargs_errorbars = Dict(:whiskerwidth=>10.0, :linewidth=>1.0)
	kwargs_scatter = Dict(:markersize=>12.0)
	for (i, transit) in enumerate(keys(cubes))
		plot_tspec!(ax, df_tspecs, transit;
			kwargs_errorbars = kwargs_errorbars,
			kwargs_scatter = kwargs_scatter,
			color = COLORS[i],
			label = transit,
		)
	end
	
	# Combined
	nudge = 25.0
	kwargs_errorbars = Dict(:whiskerwidth=>10.0, :linewidth=>3.0)
	kwargs_scatter = Dict(:color=>:white, :strokewidth=>3.0, :markersize=>16.0)
	plot_tspec!(ax, df_tspecs, "Combined";
			nudge = nudge,
			kwargs_errorbars = kwargs_errorbars,
			kwargs_scatter = kwargs_scatter,
			label = "Combined",
	)

	text!(ax, "Average precision: $(avg_prec) ppm";
		position = (4700, 16000),
		align = (:left, :center),
	)
	
	axislegend(orientation=:horizontal, valign=:top, labelsize=16)
	
	suf = basename(dirname(DATA_DIR))
	savefig(fig, "$(FIG_PATH)/tspec_$(suf).png")
	
	fig
end

# ╔═╡ 72affb58-6f0c-4b76-9956-a027b57a0c8e
let
	fig = Figure(resolution=(1_200, 400))
		
	ax = Axis(
		fig[1, 1], xlabel="Wavelength (Å)", ylabel="Transit depth (ppm)",
		limits = (4_600, 9_800, 15_500, 22_500),
		grid = (linewidth=(0, 0),),
	)
	
	# Overplot lines
	vlines!(ax, [5892.9, 7682.0, 8189.0], color=:grey, linestyle=:dash, linewidth=0.5)
	hlines!(ax, mean_wlc_depth.val, color=:grey, linewidth=3)
	hlines!.(ax, (mean_wlc_depth.val + mean_wlc_depth.err,
	mean_wlc_depth.val - mean_wlc_depth.err), linestyle=:dash, color=:grey)
	
	# Individual nights
	kwargs_errorbars = Dict(:whiskerwidth=>10.0, :linewidth=>1.0)
	kwargs_scatter = Dict(:markersize=>12.0)
	for (i, transit) in enumerate(keys(cubes))
		if occursin("LDSS3", transit)
			plot_tspec!(ax, df_tspecs, transit;
				kwargs_errorbars = kwargs_errorbars,
				kwargs_scatter = kwargs_scatter,
				color = COLORS[i],
				label = transit,
			)
		else
			plot_tspec!(ax, df_IMACS, transit;
				kwargs_errorbars = kwargs_errorbars,
				kwargs_scatter = kwargs_scatter,
				color = COLORS[i],
				label = transit,
			)
		end
	end
	
	# Combined
	nudge = 25.0
	kwargs_errorbars = Dict(:whiskerwidth=>10.0, :linewidth=>3.0)
	kwargs_scatter = Dict(:color=>:white, :strokewidth=>3.0, :markersize=>16.0)
	plot_tspec!(ax, df_IMACS, "Combined";
			nudge = nudge,
			kwargs_errorbars = kwargs_errorbars,
			kwargs_scatter = kwargs_scatter,
			label = "Combined",
	)

	text!(ax, "Average precision (IMACS): $(avg_prec_IMACS) ppm";
		position = (4700, 16500),
		align = (:left, :center),
	)
	text!(ax, "Average precision (LDSS3): $(avg_prec_LDSS3)  ppm";
		position = (4700, 16000),
		align = (:left, :center),
	)
	
	axislegend(orientation=:horizontal, valign=:top, labelsize=16)

	suf = basename(dirname(DATA_DIR))
	savefig(fig, "$(FIG_PATH)/tspec_$(suf)_uncombined.png")
	
	fig
end

# ╔═╡ 6cb26c91-8a7c-4bd3-8978-d0c23105863c
const G = Unitful.G

# ╔═╡ eaed62d7-5733-44b8-bd98-8b0fc4a18fe5
gₚ = G * Mₚ / Rₚ^2 |> u"cm/s^2"

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
# ╠═086555bd-a886-4025-a354-0388df77a591
# ╟─e8b8a0c9-0030-40f2-84e9-7fca3c5ef100
# ╟─bd80f202-b4bb-4682-aa62-e8245867208c
# ╟─0c752bd5-5232-4e82-b519-5ca23fff8a52
# ╟─c53be9cf-7722-4b43-928a-33e7b0463330
# ╠═5c4fcb25-9a26-43f1-838b-338b33fb9ee6
# ╠═1decb49e-a875-412c-938f-74b4fa0e2e85
# ╠═1e8524c4-a732-4e2f-80a9-b5e7548ef2b2
# ╠═b8abb240-65d6-4358-bd95-955be047371a
# ╠═774c4ab2-ad34-4a0d-a523-7234ac3d98e5
# ╠═7b6d3a33-cb3b-4776-86f6-3af1663b9e49
# ╟─e58ec082-d654-44e3-bcd4-906fc34171c8
# ╟─11066667-9da2-4b36-b784-c3515c04a659
# ╠═cb1b277b-aa92-44de-91ce-88122bc34bb9
# ╠═59cc9c0a-332a-4f32-9f22-824f4c8f81b3
# ╟─acde40fd-8ed4-4175-9a52-13ed91dc5495
# ╠═461097e9-a687-4ef2-a5b4-8bf4d9e1c98f
# ╠═4b9cfc02-5e18-422d-b18e-6301a659561a
# ╠═45acc116-e585-4ddf-943d-128db7736921
# ╠═ed954843-34e5-49be-8643-e2671b659e06
# ╠═b32273bc-1bb5-406a-acfe-57fd643ded51
# ╠═64f608b9-76df-402e-801c-006dc3096f94
# ╠═4ce88d88-ff08-4ea8-b397-fc61c22adc4e
# ╟─5d25caa3-916a-40b1-ba7c-ea1295afb775
# ╠═fe04e248-c47b-4913-8405-26365c6027f4
# ╠═2f377692-2abf-404e-99ea-a18c7af1a840
# ╠═c405941d-bdcc-458f-b0bf-01abf02982e0
# ╠═a915f236-8dae-4c91-8f96-fb9a805a0a7f
# ╠═8c077881-fc5f-4fad-8497-1cb6106c6ed5
# ╠═618aa031-1b09-4122-975b-562506beaae4
# ╠═09887c41-022a-4109-8c5d-0ba033c50bcb
# ╟─f92ecf4d-aab8-413e-a32d-5f77a969d1ca
# ╠═1a6067ca-645a-448b-815f-6a2966548ca6
# ╠═b555e372-2292-4f0e-b511-88b92588ad14
# ╠═940ad41b-910c-40a8-8752-e68e13ff4a1f
# ╠═094bd22d-8e42-440f-a78c-3a2787f380ea
# ╠═f37d9e45-575c-40d9-8f26-31bd6cc6d145
# ╠═b6fa6c00-14cf-47af-9593-c70514373db5
# ╠═72affb58-6f0c-4b76-9956-a027b57a0c8e
# ╠═3af0d3b0-c698-4845-a74e-c7186b03a721
# ╠═0fced6dc-fe96-4216-93b1-f8493f2402e9
# ╟─146a2be7-1c08-4d7c-802f-41f65aeae0d5
# ╠═5718672b-1bc6-4676-8703-5fc06b83f0f9
# ╠═b27f5a0a-812d-44c8-9c84-c74b0c58c794
# ╠═9141dba4-4c11-404d-b18a-b22f3466caba
# ╠═54c341d9-2065-48cf-89bd-11acf72bdf9d
# ╠═cc3aec2c-6ca3-4817-9100-3e1c01df4651
# ╠═520d2cc3-00e0-46d8-83b2-5c740fd3bdd0
# ╠═eaed62d7-5733-44b8-bd98-8b0fc4a18fe5
# ╠═cb02a053-d048-43d9-950a-de3106019520
# ╟─f8a86915-f7d8-4462-980e-7b8124b13a3f
# ╟─e9b22d93-1994-4c31-a951-1ab00dc4c102
# ╠═ef970c0c-d08a-4856-b10b-531bb5e7e53e
# ╠═6cb26c91-8a7c-4bd3-8978-d0c23105863c
# ╠═e917bd8d-7f4a-44e4-9eb9-84199dd061f5
# ╟─3510ead9-6e66-4fec-84ca-15c8a3ce4c3e
