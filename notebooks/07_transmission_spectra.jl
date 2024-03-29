### A Pluto.jl notebook ###
# v0.19.8

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
	import MarkdownLiteral: @mdx
	using AlgebraOfGraphics, CairoMakie
	using CSV, DataFrames, DataFramesMeta, DelimitedFiles, Glob, OrderedCollections
	using ImageFiltering, Measurements, StatsBase, Unitful, UnitfulAstro
	import Measurements: value, uncertainty
	using Latexify, Printf
	using Dates, NaturalSort
end

# ╔═╡ 209dae9c-4c14-4a02-bec9-36407bf1426f
begin
	const BASE_DIR = "data/detrended"
	const FIG_DIR = "figures/detrended"
	TableOfContents()
end

# ╔═╡ e8b8a0c9-0030-40f2-84e9-7fca3c5ef100
@mdx """
# Transmission Spectra

In this notebook we will load in the individual transmission spectra from each night, and combine them on a common wavelength basis.

!!! tip "Data download"
	```
	rclone sync -P ACCESS_box:WASP-50b/$(BASE_DIR) $(BASE_DIR)
	```
	* [Direct link](https://app.box.com/s/wr8tpof238cq8oj71ulaf69z9q0k7f9w)
"""

# ╔═╡ 0c752bd5-5232-4e82-b519-5ca23fff8a52
@mdx """
## Load data ⬇️

First let's load up all of the data, including the white-light transit depths from each night.
"""

# ╔═╡ c53be9cf-7722-4b43-928a-33e7b0463330
@bind DATA_DIR Select(glob("$(BASE_DIR)/out_*/WASP50"))

# ╔═╡ 288a19bf-77d1-43b8-b12d-05513190072a
fname_suff = (basename ∘ dirname)(DATA_DIR)

# ╔═╡ 96a24af1-1c91-45a9-a8f2-b4761f7f5cba
# df = cubes["Transit 1 (IMACS) sp"]["tspec"]

# ╔═╡ 1decb49e-a875-412c-938f-74b4fa0e2e85
maxmeasure(x, x_u, x_d) = x ± mean((x_u, x_d))

# ╔═╡ b8abb240-65d6-4358-bd95-955be047371a
#fpath = "data/detrended/out_sp/WASP50/w50_131219_sp_IMACS/transpec.csv"

# ╔═╡ 774c4ab2-ad34-4a0d-a523-7234ac3d98e5
#fpath = "data/detrended/out_l/WASP50/w50_161211_IMACS/transpec.csv"

# ╔═╡ 09d21338-91ef-4991-8e15-df33d720cb97
function tname(dirpath)
	if occursin("131219_IMACS", dirpath)
		transit = "Transit 1 (IMACS)"
	elseif occursin("131219_sp_IMACS", dirpath)
		transit = "Transit 1 (IMACS)"
	elseif occursin("150927_IMACS", dirpath)
		transit = "Transit 2 (IMACS)"
	elseif occursin("150927_sp_IMACS", dirpath)
		transit = "Transit 2 (IMACS)"
	elseif occursin("150927_LDSS3", dirpath)
		transit = "Transit 2 (LDSS3C)"
	elseif occursin("150927_sp_LDSS3", dirpath)
		transit = "Transit 2 (LDSS3C)"
	elseif occursin("161211_IMACS", dirpath)
		transit = "Transit 3 (IMACS)"
	elseif occursin("161211_sp_IMACS", dirpath)
		transit = "Transit 3 (IMACS)"
	end
	return transit
end

# ╔═╡ 5c4fcb25-9a26-43f1-838b-338b33fb9ee6
begin
	cubes = OrderedDict{String, OrderedDict}()

	for dirpath in sort(glob("$(DATA_DIR)/w50_*"))
		# Read tspec file
		fpath = "$(dirpath)/transpec.csv"
		transit = tname(dirpath)
		cubes[transit] = OrderedDict()

		df = cubes[transit]["tspec"] = CSV.read(fpath, DataFrame;
			normalizenames = true,
		)

		# Add wbins for LDSS3C
		if occursin("_sp_LDSS3", fpath)
			df.Wav_d, df.Wav_u = eachcol(readdlm("data/detrended/wbins/w50_bins_species.dat"; comments=true))
		elseif occursin("_LDSS3", fpath)
			df.Wav_d, df.Wav_u = eachcol(readdlm("data/detrended/wbins/w50_bins_LDSS3.dat"; comments=true))
		end

		# Compute transmission spectra specific values
		df.wav = mean([df.Wav_u, df.Wav_d])
		df.δ = maxmeasure.(df.Depth_ppm_, df.Depthup_ppm_, df.DepthDown_ppm_)

		# Extract WLC information
		df_WLC = CSV.read("$(dirpath)/white-light/results.dat", DataFrame;
			comment = "#",
			normalizenames = true,
			stripwhitespace = true,
		)
		symbol, p, p_u, p_d = eachcol(
			@subset(df_WLC, :Variable .== "p")
		)
	 	cubes[transit]["δ_WLC"] = maxmeasure(p[1], p_u[1], p_d[1])^2 * 1e6
	end

	delete!(cubes["Transit 1 (IMACS)"]["tspec"], 3)
	cubes = sort(cubes)
end

# ╔═╡ e58ec082-d654-44e3-bcd4-906fc34171c8
@mdx """
## Combined spectra 🌈

We next step through the process to produce the following figure:
"""

# ╔═╡ 4b8497f7-d285-4632-98e3-0699d284f291
label_spec!(f, label) = Label(
	f, label, valign=:top, padding=(0, 0, 0, 10), tellwidth=false, tellheight=false)

# ╔═╡ 11066667-9da2-4b36-b784-c3515c04a659
@mdx """
We start the combining process by saving the subset of the data sharing the same range of wavelength bins:
"""

# ╔═╡ cb1b277b-aa92-44de-91ce-88122bc34bb9
df_common_0 = innerjoin(
	(cube["tspec"] for (transit, cube) in cubes)...,
	on = :wav,
	makeunique = true,
);

# ╔═╡ acde40fd-8ed4-4175-9a52-13ed91dc5495
@mdx """
Conversely, we also store which points in the spectrum are not common between all nights. `Transit 2 (LDSS3C)` encompasses the spectra from all other nights, so we `antijoin` relative to this dataset:
"""

# ╔═╡ 461097e9-a687-4ef2-a5b4-8bf4d9e1c98f
dfs_unique = (
	"Transit 1 (IMACS)" => antijoin(
		cubes["Transit 1 (IMACS)"]["tspec"],
		cubes["Transit 2 (LDSS3C)"]["tspec"],
		on = :wav,
	),

	"Transit 2 (IMACS)" => antijoin(
		cubes["Transit 2 (IMACS)"]["tspec"],
		cubes["Transit 2 (LDSS3C)"]["tspec"],
		on = :wav,
	),

	"Transit 2 (LDSS3C)" => antijoin(
		cubes["Transit 2 (LDSS3C)"]["tspec"],
		cubes["Transit 1 (IMACS)"]["tspec"],
		on = :wav,
	),

	"Transit 3 (IMACS)" => antijoin(
		cubes["Transit 3 (IMACS)"]["tspec"],
		cubes["Transit 2 (LDSS3C)"]["tspec"],
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
	# Use AnalyticWeights for bias correction
	a, b = mean_and_std(x, aweights(w); corrected)
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

# ╔═╡ a8345769-504e-446b-a679-88ef7913fee5
# This is the old version we used to use (gives super tiny errorbars)
# https://stackoverflow.com/a/2415343/16402912
# It's also biased (as stated in the SO comments section and wiki link)
# https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance
# This is why we use the more sophisticated weighted average, `weightedmean2`
function weightedmean3(m)
	if length(collect(m)) == 1
		return collect(m)[1] # Convert back to Measurement from skipmissing wrapper
	end
	x = value.(m)
	x_unc = uncertainty.(m)
	w = @. 1.0 / x_unc^2
	# Use standard Weights for bias correction
	a, b = mean_and_std(x, weights(w))
	return a ± b
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
# df_tspecs = df_common

# ╔═╡ 64f608b9-76df-402e-801c-006dc3096f94
latextabular(df_tspecs, latex=false) |> PlutoUI.Text

# ╔═╡ fe04e248-c47b-4913-8405-26365c6027f4
avg_prec = getproperty.(df_tspecs[!, :Combined], :err) |> median

# ╔═╡ 5d25caa3-916a-40b1-ba7c-ea1295afb775
@mdx """
Average precision per bin: $(avg_prec) ppm
"""

# ╔═╡ 2f377692-2abf-404e-99ea-a18c7af1a840
wlc_depths = [cube["δ_WLC"] for (transit, cube) in cubes]

# ╔═╡ c405941d-bdcc-458f-b0bf-01abf02982e0
mean_wlc_depth = weightedmean2(wlc_depths)

# ╔═╡ a915f236-8dae-4c91-8f96-fb9a805a0a7f
wlc_offsets = reshape(wlc_depths .- mean_wlc_depth, 1, :)

# ╔═╡ e61c7fbd-d030-4915-bbe4-d8f4405e9c3f
@views begin
	wbins = df_tspecs[:, [:Wlow, :Wup]]
	wbins_odd = wbins[begin:2:end, :]
	wbins_even = wbins[begin+1:2:end, :]
end;

# ╔═╡ 09887c41-022a-4109-8c5d-0ba033c50bcb
function plot_tspec!(ax, df, col;
	nudge = 0.0,
	kwargs_errorbars = (),
	kwargs_scatter = (),
	color = :black,
	label = "enter label",
)
	wav, m = eachcol(select(df, "Wcen", col) |> dropmissing)
	f, ferr = value.(m), uncertainty.(m)

	errorbars!(ax, wav .+ nudge, f, ferr; color=color, kwargs_errorbars...)

	scatter!(ax, wav .+ nudge, f; color=color, kwargs_scatter..., label=label)
end

# ╔═╡ f92ecf4d-aab8-413e-a32d-5f77a969d1ca
@mdx """
## Offset test 🧪

Here we combine only the IMACS data and plot it along with the LDSS3C data set.
"""

# ╔═╡ 1a6067ca-645a-448b-815f-6a2966548ca6
begin
	df_IMACS = @chain df_tspecs begin
		select(
			_,
			["Wlow", "Wup", "Wcen",
			"Transit 1 (IMACS)", "Transit 2 (IMACS)", "Transit 3 (IMACS)"]
		)
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
		select(_, ["Wlow","Wup", "Wcen", "Transit 2 (LDSS3C)"])
		rename(_, names(_, r"Tr") .=> [:x1])
		dropmissing(_)
		@transform :Combined = :x1
		rename(_, names(_, r"x") .=> ["Transit 2 (LDSS3C)"])
	end
	insertcols!(df_LDSS3, 4,
		:ΔW => df_LDSS3.Wup .- df_LDSS3.Wlow
	)
end

# ╔═╡ 940ad41b-910c-40a8-8752-e68e13ff4a1f
latextabular(df_IMACS, latex=false) |> PlutoUI.Text

# ╔═╡ f37d9e45-575c-40d9-8f26-31bd6cc6d145
avg_prec_IMACS = getproperty.(df_IMACS[!, :Combined], :err) |> median

# ╔═╡ b6fa6c00-14cf-47af-9593-c70514373db5
avg_prec_LDSS3 = getproperty.(df_LDSS3[!, :Combined], :err) |> median

# ╔═╡ 92029dbc-2633-414f-a8ee-0ec61eda0313
df_common

# ╔═╡ 3af0d3b0-c698-4845-a74e-c7186b03a721
let
    #f = "$(DATA_DIR)/tspec_w50_IMACS.csv"
	#CSV.write(f,
	#	create_df(df_IMACS; instrument="Magellan/IMACS")
	#)
	#@info "Saved to $(f)"
	#f = "$(DATA_DIR)/tspec_w50_LDSS3.csv"
	#CSV.write(f,
	#	create_df(df_LDSS3; instrument="Clay/LDSS3C")
	#)
	#@info "Saved to $(f)"
end

# ╔═╡ 27811c9d-1ee5-49ca-bf09-04dc75dd66be
@mdx """
## A closer look at Transit 2 🔎
"""

# ╔═╡ b971503b-bfef-4ca7-98ad-b5940bbda10f
@views begin
wbins_LDSS3 = df_LDSS3[:, [:Wlow, :Wup]]
wbins_LDSS3_odd = wbins_LDSS3[begin:2:end, :]
end;

# ╔═╡ 943844ce-a78b-403d-8bae-341216308392
df_IMACS |> dropmissing#[!, "Transit 2 (IMACS)"]

# ╔═╡ ef04759d-6a2d-488b-ae3b-6595c35dd70a
wbins_LDSS3

# ╔═╡ 146a2be7-1c08-4d7c-802f-41f65aeae0d5
@mdx """
## Retrieval params 🐕

Finally, we save the final combined transmission spectrum to file for our retrieval analysis, along with planet/star parameters computed from the WLC fits:
"""

# ╔═╡ 4faac7de-8c38-4f2c-be85-569e8fc83d28
df_tspecs

# ╔═╡ 3846da6a-af41-47a0-9318-76757a1dba15
df_common

# ╔═╡ 9141dba4-4c11-404d-b18a-b22f3466caba
Rₛ = 0.873u"Rsun"

# ╔═╡ 54c341d9-2065-48cf-89bd-11acf72bdf9d
Rₚ = √(mean_wlc_depth * 1e-6 * Rₛ^2)

# ╔═╡ cc3aec2c-6ca3-4817-9100-3e1c01df4651
Rₚ |> u"Rjup"

# ╔═╡ 520d2cc3-00e0-46d8-83b2-5c740fd3bdd0
Mₚ = 1.78u"Mjup"

# ╔═╡ cb02a053-d048-43d9-950a-de3106019520
function create_df(df; instrument="add_instrument")
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

# ╔═╡ 5718672b-1bc6-4676-8703-5fc06b83f0f9
let
    f = "$(DATA_DIR)/tspec_w50_3.csv"
	CSV.write(f, create_df(df_tspecs; instrument="IMACS+LDSS3C"))
	@info "Saved to $(f)"
	#f = "$(DATA_DIR)/tspec_w50.csv"
	#CSV.write(f, create_df(df_common; instrument="IMACS+LDSS3C"))
	#@info "Saved to $(f)"
end

# ╔═╡ f8a86915-f7d8-4462-980e-7b8124b13a3f
@mdx """
## Notebook setup 🔧
"""

# ╔═╡ e9b22d93-1994-4c31-a951-1ab00dc4c102
function savefig(fig, fpath)
	mkpath(dirname(fpath))
    save(fpath, fig, pt_per_unit=1)
	@info "Saved to: $(fpath)"
end

# ╔═╡ 8644fa54-0407-4494-aef4-eb497a86c35d
let
	fig = Figure(resolution=(1100, 400))

	ax = Axis(
		fig[1, 1], xlabel="Wavelength (Å)", ylabel="Transit depth (ppm)",
		limits = (4_700, 9_769, 15_500, 22_500),
		grid = (linewidth=(0, 0),),
	)

	# Overplot lines
	#vspan!(ax, wbins_LDSS3_odd[:, 1], wbins_LDSS3_odd[:, 2], color=(:darkgrey, 0.25))
	vspan!(ax, wbins_odd[:, 1], wbins_odd[:, 2], color=(:darkgrey, 0.25))
	#vspan!(ax, wbins_even[:, 1], wbins_even[:, 2], color=(:black, 0.25))
	vlines!(ax, [5892.9, 7682.0, 8189.0], color=:grey, linestyle=:dash, linewidth=0.5)
	hlines!(ax, mean_wlc_depth.val, color=:grey, linewidth=3)
	hlines!.(ax, (mean_wlc_depth.val + mean_wlc_depth.err,
	mean_wlc_depth.val - mean_wlc_depth.err), linestyle=:dash, color=:grey)

	# Transit 2
	transit = "Transit 2 (IMACS)"
	plot_tspec!(ax, df_IMACS, transit;
		nudge = -25.0,
		kwargs_errorbars = (whiskerwidth=10.0, linewidth=3.0),
		kwargs_scatter = (color=Cycled(2), strokewidth=3.0, markersize=16.0),
		label = transit,
	)
	transit = "Transit 2 (LDSS3C)"
	plot_tspec!(ax, df_LDSS3, transit;
		nudge = 25.0,
		kwargs_errorbars = (whiskerwidth=10.0, linewidth=3.0),
		kwargs_scatter = (color=Cycled(3), strokewidth=3.0, markersize=16.0),
		label = transit,
	)

	# text!(ax, "Average precision (IMACS): $(round(avg_prec_IMACS, digits=2)) ppm";
	# 	position = (4700, 16500),
	# 	align = (:left, :center),
	# )
	# text!(ax, "Average precision (LDSS3C): $(round(avg_prec_LDSS3, digits=2))  ppm";
	# 	position = (4700, 16000),
	# 	align = (:left, :center),
	# )

	axislegend(orientation=:horizontal, valign=:top, labelsize=16)
	
	suf = basename(dirname(DATA_DIR))
	savefig(fig, "$(FIG_DIR)/tspec_transit_2_$(suf).pdf")
	
	fig
end

# ╔═╡ ef970c0c-d08a-4856-b10b-531bb5e7e53e
begin
	##############
	# PLOT CONFIGS
	##############
	const FIG_TALL = 72 .* (6, 8)
	const FIG_WIDE = 72 .* (12, 6)
	const FIG_LARGE = 72 .* (12, 12)
	const COLORS_SERIES = categorical_colors(:seaborn_colorblind, 9)
	const COLORS = parse.(Makie.Colors.Colorant,
		[
			"#66C2A5",  # Green
			"#FDBF6F",  # Yellow
			"#FF7F00",  # Orange
			"#1F78B4",  # Blue
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
				rightspinecolor = :darkgrey,
			),
			Label = (
				textsize = 18,
				font = AlgebraOfGraphics.firasans("Medium"),
			),
			Lines = (linewidth=3,),
			Scatter = (linewidth=10,),
			Text = (font = AlgebraOfGraphics.firasans("Regular"), textsize=18),
			palette = (color=COLORS, patchcolor=[(c, 0.35) for c in COLORS]),
			figure_padding = (0, 1.5, 0, 0),
			fontsize = 18,
			rowgap = 5,
			colgap = 5,
		)
	)

	COLORS
end

# ╔═╡ 8c077881-fc5f-4fad-8497-1cb6106c6ed5
function plot_tspecs!(ax; i=1, nudge_indiv=-50, nudge_comb=50, wbins=wbins_odd)
	# ax = Axis(fig[1, i];
	# 	xlabel="Wavelength (Å)", ylabel="Transit depth (ppm)",
	# 	limits,
	# 	grid = (linewidth=(0, 0),),
	# )
	
	# Individual nights
	for (i, transit) in enumerate(keys(cubes))
		plot_tspec!(ax, df_tspecs, transit;
			nudge = nudge_indiv,
			kwargs_errorbars = (whiskerwidth=10.0, linewidth=1.0),
			kwargs_scatter = (markersize=12.0,),
			color = COLORS[i],
			label = transit,
		)
	end
	
	# Combined
	plot_tspec!(ax, df_tspecs, "Combined";
			nudge = nudge_comb,
			kwargs_errorbars = (whiskerwidth=10.0, linewidth=3.0),
			kwargs_scatter = (color=:white, strokewidth=3.0, markersize=16.0),
			label = "Combined",
	)

	text!(ax, "Average precision: $(round(Int, avg_prec)) ppm";
		position = (4750, 16000),
		align = (:left, :center),
	)
end

# ╔═╡ cafc773a-ee51-4bd2-b766-182ad728e253
begin
	fig = Figure(resolution=(1100, 400))
	if occursin("sp", DATA_DIR)
		ax_NaD = Axis(fig[1, 1];
			limits = (5_755, 6_031, 15_500, 22_500),
			xticks = LinearTicks(4),
		)
		ax_K = Axis(fig[1, 2];
			limits = (7_652, 7_712, 15_500, 22_500),
			xticks = LinearTicks(4),
		)
		ax_Na8200 = Axis(fig[1, 3];
			limits = (8_069, 8_309, 15_500, 22_500),
			xticks = LinearTicks(4),
		)

		# Overplot lines
		for ax ∈ (ax_NaD, ax_K, ax_Na8200)
			#vspan!(ax, wbins_odd[:, 1], wbins_odd[:, 2], color=(:darkgrey, 0.25))
			vspan!(ax, wbins_odd[:, 1], wbins_odd[:, 2], color=(:darkgrey, 0.25))
			#vspan!(ax, wbins_even[:, 1], wbins_even[:, 2], color=(:black, 0.25))
			vlines!(ax, [5892.9, 7682.0, 8189.0], color=:grey, linestyle=:dash, linewidth=0.5)
			hlines!(ax, mean_wlc_depth.val, color=:grey, linewidth=3)
			hlines!.(ax, (mean_wlc_depth.val + mean_wlc_depth.err,
			mean_wlc_depth.val - mean_wlc_depth.err), linestyle=:dash, color=:grey)
		end

		elems = [
			[MarkerElement(marker=:circle, markersize=20, color=COLORS[i])
				for i in 1:4];
			MarkerElement(
				marker=:circle, markersize=20, color=:white, strokewidth=3,)
		]
		elem_labels = [[df.first for df ∈ dfs_unique]; "Combined"]
		
		plot_tspecs!(ax_NaD, nudge_indiv=-11.5, nudge_comb=11.5)
		label_spec!(fig[1, 1], "Na I -D")
		
		plot_tspecs!(ax_K; i=2, nudge_indiv=-2.5, nudge_comb=2.5, wbins=wbins_even)
		label_spec!(fig[1, 2], "K I")
		
		plot_tspecs!(ax_Na8200; i=2, nudge_indiv=-10, nudge_comb=10)
		label_spec!(fig[1, 3], "Na I-8200")
		
		hideydecorations!.((ax_K, ax_Na8200))

		Label(fig[2, 2], "Wavelength (Å)", tellwidth=false, tellheight=true)
		Label(fig[1, 0], "Transit depth (ppm)", tellwidth=true, tellheight=false,
		rotation=π/2)

		Legend(fig[0, 3], elems, elem_labels;
			tellwidth = false,
			tellheight = true,
			orientation = :horizontal,
		)
		colgap!(fig.layout, 30)
	else
		ax = Axis(fig[1, 1];
			limits = (4_700, 9_769, 15_500, 22_500),
			xlabel = "Wavelength (Å)",
			ylabel = "Transit depth (ppm)",
		)
		# Overplot lines
		#vspan!(ax, wbins_odd[:, 1], wbins_odd[:, 2], color=(:darkgrey, 0.25))
		vspan!(ax, wbins_odd[:, 1], wbins_odd[:, 2], color=(:darkgrey, 0.25))
		#vspan!(ax, wbins_even[:, 1], wbins_even[:, 2], color=(:black, 0.25))
		vlines!(ax, [5892.9, 7682.0, 8189.0], color=:grey, linestyle=:dash, linewidth=0.5)
		hlines!(ax, mean_wlc_depth.val, color=:grey, linewidth=3)
		hlines!.(ax, (mean_wlc_depth.val + mean_wlc_depth.err,
		mean_wlc_depth.val - mean_wlc_depth.err), linestyle=:dash, color=:grey)
		plot_tspecs!(ax)
		axislegend(orientation=:horizontal, valign=:top, labelsize=16)
	end
	
	savefig(fig, "$(FIG_DIR)/tspec_$(fname_suff).pdf")
	save("/home/mango/Desktop/tspec_wasp50b.svg", fig; pt_per_unit=1)
	fig
end

# ╔═╡ 72affb58-6f0c-4b76-9956-a027b57a0c8e
let
	fig = Figure(resolution=(1100, 400))

	ax = Axis(
		fig[1, 1], xlabel="Wavelength (Å)", ylabel="Transit depth (ppm)",
		limits = (4_700, 9_769, 15_500, 22_500),
		grid = (linewidth=(0, 0),),
	)

	# Overplot lines
	#vspan!(ax, wbins_odd[:, 1], wbins_odd[:, 2], color=(:darkgrey, 0.25))
	vspan!(ax, wbins_odd[:, 1], wbins_odd[:, 2], color=(:darkgrey, 0.25))
	#vspan!(ax, wbins_even[:, 1], wbins_even[:, 2], color=(:black, 0.25))
	vlines!(ax, [5892.9, 7682.0, 8189.0], color=:grey, linestyle=:dash, linewidth=0.5)
	hlines!(ax, mean_wlc_depth.val, color=:grey, linewidth=3)
	hlines!.(ax, (mean_wlc_depth.val + mean_wlc_depth.err,
	mean_wlc_depth.val - mean_wlc_depth.err), linestyle=:dash, color=:grey)

	# Individual nights
	kwargs_errorbars = (whiskerwidth=10.0, linewidth=1.0)
	kwargs_scatter = (markersize=12.0,)
	for (i, transit) in enumerate(keys(cubes))
		if occursin("LDSS3C", transit)
			p = plot_tspec!(ax, df_tspecs, transit;
				nudge = 25.0,
				kwargs_errorbars = (
					whiskerwidth=10.0, linewidth=3.0, color=:black
				),
				kwargs_scatter = (strokewidth=3.0, markersize=16.0),
				color = COLORS[i],
				label = transit,
			)
			translate!(p, 0, 0, 10)
		else
			plot_tspec!(ax, df_IMACS, transit;
				nudge = -50.0,
				kwargs_errorbars,
				kwargs_scatter,
				color = COLORS[i],
				label = transit,
			)
		end
	end

	# Combined
	p = plot_tspec!(ax, df_IMACS, "Combined";
			nudge = 50.0,
			kwargs_errorbars = (whiskerwidth=10.0, linewidth=3.0),
			kwargs_scatter = (color=:white, strokewidth=3.0, markersize=16.0),
			label = "Combined (IMACS)",
	)
	translate!(p, 0, 0, 100)

	text!(ax, "Average precision (IMACS): $(round(avg_prec_IMACS, digits=2)) ppm";
		position = (4700, 16500),
		align = (:left, :center),
	)
	text!(ax, "Average precision (LDSS3C): $(round(avg_prec_LDSS3, digits=2))  ppm";
		position = (4700, 16000),
		align = (:left, :center),
	)

	axislegend(orientation=:horizontal, valign=:top, labelsize=16)

	# suf = basename(dirname(DATA_DIR))
	# savefig(fig, "$(FIG_DIR)/tspec_$(suf)_uncombined.png")
	
	fig
end

# ╔═╡ 6cb26c91-8a7c-4bd3-8978-d0c23105863c
const G = Unitful.G

# ╔═╡ eaed62d7-5733-44b8-bd98-8b0fc4a18fe5
gₚ = G * Mₚ / Rₚ^2 |> u"cm/s^2"

# ╔═╡ Cell order:
# ╟─e8b8a0c9-0030-40f2-84e9-7fca3c5ef100
# ╠═209dae9c-4c14-4a02-bec9-36407bf1426f
# ╟─0c752bd5-5232-4e82-b519-5ca23fff8a52
# ╟─c53be9cf-7722-4b43-928a-33e7b0463330
# ╠═5c4fcb25-9a26-43f1-838b-338b33fb9ee6
# ╠═288a19bf-77d1-43b8-b12d-05513190072a
# ╠═96a24af1-1c91-45a9-a8f2-b4761f7f5cba
# ╠═1decb49e-a875-412c-938f-74b4fa0e2e85
# ╠═b8abb240-65d6-4358-bd95-955be047371a
# ╠═774c4ab2-ad34-4a0d-a523-7234ac3d98e5
# ╠═09d21338-91ef-4991-8e15-df33d720cb97
# ╟─e58ec082-d654-44e3-bcd4-906fc34171c8
# ╠═cafc773a-ee51-4bd2-b766-182ad728e253
# ╠═4b8497f7-d285-4632-98e3-0699d284f291
# ╠═8c077881-fc5f-4fad-8497-1cb6106c6ed5
# ╟─11066667-9da2-4b36-b784-c3515c04a659
# ╠═cb1b277b-aa92-44de-91ce-88122bc34bb9
# ╟─acde40fd-8ed4-4175-9a52-13ed91dc5495
# ╠═461097e9-a687-4ef2-a5b4-8bf4d9e1c98f
# ╠═4b9cfc02-5e18-422d-b18e-6301a659561a
# ╠═45acc116-e585-4ddf-943d-128db7736921
# ╠═a8345769-504e-446b-a679-88ef7913fee5
# ╠═ed954843-34e5-49be-8643-e2671b659e06
# ╠═b32273bc-1bb5-406a-acfe-57fd643ded51
# ╠═64f608b9-76df-402e-801c-006dc3096f94
# ╟─5d25caa3-916a-40b1-ba7c-ea1295afb775
# ╠═fe04e248-c47b-4913-8405-26365c6027f4
# ╠═2f377692-2abf-404e-99ea-a18c7af1a840
# ╠═c405941d-bdcc-458f-b0bf-01abf02982e0
# ╠═a915f236-8dae-4c91-8f96-fb9a805a0a7f
# ╠═e61c7fbd-d030-4915-bbe4-d8f4405e9c3f
# ╠═09887c41-022a-4109-8c5d-0ba033c50bcb
# ╟─f92ecf4d-aab8-413e-a32d-5f77a969d1ca
# ╠═72affb58-6f0c-4b76-9956-a027b57a0c8e
# ╠═1a6067ca-645a-448b-815f-6a2966548ca6
# ╠═b555e372-2292-4f0e-b511-88b92588ad14
# ╠═940ad41b-910c-40a8-8752-e68e13ff4a1f
# ╠═f37d9e45-575c-40d9-8f26-31bd6cc6d145
# ╠═b6fa6c00-14cf-47af-9593-c70514373db5
# ╠═92029dbc-2633-414f-a8ee-0ec61eda0313
# ╠═3af0d3b0-c698-4845-a74e-c7186b03a721
# ╟─27811c9d-1ee5-49ca-bf09-04dc75dd66be
# ╠═8644fa54-0407-4494-aef4-eb497a86c35d
# ╠═b971503b-bfef-4ca7-98ad-b5940bbda10f
# ╠═943844ce-a78b-403d-8bae-341216308392
# ╠═ef04759d-6a2d-488b-ae3b-6595c35dd70a
# ╟─146a2be7-1c08-4d7c-802f-41f65aeae0d5
# ╠═4faac7de-8c38-4f2c-be85-569e8fc83d28
# ╠═3846da6a-af41-47a0-9318-76757a1dba15
# ╠═5718672b-1bc6-4676-8703-5fc06b83f0f9
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
