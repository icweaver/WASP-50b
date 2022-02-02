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

# ╔═╡ 5939e2a3-8407-4579-8274-e3891acbabb1
begin
	import Pkg
	Pkg.activate(Base.current_project())

	using PlutoUI
	import MarkdownLiteral: @mdx
	using AlgebraOfGraphics, CairoMakie
	import CairoMakie.Makie.KernelDensity: kde
	using CSV, DataFrames, DataFramesMeta, DelimitedFiles, Glob, OrderedCollections
	using ImageFiltering, Measurements, Statistics
	using Latexify, Printf
	using Dates
	using CondaPkg
	CondaPkg.add("numpy"); CondaPkg.resolve()
	using PythonCall
end

# ╔═╡ f3a8f6fb-023c-4077-8c73-7502e56eb607
using StatsBase

# ╔═╡ 25d1284c-7260-4f3a-916a-b2814d2606af
begin
	const BASE_DIR = "data/detrended"
	const FIG_DIR = "figures/detrended"
	TableOfContents()
end

# ╔═╡ 506eeeb2-e56d-436b-91b8-605e52201563
@mdx """
# Detrended white-light curves

In this notebook we will visualize the detrended white-light curves from IMACS and LDSS3C. We used the average orbital and system parameters obtained from these detrended fits to place uniform constraints on our binned wavelength analysis.

!!! tip "Data download"
	```
	rclone sync -P ACCESS_box:WASP-50b/$(BASE_DIR) $(BASE_DIR)
	```
	* [Direct link](https://app.box.com/s/wr8tpof238cq8oj71ulaf69z9q0k7f9w)
"""

# ╔═╡ 782806a6-efd2-45a9-b898-788a276c282b
@mdx """
## Load data ⬇

First, let's load the relevant data needed for this notebook:
"""

# ╔═╡ a8d91138-73e7-4382-a032-37daec54a9c0
@bind DATA_DIR Select(glob("$(BASE_DIR)/out_*/WASP50"))

# ╔═╡ e873f9c6-fd1a-4227-9df1-70c626e4a0a1
function name(fpath, dates_to_names)
	date_instr = splitpath(split(glob(fpath)[1], "w50_")[2])[1]
	return dates_to_names[date_instr]
end

# ╔═╡ e72dba55-6a33-462f-aeac-5f62b25cb46a
dates_to_names = Dict(
	"131219_IMACS" => "Transit 1 (IMACS)",
	"131219_sp_IMACS" => "Transit 1 (IMACS)",
	"150927_IMACS" => "Transit 2 (IMACS)",
	"150927_sp_IMACS" => "Transit 2 (IMACS)",
	"150927_LDSS3C" => "Transit 2 (LDSS3C)",
	"150927_sp_LDSS3C" => "Transit 2 (LDSS3C)",
	"161211_IMACS" => "Transit 3 (IMACS)",
	"161211_sp_IMACS" => "Transit 3 (IMACS)",
 )

# ╔═╡ 579e62da-7ffb-4639-bd73-3826ade1cfa2
@mdx """
The data cube is organized by night as follows:

```julia
cubes
├── Transit 1 (IMACS)
├── Transit 2 (IMACS)
├── Transit 2 (LDSS3C)
└── Transit 3 (IMACS)
    ├── samples
    │   └── Dict
    ├── models
    │   └── Dict
    └── results
        ├── Dict
        └── Table
```

where `samples` is a dictionary of the Bayesian model averaged posterior samples for each fitted parameter, `models` is a dictionary of the final detrended light curve and associated Gaussian Process parameters, and `results` is a table summarizing the mean and +/- 1σ uncertainty associated with each sample.
"""

# ╔═╡ a8cf11e2-796e-45ff-bdc9-e273b927700e
@mdx """
## Transit curves 🌅
"""

# ╔═╡ ae82d3c1-3912-4a5e-85f5-6383af42291e
@mdx """
Plotting the data from the `models` cube returns the following detrended white-light curves:
"""

# ╔═╡ 0dd63eaf-1afd-4caf-a74b-7cd217b3c515
# Returns value `v` from `Variables` columns in results.dat file
val(df, v) = @subset(df, :Variable .== v)[1, "Value"]

# ╔═╡ d43ec3eb-1d5e-4a63-b5e8-8dcbeb57ae7c
# Computes orbital phase
function ϕ(t, t₀, P)
	phase = ((t - t₀) / P) % 1.0
    phase ≥ 0.5 && (phase -= 1.0)
	return phase
end

# ╔═╡ b28bb1b6-c148-41c4-9f94-0833e365cad4
@mdx """
## Summary table 📋
"""

# ╔═╡ 30b84501-fdcd-4d83-b929-ff354de69a17
@mdx """
We summarize the Bayesian Model Averag (BMA) results for selected parameters for each night below, and average together each parameter from each night, weighted by its maximum uncertainty per night:
"""

# ╔═╡ 694a6067-390e-4b79-b8b1-87e0c6d15f48
.1382 - 0.13746

# ╔═╡ 22c72eeb-8e32-4d7c-86c8-ab117735769e
@mdx """
The standard version we used to use gives relatively errorbars (thanks for bringing this to my attention during my TAC, Dave!) and is also particularly biased for small sample sizes. For these reasons, we have opted for the more sophisticated machinery of [reliability weighting](https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Reliability_weights) going forward.

* Old version:
```python
# From https://stackoverflow.com/a/2415343/16402912
def weightedmean(values, weights):
	# weights = 1./err_bar^2. Where err_bar=std & err_bar^2 = variance
	average = np.average(values, weights=weights)
	# Fast and numerically precise
	variance = np.average((values-average)**2, weights=weights)
	return average, np.sqrt(variance)
```

AKA
```julia
function weightedmean3(m)
	x = value.(m)
	x_unc = uncertainty.(m)
	w = @. 1.0 / x_unc^2
	a, b = mean_and_std(x, weights(w))
	return a ± b
end
```

* New version: `weightedmean2`
"""

# ╔═╡ a96ee19b-39fe-4f5b-bc75-124b4e422713
begin
	@pyexec """
	global weightedmean, np
	import numpy as np
	def weightedmean(vals, uncs):
		vals = np.array(vals)
		uncs = np.array(uncs)
		weights = 1.0 / uncs**2
		average = np.average(vals, weights=weights)
		variance = np.average((vals-average)**2, weights=weights)
		return average, np.sqrt(variance)
	"""
	weightedmean_py(vals, uncs) = @pyeval("weightedmean")(vals, uncs)
end

# ╔═╡ 47278372-b311-4ea7-bfa4-82b8f95c97fa
import Measurements: value, uncertainty

# ╔═╡ bbc14e57-57fe-4811-91d4-d07b102cfa5d
@doc raw"""
Given a collection of $N$ dependent observations $\boldsymbol x = [x_1, x_2, \dots, x_N]$, the biased-corrected weighted mean estimator $\hat x \equiv \mu^* \pm s_\mathrm{w}$ is given by:

```math
\begin{align}
\mu^* &= \sum w_i x_i / \sum w_i\,, \\
s_\mathrm{w} &= \sqrt{
    \frac
    {\sum w_i(x_i - \mu^*)^2}
    {\sum w_i - \sum w_i^2 / \sum w_i}
}\,,
\end{align}
```

where $\sum$ is taken to be the sum over all indices $i$ for convenience, and $w_i$ is defined to be the inverse variance $w_i \equiv 1/\sigma_i^2$ for each measurement $x_i$. More [here](https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Reliability_weights).
"""
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

# ╔═╡ 7cbadc33-5b32-4816-8092-09054c64073f
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

# ╔═╡ c936b76a-636b-4f10-b556-fa19808c1562
@mdx """
### Print table
"""

# ╔═╡ 68ec4343-5f6c-4dfd-90b5-6393b4c819b9
@mdx """
## Corner plots 📐
"""

# ╔═╡ e452d2b1-1010-4ce3-8d32-9e9f1d0dfa0b
@mdx """
To visualize the spread of each parameter, we define the dimensionless metric:

```math
Δx \\equiv \\frac{x - \\overline{\\text{BMA}}_\\mu}{\\overline{\\text{BMA}}_\\sigma}\\quad,
```

where ``x`` is a sample from the posterior distribution for the given parameter, ``\\overline{\\text{BMA}}_\\mu`` is the BMA averaged across nights, and ``\\overline{\\text{BMA}}_\\sigma`` is the maximum uncertainty, also averaged across nights. These values correspond to the `Combined` column in the table above. ``\\Delta x`` is then a measure of the displacement of each sample in its posterior distribution from its corresponding average BMA value, scaled by the uncertainty in the average BMA.
"""

# ╔═╡ 706f1fb6-2895-48f6-a315-842fbf35da18
function scale_samples(samples, param, BMA)
	m = @subset(BMA, :Parameter .== param)[!, "Combined"][1]
	Δx = (samples .- m.val) ./ m.err
	return Δx
end

# ╔═╡ 7fdd0cb7-7af3-4ca1-8939-9d9b7d6e9527
function plot_x!(ax, x; h=1, errorbar_kwargs=(), scatter_kwargs=())
	errorbars!(ax, [x.val], [h], [x.err], [x.err], direction=:x; errorbar_kwargs...)
	scatter!(ax, [x.val], [h]; scatter_kwargs...)
end

# ╔═╡ 3225ad37-6264-48b8-b1fc-c6f0136f8beb
get_m_scaled(m, param, BMA) = scale_samples(m, param, BMA)

# ╔═╡ 8f795cda-b905-4c89-ade5-05e35b10d4b3
get_m(param, transit, BMA) = BMA[BMA.Parameter .== param, transit][1]

# ╔═╡ 266b5871-22e5-4b4a-80b1-468a9ddd9193
function get_x(param, transit, BMA)
	m = get_m(param, transit, BMA)
	return get_m_scaled(m, param, BMA)
end

# ╔═╡ 8c9428a9-f2ad-4329-97a3-97ffa6b40f28
function plot_scaled_density!(ax, dist)
	k = kde(dist)
	k_scaled = k.density ./ maximum(k.density)
	lines!(ax, k.x, k_scaled)
	band!(ax, k.x, zero(k.x), k_scaled)
end

# ╔═╡ ed935d16-ddce-4334-a880-005732b38936
# Params to show in corner plot
const PARAMS = OrderedDict(
	"p" => "Rp/Rs",
	"t0" => "t₀",
	"P" => "P",
	"rho" => "ρs",
	"aR" => "a/Rs",
	"inc" => "i",
	"b" => "b",
	"q1" => "u",
);

# ╔═╡ 6fcd1377-8364-45a3-9ff6-89d61df1ef42
# Number of levels `n` to show in contour plots
compute_levels(A, n) = reverse(
	range(maximum(A), step=-maximum(A)/(n+1), length=(n+1))
)

# ╔═╡ 2cbc6ddb-210e-41e8-b745-5c41eba4e778
function plot_corner!(fig, samples, params; n_levels=4, color=:blue)
	for (j, p1) in enumerate(params), (i, p2) in enumerate(params)
		# 1D plot
		if i == j
			plot_scaled_density!(fig[i, i], samples[p1])
			# density!(fig[i, i], samples[p1];
			# 	color = (color, 0.125),
			# 	strokewidth = 3,
			# 	strokecolor = color,
			# 	strokearound = true,
			# 	offset = 1.0 - maximum(k.density),
			# )
		end
		# 2D plot
		if i > j
			Z = kde((samples[p1], samples[p2]), npoints=(2^5, 2^5),)
			contourf!(fig[i, j], Z;
				levels = compute_levels(Z.density, n_levels),
				colormap = cgrad(
					range(Makie.Colors.colorant"white", color, length=n_levels),
					alpha=0.5,
				),
			)
			contour!(fig[i, j], Z;
				levels = compute_levels(Z.density, n_levels),
				color = color,
				linewidth = 3,
			)
		end
	end
end

# ╔═╡ 807e913f-d8c3-41e9-acb3-4c024dedd67b
@mdx """
## A closer look at Transit 2
"""

# ╔═╡ 2bfaa81e-f0bd-4527-a5fe-091589b1f76e


# ╔═╡ 30ae3744-0e7e-4c16-b91a-91eb518fba5b
@mdx """
## Notebook setup
"""

# ╔═╡ bf9c0b95-fe17-425d-8904-8298f7e5451c
function savefig(fig, fpath)
	mkpath(dirname(fpath))
    save(fpath, fig)
	@info "Saved to: $(fpath)"
end

# ╔═╡ db539901-f0b0-4692-a8d2-6c72dff41196
begin
	@pyexec """
	global np, pickle, load_npz, load_pickle
	import numpy as np
	import pickle

	def load_npz(fpath, allow_pickle=False):
		return np.load(fpath, allow_pickle=allow_pickle)[()]

	def load_pickle(fpath):
		with open(fpath, "rb") as f:
			data = pickle.load(f)
		return data
	"""
	load_npz(s; allow_pickle=false) = @pyeval("load_npz")(s, allow_pickle=allow_pickle)
	load_pickle(s) = @pyeval("load_pickle")(s)
end;

# ╔═╡ f539e06d-a1b5-413a-b90e-91cb0bbd5a4c
load_data(fpath_sample, fpath_model, fpath_result; allow_pickle=true) = Dict(
	"samples" => load_pickle(fpath_sample),
	"models" => load_npz(fpath_model, allow_pickle=allow_pickle),
	"results" => CSV.File(
		fpath_result,
		comment = "#",
		normalizenames = true,
	) |> DataFrame
)

# ╔═╡ 2191791b-df62-4f1b-88bf-060cc47896b2
cubes = OrderedDict(
	name(fpath_sample, dates_to_names) => load_data(
		fpath_sample, fpath_model, fpath_result
	)
	for (fpath_sample, fpath_model, fpath_result) ∈ zip(
		sort(glob("$(DATA_DIR)/w50*/white-light/BMA_posteriors.pkl")),
		sort(glob("$(DATA_DIR)/w50*/white-light/BMA_WLC.npy")),
		sort(glob("$(DATA_DIR)/w50*/white-light/results.dat")),
	)
)

# ╔═╡ 88a92472-f20d-4908-8bfd-fcbf11f80b8a
cubes["Transit 1 (IMACS)"]["models"]["LC_det_err"]

# ╔═╡ ee9347b2-e97d-4f66-9c21-7487ca2c2e30
begin
	summary_tables = DataFrame[]
	
	for (transit, cube) in cubes
		param_idxs = [
			findfirst(cube["results"][!, :Variable] .== param)
			for param in keys(PARAMS)
		]
		summary_table = cube["results"][param_idxs, :]
		push!(summary_tables, summary_table)
	end
end

# ╔═╡ de0a4468-56aa-4748-80a0-6c9ab6b8579e
BMA_matrix = let
	x = hcat((
	summary[!, "Value"] .± maximum((summary[!, "SigmaUp"], summary[!, "SigmaDown"]))
	for summary in summary_tables
	)...) |> x-> hcat(x, weightedmean2.(eachrow(x)))
end

# ╔═╡ 19fcaa15-6f01-46a6-8225-4b5cafd89cc1
BMA = DataFrame(
	[PARAMS.vals BMA_matrix],
	["Parameter", keys(cubes)..., "Combined"]
);

# ╔═╡ c7a179a3-9966-452d-b430-a28b2f004bc5
latextabular(BMA, latex=false) |> PlutoUI.Text

# ╔═╡ d279e93e-8665-41b2-bd5c-723458fabe86
BMA |> x -> latexify(x, env=:table) |> PlutoUI.Text

# ╔═╡ 18c90ed4-2f07-493f-95b2-e308cd7a03a9
let
	fig  = Figure()
	ax = Axis(fig[1, 1])

	x = get_x("Rp/Rs", "Transit 2 (IMACS)", BMA)
	plot_x!(ax, x; scatter_kwargs=(color=:black, markersize=20))

	x = get_x("Rp/Rs", "Transit 2 (LDSS3C)", BMA)
	plot_x!(ax, x; h=3, scatter_kwargs=(color=:red, markersize=20))
		
	fig
end

# ╔═╡ 8969ff2e-d160-4eb3-8f77-f8c2006518ff
get_x("Rp/Rs", "Transit 2 (IMACS)", BMA)

# ╔═╡ 25dd0c88-089b-406b-ac0f-6f21c57fe986
@with_terminal begin
	map(enumerate(BMA_matrix[:, end])) do (i, x)
		#println(x.val)
		println(PARAMS.vals[i], ": ", round(x.val, digits=10))
	end
end

# ╔═╡ 56d0de38-5639-4196-aafe-79a9ab933980
begin
	samples_cube = Dict()
	
	for (transit, cube) in cubes
		samples_cube[transit] = Dict(
			param => scale_samples(
				pyconvert(Dict{String, Vector}, cube["samples"])[param],
				PARAMS[param],
				BMA
			)
			for param in keys(PARAMS)
		)
	end
end

# ╔═╡ c7791371-9010-44f6-9aca-1d50f8e43ad4
Makie.ColorSchemes.Paired_8[3] |> Makie.Colors.hex

# ╔═╡ 1f7b883c-0192-45bd-a206-2a9fde1409ca
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
			"#B2DF8A",  # Green
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
			Label = (textsize=18,  padding=(0, 10, 0, 0)),
			Lines = (linewidth=3,),
			Scatter = (linewidth=10,),
			palette = (color=COLORS, patchcolor=[(c, 0.35) for c in COLORS]),
			fontsize = 18,
			rowgap = 5,
			colgap = 5,
		)
	)

	COLORS
end

# ╔═╡ 89c48710-651e-45ff-8fcb-e4173559defd
function plot_lc!(gl, i, transit, cube; ax_top_kwargs=(), ax_bottom_kwargs=())
	ax_top = Axis(gl[1, 1]; ax_top_kwargs...)
	ax_bottom = Axis(gl[2, 1]; ax_bottom_kwargs...)

	t₀ = val(cube["results"], "t0")
	P = val(cube["results"], "P")
	t = ϕ.(pyconvert(Vector, cube["models"]["t"]), t₀, P)
	t_interp = ϕ.(pyconvert(Vector, cube["models"]["t_interp"]), t₀, P)
	LC_det = pyconvert(Vector, cube["models"]["LC_det"])
	LC_transit_model = pyconvert(Vector, cube["models"]["LC_transit_model"])
	LC_det_model_interp = pyconvert(Vector, cube["models"]["LC_det_model_interp"])
	resids = LC_det - LC_transit_model
	resids_σ = round(Int, std(resids) * 1e6)
	resids .*= 1e6

	color = COLORS[i]
	color_dark = 0.75 * COLORS[i]

	# Top panel
	scatter!(ax_top,
		t,
		LC_det,
		color = color,
	)
	lines!(ax_top,
		t_interp,
		LC_det_model_interp,
		color = color_dark,
	)
	t_x = 0.06
	text!(ax_top, "$transit";
		position = (t_x, 0.9765),
		align = (:right, :bottom),
		color = color_dark,
	)

	# Bottom panel
	scatter!(ax_bottom, t, resids, color=color)
	lines!(ax_bottom, t_interp, zero(t_interp), color=color_dark)
	text!(ax_bottom, "$resids_σ ppm";
		position = (t_x, 2300),
		align = (:right, :top),
		color = color_dark,
	)

	linkxaxes!(ax_top, ax_bottom)
	hidexdecorations!(ax_top)

	rowsize!(gl, 2, Relative(1/3))
end

# ╔═╡ 4be0d7b7-2ea5-4c4d-92b9-1f8109014e12
begin
	fig = Figure(resolution=FIG_LARGE)
	
	grid = CartesianIndices((2, 2))
	xlims = -0.065, 0.065
	ax_top_kwargs = (
		ylabel = "Relative flux",
		#limits = (xlims..., 0.975, 1.005),
	)
	ax_bottom_kwargs = (
		xlabel = "Phase",
		ylabel = "Residuals (ppm)",
		#limits=(xlims..., -3000, 3000),
	)
	for (i, (transit, cube)) ∈ enumerate(cubes)
		g_idxs = grid[i]
		gl = fig[g_idxs.I...] = GridLayout()
		plot_lc!(gl, i, transit, cube; ax_top_kwargs, ax_bottom_kwargs)
	end
	
	axs = reshape(fig.content, 4, 2)
	linkxaxes!(axs...)
	linkyaxes!.(axs[1, 1], axs[3, 1], axs[1, 2], axs[3, 2])
	linkyaxes!.(axs[2, 1], axs[4, 1], axs[2, 2], axs[4, 2])
	hidexdecorations!.(axs[2, :])
	hideydecorations!.(axs[:, 2])

	if occursin("sp", DATA_DIR)
		savefig(fig, "$(FIG_DIR)/detrended_wlcs_sp.png")
	else
		savefig(fig, "$(FIG_DIR)/detrended_wlcs.png")
	end
	
	fig
end

# ╔═╡ d5ff9b30-00dd-41d3-9adf-ff7905d71ae8
let
	n_params = length(PARAMS) # Number of fitted parameters
	
	# Create empty corner plot grid
	fig = Figure(resolution=(1_400, 1_400))
	
	for j in 1:n_params, i in 1:n_params
		# Create subplot apply global settings
		ax = Axis(fig[i, j];
			aspect = 1,
			xticklabelrotation = π/4,
			xticks = LinearTicks(3),
			yticks = LinearTicks(3),
			limits = ((-15, 15), (-15, 15)),
		)
		# Hide upper triangle
		j > i && (hidedecorations!(ax); hidespines!(ax))
		# Hide y ticks on diagonals
		j == i && (hideydecorations!(ax); ylims!(0, 1.2))
		# Hide x ticks on all diagonal elements except the bottom one
		j == i && i != n_params && (
			hidexdecorations!(ax, grid=false);
		)
		# Hide ticks on interior lower triangle
		j < i && i != n_params && j != 1 && (
			hideydecorations!(ax, grid=false);
			hidexdecorations!(ax, grid=false);
		)
		# Hide remaining xyticks
		j < i && j == 1 && i != n_params && (
			hidexdecorations!(ax, grid=false);
		)
		j < i && i == n_params && j != 1 && (
			hideydecorations!(ax, grid=false);
		)
	end
			
	# Plot corners from each night
	elems = MarkerElement[]
	elem_labels = String[]
	for (i, (transit, cube)) in enumerate(cubes)
		c = COLORS[i]
		plot_corner!(fig,
			samples_cube[transit], keys(PARAMS), color=c
		)
		push!(elems,
			MarkerElement(marker='◇', color=c, strokecolor=c)
		)
		push!(elem_labels, transit)
	end
	
	# Align axes limits and apply labels
	axs = reshape(copy(fig.content), n_params, n_params)
	[linkxaxes!(reverse(axs[:, j])...) for j in 1:n_params]
	for (j, (param, param_latex)) in enumerate(PARAMS)
		axs[end, j].xlabel = "Δ"*param_latex
		axs[j, begin].ylabel = "Δ"*param_latex
	end
	
	Legend(fig[1, end], elems, elem_labels;
		halign = :right,
		valign = :top,
		patchsize = (25, 25),
		tellwidth = false,
		tellheight = false,
		rowgap = 10,
		labelsize = 25,
		nbanks = 2,
		orientation = :horizontal,
		markersize = 35,
		markerstrokewidth = 1,
	)

	if occursin("sp", DATA_DIR)
		savefig(fig, "$(FIG_DIR)/detrended_wlcs_corner_sp.png")
	else
		savefig(fig, "$(FIG_DIR)/detrended_wlcs_corner.png")
	end
	
	fig
end

# ╔═╡ Cell order:
# ╟─506eeeb2-e56d-436b-91b8-605e52201563
# ╠═25d1284c-7260-4f3a-916a-b2814d2606af
# ╟─782806a6-efd2-45a9-b898-788a276c282b
# ╟─a8d91138-73e7-4382-a032-37daec54a9c0
# ╠═2191791b-df62-4f1b-88bf-060cc47896b2
# ╠═f539e06d-a1b5-413a-b90e-91cb0bbd5a4c
# ╠═e873f9c6-fd1a-4227-9df1-70c626e4a0a1
# ╠═e72dba55-6a33-462f-aeac-5f62b25cb46a
# ╟─579e62da-7ffb-4639-bd73-3826ade1cfa2
# ╟─a8cf11e2-796e-45ff-bdc9-e273b927700e
# ╟─ae82d3c1-3912-4a5e-85f5-6383af42291e
# ╠═4be0d7b7-2ea5-4c4d-92b9-1f8109014e12
# ╠═88a92472-f20d-4908-8bfd-fcbf11f80b8a
# ╠═89c48710-651e-45ff-8fcb-e4173559defd
# ╠═0dd63eaf-1afd-4caf-a74b-7cd217b3c515
# ╠═d43ec3eb-1d5e-4a63-b5e8-8dcbeb57ae7c
# ╟─b28bb1b6-c148-41c4-9f94-0833e365cad4
# ╟─30b84501-fdcd-4d83-b929-ff354de69a17
# ╠═c7a179a3-9966-452d-b430-a28b2f004bc5
# ╠═694a6067-390e-4b79-b8b1-87e0c6d15f48
# ╠═19fcaa15-6f01-46a6-8225-4b5cafd89cc1
# ╠═de0a4468-56aa-4748-80a0-6c9ab6b8579e
# ╠═bbc14e57-57fe-4811-91d4-d07b102cfa5d
# ╠═25dd0c88-089b-406b-ac0f-6f21c57fe986
# ╟─22c72eeb-8e32-4d7c-86c8-ab117735769e
# ╠═a96ee19b-39fe-4f5b-bc75-124b4e422713
# ╠═7cbadc33-5b32-4816-8092-09054c64073f
# ╠═47278372-b311-4ea7-bfa4-82b8f95c97fa
# ╠═f3a8f6fb-023c-4077-8c73-7502e56eb607
# ╠═ee9347b2-e97d-4f66-9c21-7487ca2c2e30
# ╟─c936b76a-636b-4f10-b556-fa19808c1562
# ╠═d279e93e-8665-41b2-bd5c-723458fabe86
# ╟─68ec4343-5f6c-4dfd-90b5-6393b4c819b9
# ╟─e452d2b1-1010-4ce3-8d32-9e9f1d0dfa0b
# ╠═56d0de38-5639-4196-aafe-79a9ab933980
# ╠═706f1fb6-2895-48f6-a315-842fbf35da18
# ╠═18c90ed4-2f07-493f-95b2-e308cd7a03a9
# ╠═7fdd0cb7-7af3-4ca1-8939-9d9b7d6e9527
# ╠═8969ff2e-d160-4eb3-8f77-f8c2006518ff
# ╠═266b5871-22e5-4b4a-80b1-468a9ddd9193
# ╠═3225ad37-6264-48b8-b1fc-c6f0136f8beb
# ╠═8f795cda-b905-4c89-ade5-05e35b10d4b3
# ╠═d5ff9b30-00dd-41d3-9adf-ff7905d71ae8
# ╠═8c9428a9-f2ad-4329-97a3-97ffa6b40f28
# ╠═2cbc6ddb-210e-41e8-b745-5c41eba4e778
# ╠═ed935d16-ddce-4334-a880-005732b38936
# ╠═6fcd1377-8364-45a3-9ff6-89d61df1ef42
# ╟─807e913f-d8c3-41e9-acb3-4c024dedd67b
# ╠═2bfaa81e-f0bd-4527-a5fe-091589b1f76e
# ╟─30ae3744-0e7e-4c16-b91a-91eb518fba5b
# ╟─bf9c0b95-fe17-425d-8904-8298f7e5451c
# ╠═db539901-f0b0-4692-a8d2-6c72dff41196
# ╠═c7791371-9010-44f6-9aca-1d50f8e43ad4
# ╠═1f7b883c-0192-45bd-a206-2a9fde1409ca
# ╠═5939e2a3-8407-4579-8274-e3891acbabb1
