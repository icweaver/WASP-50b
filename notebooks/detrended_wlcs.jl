### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 691eddff-f2eb-41a8-ab05-63afb46d15f2
begin
	import PlutoUI as pl
	using CairoMakie
	using Colors
	using CSV
	using DataFrames
	using Glob
	using Measurements
	using PyCall
	using Statistics, KernelDensity
	using Latexify, OrderedCollections
	using DataFramesMeta
end

# ╔═╡ 506eeeb2-e56d-436b-91b8-605e52201563
md"""
# Detrended white light curves

In this notebook we will visualize the detrended white light curves from IMACS and LDSS3. We used the average orbital and system parameters obtained from these detrended fits to place uniform constraints on the binned wavelength analysis **<ADD LINK>.**

!!! note "TODO"
	Add LDSS3 plots

$(pl.TableOfContents())
"""

# ╔═╡ 782806a6-efd2-45a9-b898-788a276c282b
md"""
## Load data

First, let's load the relevant data needed for this notebook:
"""

# ╔═╡ 3f0f5777-00f1-443d-8ced-d901550010d3
const DATA_DIR = "data/detrended/out_l_C/WASP50"

# ╔═╡ 579e62da-7ffb-4639-bd73-3826ade1cfa2
md"""
The data cube is organized by night as follows:

```julia
cubes
├── Transit_1
├── Transit_2
└── Transit_3
    ├── samples
    │   └── Dict
    ├── models
    │   └── Dict
    └── results
        ├── Dict
        └── Table
```

where `samples` is a dictionary of the Bayesian model averaged posterior samples for each fitted parameter, `models` is a dictionary of the final detrended light curve and associated Gaussian Process parameters, and `results` is a table summarizing the mean and +/- 1σ uncertainty associated with each sample.

!!! note "TODO"
	Update with LDSS3 paths 
"""

# ╔═╡ 583b3377-f4f6-4170-8658-d3ba17a5b86d
md"""
### Helper functions
"""

# ╔═╡ 39dbca86-a4b9-11eb-1c64-9ddf1a9990ab
begin
	py"""
	import numpy as np
	import pickle
	
	def load_npz(fpath, allow_pickle=False):
		return np.load(fpath, allow_pickle=allow_pickle)[()]
	
	def load_pickle(fpath):
		with open(fpath, "rb") as f:
			data = pickle.load(f)
		return data
	"""
	load_npz(s; allow_pickle=false) = py"load_npz"(s, allow_pickle=allow_pickle)
	load_pickle(s) = py"load_pickle"(s)
end;

# ╔═╡ 2191791b-df62-4f1b-88bf-060cc47896b2
begin
	# Load IMACS
	cubes = OrderedDict(
		"Transit $i IMACS" => Dict(
			
			"samples" => load_pickle(fpath_sample),
			
			"models" => load_npz(fpath_model, allow_pickle=true),
			
			"results" => CSV.File(
				fpath_result,
				comment = "#",
				normalizenames=true,
			) |> DataFrame
		)
		
		for (i, (fpath_sample, fpath_model, fpath_result)) in enumerate(zip(
			sort(glob("$(DATA_DIR)/w50*IMACS*/white-light/BMA_posteriors.pkl")),
			sort(glob("$(DATA_DIR)/w50*IMACS*/white-light/BMA_WLC.npy")),
			sort(glob("$(DATA_DIR)/w50*IMACS*/white-light/results.dat")),
		))
	)
	
	
	# Load LDSS3
# 	for (fpath_sample, fpath_model, fpath_result) in zip(
# 			sort(glob("$(DATA_DIR)/w50*LDSS3*/white-light/BMA_posteriors.pkl")),
# 			sort(glob("$(DATA_DIR)/w50*LDSS3*/white-light/BMA_WLC.npy")),
# 			sort(glob("$(DATA_DIR)/w50*LDSS3*/white-light/results.dat")),
# 		)
				
# 		cubes["Transit 2 LDSS3 $(isflat(fpath_sample))"] = OrderedDict(
			
# 			"samples" => load_pickle(fpath_sample),
			
# 			"models" => load_npz(fpath_model, allow_pickle=true),
			
# 			"results" => CSV.File(
# 				fpath_result,
# 				comment = "#",
# 				normalizenames=true,
# 			) |> DataFrame
# 		)
# 	end
end

# ╔═╡ d79dbe8c-effc-4537-b0a1-6a3bcb5db2e5
cubes |> keys

# ╔═╡ d0172076-d245-4246-b35f-49cc6f1cde0b
# Return `flat` or `noflat` for data organization
isflat(fpath) = split(split(fpath, "LDSS3_")[2], '/')[1]

# ╔═╡ a8cf11e2-796e-45ff-bdc9-e273b927700e
md"""
## Transit curves ⚪
"""

# ╔═╡ ae82d3c1-3912-4a5e-85f5-6383af42291e
md"""
Plotting the data from the `models` cube returns the following detrended white light curves:
"""

# ╔═╡ 0dd63eaf-1afd-4caf-a74b-7cd217b3c515
# Returns value `v` from `Variables` columns in results.dat file 
val(df, v) = @where(df, :Variable .== v)[1, "Value"]

# ╔═╡ d43ec3eb-1d5e-4a63-b5e8-8dcbeb57ae7c
# Computes orbital phase
function ϕ(t, t₀, P)
	phase = ((t - t₀) / P) % 1.0
    phase ≥ 0.5 && (phase -= 1.0)
	return phase
end

# ╔═╡ b28bb1b6-c148-41c4-9f94-0833e365cad4
md"""
## Summary table 📓
"""

# ╔═╡ 30b84501-fdcd-4d83-b929-ff354de69a17
md"""
We summarize the Bayesian Model Averag (BMA) results for selected parameters for each night below:
"""

# ╔═╡ 1bc1110f-e57d-4f31-a309-9b4e1aed1c0a
md"""
Finally, we average together each parameter from each night, weighted by its maximum uncertainty per night:
"""

# ╔═╡ c936b76a-636b-4f10-b556-fa19808c1562
md"""
### Save to file
"""

# ╔═╡ 68ec4343-5f6c-4dfd-90b5-6393b4c819b9
md"""
## Corner plots 📐
"""

# ╔═╡ e452d2b1-1010-4ce3-8d32-9e9f1d0dfa0b
md"""
To visualize the spread of each parameter, we define the dimensionless metric:

```math
Δx \equiv \frac{x - \overline{\text{BMA}}_μ}{\overline{\text{BMA}}_σ}\quad,
```

where $x$ is a sample from the posterior distribution for the given parameter, ``\overline{\text{BMA}}_μ`` is the BMA averaged across nights, and ``\overline{\text{BMA}}_σ`` is the maximum uncertainty, also averaged across nights. These values correspond to the `Combined` column in the table above. ``\Delta x`` is then a measure of the displacement of each sample in its posterior distribution from its corresponding average BMA value, scaled by the uncertainty in the average BMA.
"""

# ╔═╡ 706f1fb6-2895-48f6-a315-842fbf35da18
function scale_samples(samples, param, BMA)
	m = @where(BMA, :Parameter .== param)[!, "Combined"][1]
	return (samples .- m.val) ./ m.err
end

# ╔═╡ ed935d16-ddce-4334-a880-005732b38936
# Params to show in corner plot
const PARAMS = OrderedDict(
	"p" => "Rₚ/Rₛ",
	"t0" => "t₀",
	"P" => "P",
	"rho" => "ρₛ",
	# "aR" => "a/Rₛ",
	# "inc" => "i",
	# "b" => "b",
	# "q1" => "u",
);

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
	
	summary_tables
end

# ╔═╡ de0a4468-56aa-4748-80a0-6c9ab6b8579e
BMA_matrix = hcat((
	summary[!, "Value"] .± maximum((summary[!, "SigmaUp"], summary[!, "SigmaDown"]))
	for summary in summary_tables
)...) |> x-> hcat(x, mean(x, dims=2));

# ╔═╡ 19fcaa15-6f01-46a6-8225-4b5cafd89cc1
BMA = DataFrame(
	[PARAMS.vals BMA_matrix],
	["Parameter", keys(cubes)..., "Combined"]
);

# ╔═╡ c7a179a3-9966-452d-b430-a28b2f004bc5
latexify(BMA)

# ╔═╡ d279e93e-8665-41b2-bd5c-723458fabe86
# Will probably just copy-paste directly into paper
pl.with_terminal() do
	BMA |> x -> latexify(x, env=:table) |> print
end;

# ╔═╡ 56d0de38-5639-4196-aafe-79a9ab933980
begin
	samples_cube = Dict()
	
	for (transit, cube) in cubes
		samples_cube[transit] = Dict(
			param => scale_samples(cube["samples"][param], PARAMS[param], BMA)
			for param in keys(PARAMS)
		)
	end
end

# ╔═╡ c5e10e47-ec64-4911-8107-487d1ef3f134
md"""
### Helper functions
"""

# ╔═╡ 6fcd1377-8364-45a3-9ff6-89d61df1ef42
# Number of levels `n` to show in contour plots
levels(A, n) = reverse(
	range(maximum(A), step=-maximum(A)/(n+1), length=(n+1))
)

# ╔═╡ 2cbc6ddb-210e-41e8-b745-5c41eba4e778
function plot_corner!(fig, samples, params; n_levels=4, color=:blue)
	for (j, p1) in enumerate(params), (i, p2) in enumerate(params)
		# 1D plot
		if i == j
			density!(fig[i, i], samples[p1];
				color = (color, 0.125),
				strokewidth = 3,
				strokecolor = color,
				strokearound = true,
			)
		end
		# 2D plot
		if i > j
			Z = kde((samples[p1], samples[p2]), npoints=(2^4, 2^4),)
			contourf!(fig[i, j], Z;
				levels = levels(Z.density, n_levels),
				colormap = cgrad(
					range(colorant"white", color, length=n_levels),
					alpha=0.5
				),
			)
			contour!(fig[i, j], Z;
				levels = levels(Z.density, n_levels),
				color = color,
				linewidth = 3,
			)
		end
	end
end

# ╔═╡ 82a23101-9e1f-4eae-b529-e750a44c98b1
md"""
!!! note
	The sampled values have been scaled so that the distance of each sample from its literature "truth" value is in units of that truth value's reported uncertainty. We show this operation below.
"""

# ╔═╡ 30ae3744-0e7e-4c16-b91a-91eb518fba5b
md"""
## Plot configs
"""

# ╔═╡ 940ebaf2-659a-4319-bbe6-e0290752f1fb
#const COLORS = cgrad(:seaborn_colorblind6, categorical=true)

# ╔═╡ feee4fd1-e16d-4d9b-8bc0-2c4f7afb0c43
const COLORS = parse.(Colorant,
	[
		"#5daed9",  # Cyan
		"plum",
		"#f7ad4d",  # Yellow
		"mediumaquamarine",
		"slategray",
		"#126399",  # Blue
		"#956cb4",  # Purple
		"#ff7f00",  # Orange
		"#029e73",  # Green
	]
)

# ╔═╡ 4be0d7b7-2ea5-4c4d-92b9-1f8109014e12
let
	fig = Figure(resolution=(700, 800))
	
	i = 1
	for (i, (transit, cube)) in enumerate(cubes)
		t₀ = val(cube["results"], "t0")
		P = val(cube["results"], "P")
		
		scatter(fig[i, 1],
			ϕ.(cube["models"]["t"], t₀, P),
			cube["models"]["LC_det"],
			color = COLORS[i],
			strokewidth = 0,
			#axis = (title=transit,),
		)
		
		lines!(fig[i, 1],
			ϕ.(cube["models"]["t_interp"], t₀, P),
			cube["models"]["LC_det_model_interp"],
			color = 0.75*COLORS[i],
			linewidth = 3,
		)

		text!(fig[i, 1], transit;
			position = Point2f0(0.05, 0.980),
			space = :data,
			textsize = 0.0025,
			align = (:center, :bottom),
			color = 0.75*COLORS[i],
		)
	end
	
	axs = reshape(copy(fig.content), (length(cubes), 1))
	linkaxes!(axs...)
	#ylims!(axs[end], 0.96, 1.02)
	hidexdecorations!.(axs[begin:end-1], grid=false)
	axs[end].xlabel = "Phase"
	axs[end].ylabel = "Relative flux"
	
	fig |> pl.as_png
end

# ╔═╡ d5ff9b30-00dd-41d3-9adf-ff7905d71ae8
let
	n_params = length(PARAMS) # Number of fitted parameters
	
	# Create empty corner plot grid
	fig = Figure(resolution=(1_400, 1_400))
	
	for j in 1:n_params, i in 1:n_params
		# Create subplot apply global settings
		ax = Axis(fig[i, j], axis=(aspect=1, xticklabelrotation=π/4),)
		#ax.xticklabelrotation = π/4
		ax.aspect = 1.0
		# Hide upper triangle
		j > i && (hidedecorations!(ax); hidespines!(ax))
		# Hide y ticks on diagonals
		j == i && hideydecorations!(ax)
		# Hide x ticks on all diagonals except the bottom one
		j == i && i != n_params && hidexdecorations!(ax, grid=false)
		# Hide ticks on interior lower triangle
		j < i && i != n_params && j != 1 && (
			hideydecorations!(ax, grid=false);
			hidexdecorations!(ax, grid=false))
		# Hide remaining xyticks
		j < i && j == 1 && i != n_params && hidexdecorations!(ax, grid=false)
		j < i && i == n_params && j != 1 && hideydecorations!(ax, grid=false) 
	end
			
	# Plot corners from each night
	diamond = Point2f0[(0.5, 0), (1, 0.5), (0.5, 1), (0, 0.5), (0.5, 0.0)]
	elems = LineElement[]
	elem_labels = String[]
	for (i, (transit, cube)) in enumerate(cubes)
		c = COLORS[i]
		plot_corner!(fig,
			samples_cube[transit], keys(PARAMS), color=c
		)
		push!(elems,
			LineElement(color=c, linestyle=nothing, linepoints=diamond,)
		)
		push!(elem_labels, transit)
	end
	
	# Align axes limits and apply labels
	axs = reshape(copy(fig.content), n_params, n_params)
	[linkxaxes!(reverse(axs[:, j])...) for j in 1:n_params]
	for (j, (param, param_latex)) in enumerate(PARAMS)
		axs[end, j].xlabel = "Δ"*param_latex
		axs[end, j].xlabelsize = 26
		axs[j, begin].ylabel = "Δ"*param_latex
		axs[j, begin].ylabelsize = 26
	end
	
	#Legend(fig[1+1, n_params-1],
	Legend(fig[1, 3],
		elems,
		elem_labels,
		patchsize = (25, 25),
		tellwidth = false,
		tellheight = false,
		rowgap = 10,
		labelsize = 25,
	)
	
	fig |> pl.as_png
end

# ╔═╡ a9747b5e-adf9-48dd-96c2-f184d873d1ac
set_theme!(palette = (color=COLORS,),)

# ╔═╡ baeadfce-535a-46c3-8cb9-79cf6bde8555
md"""
## Packages
"""

# ╔═╡ Cell order:
# ╟─506eeeb2-e56d-436b-91b8-605e52201563
# ╟─782806a6-efd2-45a9-b898-788a276c282b
# ╠═3f0f5777-00f1-443d-8ced-d901550010d3
# ╠═d79dbe8c-effc-4537-b0a1-6a3bcb5db2e5
# ╠═2191791b-df62-4f1b-88bf-060cc47896b2
# ╟─579e62da-7ffb-4639-bd73-3826ade1cfa2
# ╟─583b3377-f4f6-4170-8658-d3ba17a5b86d
# ╠═39dbca86-a4b9-11eb-1c64-9ddf1a9990ab
# ╠═d0172076-d245-4246-b35f-49cc6f1cde0b
# ╟─a8cf11e2-796e-45ff-bdc9-e273b927700e
# ╟─ae82d3c1-3912-4a5e-85f5-6383af42291e
# ╠═4be0d7b7-2ea5-4c4d-92b9-1f8109014e12
# ╠═0dd63eaf-1afd-4caf-a74b-7cd217b3c515
# ╠═d43ec3eb-1d5e-4a63-b5e8-8dcbeb57ae7c
# ╟─b28bb1b6-c148-41c4-9f94-0833e365cad4
# ╟─30b84501-fdcd-4d83-b929-ff354de69a17
# ╠═ee9347b2-e97d-4f66-9c21-7487ca2c2e30
# ╟─1bc1110f-e57d-4f31-a309-9b4e1aed1c0a
# ╠═c7a179a3-9966-452d-b430-a28b2f004bc5
# ╠═19fcaa15-6f01-46a6-8225-4b5cafd89cc1
# ╠═de0a4468-56aa-4748-80a0-6c9ab6b8579e
# ╟─c936b76a-636b-4f10-b556-fa19808c1562
# ╠═d279e93e-8665-41b2-bd5c-723458fabe86
# ╟─68ec4343-5f6c-4dfd-90b5-6393b4c819b9
# ╟─e452d2b1-1010-4ce3-8d32-9e9f1d0dfa0b
# ╠═56d0de38-5639-4196-aafe-79a9ab933980
# ╠═706f1fb6-2895-48f6-a315-842fbf35da18
# ╠═d5ff9b30-00dd-41d3-9adf-ff7905d71ae8
# ╠═ed935d16-ddce-4334-a880-005732b38936
# ╟─c5e10e47-ec64-4911-8107-487d1ef3f134
# ╠═2cbc6ddb-210e-41e8-b745-5c41eba4e778
# ╠═6fcd1377-8364-45a3-9ff6-89d61df1ef42
# ╟─82a23101-9e1f-4eae-b529-e750a44c98b1
# ╟─30ae3744-0e7e-4c16-b91a-91eb518fba5b
# ╠═940ebaf2-659a-4319-bbe6-e0290752f1fb
# ╠═feee4fd1-e16d-4d9b-8bc0-2c4f7afb0c43
# ╠═a9747b5e-adf9-48dd-96c2-f184d873d1ac
# ╟─baeadfce-535a-46c3-8cb9-79cf6bde8555
# ╠═691eddff-f2eb-41a8-ab05-63afb46d15f2
