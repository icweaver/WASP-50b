### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 6a9b28b5-bc77-4c9f-86e4-a564d57feb72
using Latexify, OrderedCollections

# ╔═╡ 4cba0dcb-5d6e-42d2-baff-134dbe8500ba
using DataFramesMeta

# ╔═╡ 691eddff-f2eb-41a8-ab05-63afb46d15f2
begin
	import PlutoUI as pl
	import JSON
	using CairoMakie
	using Colors
	using CSV
	using DataFrames
	using Glob
	using Measurements
	using PyCall
	using Statistics, KernelDensity
end

# ╔═╡ 506eeeb2-e56d-436b-91b8-605e52201563
md"""
# Detrended white light curves

In this notebook we will visualize the detrended white light curves. We used the average orbital and system parameters obtained from these detrended fits to place uniform constraints on the binned wavelength analysis.

$(pl.TableOfContents())
"""

# ╔═╡ 782806a6-efd2-45a9-b898-788a276c282b
md"""
## Load data

First, let's load the relevant data needed for this notebook:
"""

# ╔═╡ 3f0f5777-00f1-443d-8ced-d901550010d3
const DATA_DIR = "data/detrended_wlcs/out_l/WASP50"

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
cubes = Dict(
	"Transit $i" => Dict(
		
		"samples" => load_pickle(fpath_sample),
		
		"models" => load_npz(fpath_model, allow_pickle=true),
		
		"results" => CSV.File(
			fpath_result,
			comment = "#",
			delim = ' ',			
			ignorerepeated=true,
		) |> DataFrame
	)
	
	for (i, (fpath_sample, fpath_model, fpath_result)) in enumerate(zip(
		sort(glob("$(DATA_DIR)/w50_*/white-light/BMA_posteriors.pkl")),
		sort(glob("$(DATA_DIR)/w50_*/white-light/BMA_WLC.npy")),
		sort(glob("$(DATA_DIR)/w50_*/white-light/results.dat")),
	))
)

# ╔═╡ a8cf11e2-796e-45ff-bdc9-e273b927700e
md"""
## Transit curves ⚪
"""

# ╔═╡ ae82d3c1-3912-4a5e-85f5-6383af42291e
md"""
Plotting the data from the `models` cube returns the following detrended white light curves:
"""

# ╔═╡ b28bb1b6-c148-41c4-9f94-0833e365cad4
md"""
## Summary table 📓
"""

# ╔═╡ 30b84501-fdcd-4d83-b929-ff354de69a17
md"""
We summarize the Bayesian modeled results for selected parameters for each night:
"""

# ╔═╡ 1bc1110f-e57d-4f31-a309-9b4e1aed1c0a
md"""
Finally, we average together each parameter from each night, weighted by its maximum uncertainty per night:
"""

# ╔═╡ 68ec4343-5f6c-4dfd-90b5-6393b4c819b9
md"""
## Corner plots 📐
"""

# ╔═╡ e452d2b1-1010-4ce3-8d32-9e9f1d0dfa0b
md"""
The `samples` cube returns the following corner plot for the fitted `PARAMS` below for each night:
"""

# ╔═╡ ed935d16-ddce-4334-a880-005732b38936
const PARAMS = OrderedDict(
	"p" => "Rₚ/Rₛ",
	"t0" => "t₀",
	"P" => "P",
	# "rho" => "ρₛ",
	# "aR" => "a/Rₛ",
	# "inc" => "i",
	# "b" => "b",
	# "q1" => "u"
)

# ╔═╡ ee9347b2-e97d-4f66-9c21-7487ca2c2e30
summary_tables = filter.(
	row -> row.Variable in keys(PARAMS),
	(cube["results"] for (transit, cube) in cubes)
)

# ╔═╡ de0a4468-56aa-4748-80a0-6c9ab6b8579e
m_results = hcat((
	summary[!, "Value"] .± maximum((summary[!, "SigmaUp"], summary[!, "SigmaDown"]))
	for summary in summary_tables
)...) |> x-> hcat(x, mean(x, dims=2));

# ╔═╡ 19fcaa15-6f01-46a6-8225-4b5cafd89cc1
BMA = DataFrame(
	[PARAMS.vals m_results],
	[:Parameter, :Transit_1, :Transit_2, :Transit_3, :Combined]
)

# ╔═╡ d279e93e-8665-41b2-bd5c-723458fabe86
pl.with_terminal() do
BMA |> latexify |> print
end

# ╔═╡ c5e10e47-ec64-4911-8107-487d1ef3f134
md"""
### Helper functions
"""

# ╔═╡ 6fcd1377-8364-45a3-9ff6-89d61df1ef42
levels(A, n) = reverse(
	range(maximum(A), step=-maximum(A)/(n+1), length=(n+1))
)

# ╔═╡ 2cbc6ddb-210e-41e8-b745-5c41eba4e778
function plot_corner!(fig, cube, params; n_levels=4, color=:blue)
	for (j, p1) in enumerate(params), (i, p2) in enumerate(params)
		if i == j
			density!(fig[i, i], cube[p1];
				color = (color, 0.125),
				strokewidth = 3,
				strokecolor = color,
				strokearound = true,
			)
		end
		if i > j
			Z = kde((cube[p1], cube[p2]), npoints=(2^4, 2^4),)
			contourf!(fig[i, j], Z;
				levels = levels(Z.density, n_levels),
				colormap = cgrad(range(colorant"white", color), alpha=0.5),
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

# ╔═╡ ddd9a95b-735a-4995-8893-542128bb56d6
truths = JSON.parsefile("$(DATA_DIR)/truth.json")

# ╔═╡ a0b7fc53-6d0a-47c1-94b9-6b742c014a93
Measurements.uncertainty.(BMA[!, "Combined"])

# ╔═╡ bfdf906b-1e54-4a77-820e-0ea43fe67210
@where(BMA, :Parameter .== "P").Combined[1]

# ╔═╡ 82647208-82e1-4e5a-b6f0-8c77db8b230a
String.(keys(PARAMS))

# ╔═╡ b676f75f-51e1-4ea8-acda-67d62607465e


# ╔═╡ 931ce3d5-c4ed-496c-883b-d7ee33e957cc
function adjust_dict(dict, params, truths)
	d = filter!(p -> p.first ∈ keys(params), dict)
	for param in params_str
		truth_val = truths[param]["truth"][1]
		truth_val_err = maximum((truths[param]["truth"][2:3])) 
		d[param] = @. ((d[param] - truth_val) / truth_val_err)
	end
	return d
end

# ╔═╡ 30ae3744-0e7e-4c16-b91a-91eb518fba5b
md"""
## Plot configs
"""

# ╔═╡ 940ebaf2-659a-4319-bbe6-e0290752f1fb
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

# ╔═╡ 4be0d7b7-2ea5-4c4d-92b9-1f8109014e12
let
	fig = Figure(resolution=(700, 800))
	
	i = 1
	for (i, (transit, cube)) in enumerate(cubes)
			scatter(
			fig[i, 1], cube["models"]["t"] .- 2.45e6,
			cube["models"]["LC_det"],
			color = COLORS[i],
			strokewidth = 0,
			axis = (title=transit,),
		)
		lines!(
			fig[i, 1],
			cube["models"]["t_interp"] .- 2.45e6,
			cube["models"]["LC_det_model_interp"],
			color = 0.75*COLORS[i],
			linewidth = 3,
		)
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
		ax = Axis(fig[i, j], axis=(aspect=1, xticklabelrotation=π/4),)
		ax.xticklabelrotation = π/4
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
	for (i, (transits, cube)) in enumerate(cubes)
		plot_corner!(fig, cube["samples"], keys(PARAMS), color=COLORS[i])
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
	
	fig
end

# ╔═╡ baeadfce-535a-46c3-8cb9-79cf6bde8555
md"""
## Packages
"""

# ╔═╡ Cell order:
# ╟─506eeeb2-e56d-436b-91b8-605e52201563
# ╟─782806a6-efd2-45a9-b898-788a276c282b
# ╠═3f0f5777-00f1-443d-8ced-d901550010d3
# ╠═2191791b-df62-4f1b-88bf-060cc47896b2
# ╟─579e62da-7ffb-4639-bd73-3826ade1cfa2
# ╟─583b3377-f4f6-4170-8658-d3ba17a5b86d
# ╠═39dbca86-a4b9-11eb-1c64-9ddf1a9990ab
# ╟─a8cf11e2-796e-45ff-bdc9-e273b927700e
# ╟─ae82d3c1-3912-4a5e-85f5-6383af42291e
# ╠═4be0d7b7-2ea5-4c4d-92b9-1f8109014e12
# ╟─b28bb1b6-c148-41c4-9f94-0833e365cad4
# ╟─30b84501-fdcd-4d83-b929-ff354de69a17
# ╠═ee9347b2-e97d-4f66-9c21-7487ca2c2e30
# ╟─1bc1110f-e57d-4f31-a309-9b4e1aed1c0a
# ╠═19fcaa15-6f01-46a6-8225-4b5cafd89cc1
# ╠═de0a4468-56aa-4748-80a0-6c9ab6b8579e
# ╟─68ec4343-5f6c-4dfd-90b5-6393b4c819b9
# ╟─e452d2b1-1010-4ce3-8d32-9e9f1d0dfa0b
# ╠═d5ff9b30-00dd-41d3-9adf-ff7905d71ae8
# ╠═ed935d16-ddce-4334-a880-005732b38936
# ╠═6a9b28b5-bc77-4c9f-86e4-a564d57feb72
# ╠═d279e93e-8665-41b2-bd5c-723458fabe86
# ╟─c5e10e47-ec64-4911-8107-487d1ef3f134
# ╠═2cbc6ddb-210e-41e8-b745-5c41eba4e778
# ╠═6fcd1377-8364-45a3-9ff6-89d61df1ef42
# ╟─82a23101-9e1f-4eae-b529-e750a44c98b1
# ╠═ddd9a95b-735a-4995-8893-542128bb56d6
# ╠═a0b7fc53-6d0a-47c1-94b9-6b742c014a93
# ╠═4cba0dcb-5d6e-42d2-baff-134dbe8500ba
# ╠═bfdf906b-1e54-4a77-820e-0ea43fe67210
# ╠═82647208-82e1-4e5a-b6f0-8c77db8b230a
# ╠═b676f75f-51e1-4ea8-acda-67d62607465e
# ╠═931ce3d5-c4ed-496c-883b-d7ee33e957cc
# ╟─30ae3744-0e7e-4c16-b91a-91eb518fba5b
# ╠═940ebaf2-659a-4319-bbe6-e0290752f1fb
# ╟─baeadfce-535a-46c3-8cb9-79cf6bde8555
# ╠═691eddff-f2eb-41a8-ab05-63afb46d15f2
