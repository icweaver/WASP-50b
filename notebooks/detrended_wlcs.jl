### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ 691eddff-f2eb-41a8-ab05-63afb46d15f2
begin
	import PlutoUI as pl
	import JSON
	using CairoMakie
	using Colors
	using CSV
	using Glob
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

# ╔═╡ a8cf11e2-796e-45ff-bdc9-e273b927700e
md"""
## Transit curves ⚪
"""

# ╔═╡ ae82d3c1-3912-4a5e-85f5-6383af42291e
md"""
Plotting the data from the `models` cube returns the following detrended white light curves:
"""

# ╔═╡ 68ec4343-5f6c-4dfd-90b5-6393b4c819b9
md"""
## Corner plots 📐
"""

# ╔═╡ e452d2b1-1010-4ce3-8d32-9e9f1d0dfa0b
md"""
The `samples` cube returns the following corner plot for the fitted `PARAMS` below for each night:
"""

# ╔═╡ b7eac49f-f140-43ca-876a-e480b593e885
const PARAMS = ["p", "t0", "P", "rho", "aR", "inc", "b", "q1"]

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
				#strokewidth = 0,
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

# ╔═╡ d81d5ffd-f51d-4b0c-bf33-0bef9899d549
truth_vals(d, param, truths) = (
	d.Value[data_WLC.Variable .== param][1] - truths[param]["truth"][1],
	d.SigmaUp[data_WLC.Variable .== param][1],
	d.SigmaDown[data_WLC.Variable .== param][1]
)

# ╔═╡ 931ce3d5-c4ed-496c-883b-d7ee33e957cc
function adjust_dict(dict, params, truths)
	d = filter!(p -> p.first ∈ params, dict)
	for param in params
		truth_val = truths[param]["truth"][1]
		truth_val_err = maximum((truths[param]["truth"][2:3])) 
		d[param] = @. ((d[param] - truth_val) / truth_val_err)
	end
	
	return d
end

# ╔═╡ 2191791b-df62-4f1b-88bf-060cc47896b2
cubes = Dict(
	"Transit $i" => Dict(
		
		"samples" => adjust_dict(load_pickle(fpath_sample), PARAMS, truths),
		
		"models" => load_npz(fpath_model, allow_pickle=true),
		
		"results" => CSV.File(
			fpath_result,
			comment = "#",
			delim = ' ',			
			ignorerepeated=true,
		)
	)
	
	for (i, (fpath_sample, fpath_model, fpath_result)) in enumerate(zip(
		sort(glob("$(DATA_DIR)/w50_*/white-light/BMA_posteriors.pkl")),
		sort(glob("$(DATA_DIR)/w50_*/white-light/BMA_WLC.npy")),
		sort(glob("$(DATA_DIR)/w50_*/white-light/results.dat")),
	))
)

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
	fig = Figure(resolution=(1_400, 1_400))
	
	# Number of fitted parameters
	n_params = length(PARAMS)
	
	for j in 1:n_params, i in 1:n_params
		ax = Axis(fig[i, j], axis=(aspect=1, xticklabelrotation=π/4),)
		ax.xticklabelrotation = π/4
		ax.aspect = 1.0
		j > i && (hidedecorations!(ax); hidespines!(ax))
		i == j && i != n_params && (hidedecorations!(ax))
	end
			
	for (i, (transits, cube)) in enumerate(cubes)
		plot_corner!(fig, cube["samples"], PARAMS, color=COLORS[i])
	end
	
	axs = reshape(copy(fig.content), n_params, n_params)
	
	[linkxaxes!(reverse(axs[:, j])...) for j in 1:n_params]
		
	hidexdecorations!.(axs[begin:end-1, :])
	hideydecorations!.(axs[:, begin+1:end])
	
	for j in 1:n_params, i in 1:n_params
		i > j && (
			axs[i, j].xgridvisible = true;
			axs[i, j].ygridvisible = true
		)
		if i == j
			axs[i, j].xgridvisible = true
		end
	end
	
	for (j, param) in enumerate(PARAMS)
		axs[end, j].xlabel = "Δ"*truths[param]["symbol"]
		axs[end, j].xlabelsize = 26
		axs[j, begin].ylabel = "Δ"*truths[param]["symbol"]
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
# ╠═579e62da-7ffb-4639-bd73-3826ade1cfa2
# ╠═583b3377-f4f6-4170-8658-d3ba17a5b86d
# ╠═39dbca86-a4b9-11eb-1c64-9ddf1a9990ab
# ╟─a8cf11e2-796e-45ff-bdc9-e273b927700e
# ╟─ae82d3c1-3912-4a5e-85f5-6383af42291e
# ╠═4be0d7b7-2ea5-4c4d-92b9-1f8109014e12
# ╟─68ec4343-5f6c-4dfd-90b5-6393b4c819b9
# ╟─e452d2b1-1010-4ce3-8d32-9e9f1d0dfa0b
# ╠═b7eac49f-f140-43ca-876a-e480b593e885
# ╠═d5ff9b30-00dd-41d3-9adf-ff7905d71ae8
# ╟─c5e10e47-ec64-4911-8107-487d1ef3f134
# ╠═2cbc6ddb-210e-41e8-b745-5c41eba4e778
# ╠═6fcd1377-8364-45a3-9ff6-89d61df1ef42
# ╟─82a23101-9e1f-4eae-b529-e750a44c98b1
# ╠═ddd9a95b-735a-4995-8893-542128bb56d6
# ╠═d81d5ffd-f51d-4b0c-bf33-0bef9899d549
# ╠═931ce3d5-c4ed-496c-883b-d7ee33e957cc
# ╟─30ae3744-0e7e-4c16-b91a-91eb518fba5b
# ╠═940ebaf2-659a-4319-bbe6-e0290752f1fb
# ╟─baeadfce-535a-46c3-8cb9-79cf6bde8555
# ╠═691eddff-f2eb-41a8-ab05-63afb46d15f2
