### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ 63754778-2464-423c-a4ed-9f202be4db15
using CSV

# ╔═╡ 691eddff-f2eb-41a8-ab05-63afb46d15f2
begin
	import PlutoUI as pl
	import JSON
	using CairoMakie
	using Colors
	using Glob
	using PyCall
	using Statistics, KernelDensity
end

# ╔═╡ 506eeeb2-e56d-436b-91b8-605e52201563
md"""
# Detrended white light curves

Gonna do some stuff

$(pl.TableOfContents())
"""

# ╔═╡ a8cf11e2-796e-45ff-bdc9-e273b927700e
md"""
## Transit curves ⚪
"""

# ╔═╡ 3f0f5777-00f1-443d-8ced-d901550010d3
const DATA_DIR = "data/detrended_wlcs/out_l/WASP50"

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
end

# ╔═╡ b3e47bc9-21dd-4e32-a03e-25ec32f9f4b9
# data_WLC = CSV.File(
# 	"$(DATA_DIR)/w50_131219/white-light/results.dat",
# 	comment = "#",
# 	delim = ' ',
# 	ignorerepeated=true,
	
# )

# ╔═╡ 68ec4343-5f6c-4dfd-90b5-6393b4c819b9
md"""
## Corner plots 📐
"""

# ╔═╡ b7eac49f-f140-43ca-876a-e480b593e885
const PARAMS = ["p", "t0", "P", "rho"] #, "aR", "inc", "b", "q1"]

# ╔═╡ ddd9a95b-735a-4995-8893-542128bb56d6
truths = JSON.parsefile("$(DATA_DIR)/truth.json")

# ╔═╡ 931ce3d5-c4ed-496c-883b-d7ee33e957cc
function adjust_dict(dict, params, truths)
	d = filter!(p -> p.first ∈ params, dict)
	for param in params
		d[param] = d[param] .- truths[param]["truth"][1]
	end
	
	return d
end

# ╔═╡ 2191791b-df62-4f1b-88bf-060cc47896b2
cubes = Dict(
	"Transit $i" => Dict(
		
		"models" => load_npz(fpath_model, allow_pickle = true),
		
		"samples" => adjust_dict(load_pickle(fpath_sample), PARAMS, truths),
		
		"results" => CSV.File(
			fpath_result,
			comment = "#",
			delim = ' ',			
			ignorerepeated=true,
		)
	)
	
	for (i, (fpath_model, fpath_result, fpath_sample)) in enumerate(
				zip(
					sort(glob("$(DATA_DIR)/w50_*/white-light/BMA_WLC.npy")),
					sort(glob("$(DATA_DIR)/w50_*/white-light/results.dat")),
					sort(glob("$(DATA_DIR)/w50_*/white-light/BMA_posteriors.pkl"))
				)
		)
)

# ╔═╡ 831c5bbd-1b55-4b26-99f0-b9ae1959abef
# cubes_dist = Dict(
# 	"Transit $i" => adjust_dict(load_pickle(fpath), PARAMS, truths)
# 	for (i, fpath) in enumerate(
# 		glob("$(DATA_DIR)/w50_*/white-light/BMA_posteriors.pkl") |> sort
# 	)
# )

# ╔═╡ d81d5ffd-f51d-4b0c-bf33-0bef9899d549
truth_vals(d, param, truths) = (
	d.Value[data_WLC.Variable .== param][1] - truths[param]["truth"][1],
	d.SigmaUp[data_WLC.Variable .== param][1],
	d.SigmaDown[data_WLC.Variable .== param][1]
)

# ╔═╡ 6fcd1377-8364-45a3-9ff6-89d61df1ef42
levels(A, n) = reverse(
	range(maximum(A), step=-maximum(A)/(n+1), length=(n+1))
)

# ╔═╡ 2cbc6ddb-210e-41e8-b745-5c41eba4e778
function plot_corner!(fig, cube, params; n_levels=4, color=:blue)
	for (j, p1) in enumerate(params), (i, p2) in enumerate(params)
		if i == j
			density!(
			fig[i, i],
			cube[p1],
			color = (color, 0.125),
			#strokewidth = 0,
			strokewidth = 3,
			strokecolor = color,
			strokearound = true,
		)
		end
		if i > j
			Z = kde((cube[p1], cube[p2]), npoints=(2^4, 2^4),)
			contourf!(
				fig[i, j],
				Z,
				levels = levels(Z.density, n_levels),
				#levels = levels(Z.weights, 5),
				colormap = cgrad(range(colorant"white", color), alpha=0.5),
			)
			contour!(
				fig[i, j],
				Z,
				levels = levels(Z.density, n_levels),
				#levels = levels(Z.weights, 5),
				color = color,
				linewidth = 3,
			)
		end
	end
end

# ╔═╡ ce4e351d-621a-4031-ac31-e38ddccfa3ed
let
	fig = Figure()

	ax, l = lines(fig[1, 1], [1, 2, 3])
	
	#current_axis()
	
	l
end

# ╔═╡ 900bcec9-00b6-4c86-8782-0f0b8473620c
to_value

# ╔═╡ 940ebaf2-659a-4319-bbe6-e0290752f1fb
const COLORS =  [
	# "#fdbf6f",  # Yellow
	colorant"#f7ad4d",  # Yellow
	colorant"#ff7f00",  # Orange
	# "#a6cee3",  # Cyan
	colorant"#5daed9",  # Cyan
	# "#75bfe6",  # Cyan
	# "#1f78b4",  # Blue
	colorant"#126399",  # Blue
	colorant"plum",
	colorant"#956cb4",  # Purple
	colorant"mediumaquamarine",
	colorant"#029e73",  # Green
	colorant"slategray",
]

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
begin
	fig = Figure(resolution=(1_600, 1_600))
	
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
	
	# Axis(fig[0, 1])
	# for j in 2:n_params
	# 	Axis(fig[1, j])
	# end
	
	# lines!(fig[1, 1], [-1, 0, 1], [1, 1, 1])
	# text!(
	# 	fig[1, 1],
	# 	"ok",
	# 	position = Point(0, 2),
	# 	textsize = 0.2,
	# 	align = (:center, :baseline),
	# )

	axs = reshape(copy(fig.content), n_params, n_params)
	
	#lines!(fig[3, 2], [1, 2, 3])
	# for i in 1:n_params
	# 	v, vu, vd = truth_vals(data_WLC, PARAMS[i], truths)
	# 	poly!(fig[i+1, i], Rect(0, v, maximum((vu, vd)), 1.5e6), color=(:red, 0.5))
	# end
	
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
	
	fig
	
	
end

# ╔═╡ baeadfce-535a-46c3-8cb9-79cf6bde8555
md"""
## Packages 📦
"""

# ╔═╡ Cell order:
# ╟─506eeeb2-e56d-436b-91b8-605e52201563
# ╟─a8cf11e2-796e-45ff-bdc9-e273b927700e
# ╠═3f0f5777-00f1-443d-8ced-d901550010d3
# ╠═2191791b-df62-4f1b-88bf-060cc47896b2
# ╠═39dbca86-a4b9-11eb-1c64-9ddf1a9990ab
# ╠═b3e47bc9-21dd-4e32-a03e-25ec32f9f4b9
# ╠═4be0d7b7-2ea5-4c4d-92b9-1f8109014e12
# ╟─68ec4343-5f6c-4dfd-90b5-6393b4c819b9
# ╠═b7eac49f-f140-43ca-876a-e480b593e885
# ╠═ddd9a95b-735a-4995-8893-542128bb56d6
# ╠═931ce3d5-c4ed-496c-883b-d7ee33e957cc
# ╠═831c5bbd-1b55-4b26-99f0-b9ae1959abef
# ╠═63754778-2464-423c-a4ed-9f202be4db15
# ╠═d81d5ffd-f51d-4b0c-bf33-0bef9899d549
# ╠═d5ff9b30-00dd-41d3-9adf-ff7905d71ae8
# ╠═2cbc6ddb-210e-41e8-b745-5c41eba4e778
# ╠═6fcd1377-8364-45a3-9ff6-89d61df1ef42
# ╠═ce4e351d-621a-4031-ac31-e38ddccfa3ed
# ╠═900bcec9-00b6-4c86-8782-0f0b8473620c
# ╠═940ebaf2-659a-4319-bbe6-e0290752f1fb
# ╟─baeadfce-535a-46c3-8cb9-79cf6bde8555
# ╠═691eddff-f2eb-41a8-ab05-63afb46d15f2
