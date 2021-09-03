### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 239a91a6-f68a-11eb-14fd-0ba8d08b08f9
begin
	import Pkg
	Pkg.activate(Base.current_project())
	
	using AlgebraOfGraphics
	using CSV
	using CairoMakie
	using CCDReduction: fitscollection
	using Colors
	using DataFrames
	using DataFramesMeta
	using Dates
	using DelimitedFiles
	using Glob
	using ImageFiltering
	using KernelDensity
	using Latexify
	using LaTeXStrings
	using LombScargle
	using Measurements
	using Measurements: value, uncertainty
	using NaturalSort
	using OrderedCollections
	using Printf
	using Statistics
	using PlutoUI: TableOfContents, Select, Slider, as_svg, with_terminal
	using Unitful
	
	# Python setup
	#ENV["PYTHON"] = "~/miniconda3/envs/WASP-50b/bin/python"
	Pkg.build("PyCall")
	using PyCall
	
	##############
	# PLOT CONFIGS
	##############
	const FIG_TALL = (900, 1_200)
	const FIG_WIDE = (1_350, 800)
	
	set_aog_theme!()
	update_theme!(
		Theme(
			Axis = (xlabelsize=18, ylabelsize=18,),
			Label = (textsize=18,  padding=(0, 10, 0, 0)),
			Lines = (linewidth=3, cycle=Cycle([:color, :linestyle], covary=true)),
			Scatter = (linewidth=10,),
			fontsize = 18,
			rowgap = 0,
			colgap = 0,
		)
	)
	
	COLORS = Makie.wong_colors()
end

# ╔═╡ 3dfe9e7d-3d77-4b49-b25b-3e7049906d26
using NamedArrays

# ╔═╡ 0132b4ab-0447-4546-b412-ec598b20d21d
md"""
# Retrievals

[Intro text here]

$(TableOfContents(title="📖 Table of Contents"))
"""

# ╔═╡ 60dc161c-2aa2-4264-884d-6da3ead0e57b
base_dir = "./data/retrievals/all_WASP50/NofitR0"

# ╔═╡ d7ce97c1-82f2-46f1-a5ac-73e38e032fc8
fit_R0 = basename(base_dir)

# ╔═╡ 0f65d095-09af-44d2-907b-c30e2c16b609
species = ["K", "TiO", "Na_TiO", "K_TiO", "Na_K_TiO"]

# ╔═╡ 093156c7-9da7-4814-9260-5173f27fa497
model_names = OrderedDict(
	"clear" => "NoHet_FitP0_NoClouds_NoHaze_$(fit_R0)",
	"clear+cloud" => "NoHet_FitP0_Clouds_NoHaze_$(fit_R0)",
	"clear+haze" => "NoHet_FitP0_NoClouds_Haze_$(fit_R0)",
	"clear+cloud+haze" => "NoHet_FitP0_Clouds_Haze_$(fit_R0)",
	"clear+spot" => "Het_FitP0_NoClouds_NoHaze_$(fit_R0)",
	#"clear+spot+cloud" => "Het_FitP0_Clouds_NoHaze_$(fit_R0)",
	"clear+spot+haze" => "Het_FitP0_NoClouds_Haze_$(fit_R0)",
	"clear+spot+cloud+haze" => "Het_FitP0_Clouds_Haze_$(fit_R0)",
)

# ╔═╡ 41370a85-7abc-42ac-b82e-f6d739d8b5a8
md"""
## Table
"""

# ╔═╡ 0064b4ce-c41d-4b4e-bba4-ef4be0430edc
function print_results(x, y, arr, idx)
	x[idx[1]], y[idx[2]], arr[idx]
end

# ╔═╡ 869c4e1e-ef11-4048-bb83-6710ce0b3c8e
md"""
## Plot
"""

# ╔═╡ 812210c9-e294-4d61-bdf6-a03284199188
function plot_evidences(nm)
	arr = nm.array
	n_subgroups, n_groups = size(arr)
	group_labels = nm.dicts[2].keys
	subgroup_labels = nm.dicts[1].keys
	
	tbl = (
		x = vcat((fill(n, n_subgroups) for n ∈ 1:n_groups)...),
		height = value.(vcat((eachcol(arr) .|> copy)...)),
		grp = vcat((1:n_subgroups for _ ∈ 1:n_groups)...),
	)
	
	fig = Figure()
	ax = Axis(
		fig[1, 1],
		xticks = (1:n_groups, group_labels),
	)
	
	barplot!(ax, tbl.x, tbl.height;
		dodge = tbl.grp,
		color = COLORS[tbl.grp]
	)
	
	labels = String.(subgroup_labels)
	elements = [PolyElement(polycolor = COLORS[i]) for i in 1:length(labels)]
	
	Legend(fig[0, :], elements, labels, nbanks=4, tellheight=true)
	
	save("../../ACCESS_WASP-50b/figures/detrended/retrieval_evidences.pdf", fig)
	
	fig
end

# ╔═╡ db524678-9ee2-4934-b1bb-6a2f13bf0fa6
dirpath = 
"data/retrievals/" *
"spot_lower_fit_R0/"*
"WASP50_E1_NoHet_FitP0_NoClouds_NoHaze_fitR0_Na_K"

# ╔═╡ 38f37547-9663-464d-9983-4fd3bdbc79b6
retr = CSV.File("$(dirpath)/retr_Magellan_IMACS.txt")

# ╔═╡ b245f244-897e-4eb2-a6a6-cef7c85ca390
retr_model = CSV.File("$(dirpath)/retr_model.txt") 

# ╔═╡ 66ad752e-29a2-4d1d-8985-4bda58138e31
retr_model_sampled = CSV.File("$(dirpath)/retr_model_sampled_Magellan_IMACS.txt")

# ╔═╡ e801501c-a882-4f2d-bbc1-40028c1c91d8
begin
	fig = Figure()
	ax = Axis(fig[1, 1])
	
	errorbars!(ax, retr.wav, retr.flux, retr.flux_err)
	scatter!(ax, retr.wav, retr.flux)
	
	lines!(ax, retr_model.wav, retr_model.flux)
	scatter!(ax, retr_model_sampled.wav, retr_model_sampled.flux)
	
	xlims!(ax, 0.2, 1.1)
	
	fig
end

# ╔═╡ 41a233c7-5357-453c-b7ad-36fdf9f709cb
md"""
## Helper functions
"""

# ╔═╡ 44b3b8cd-4b83-4b27-a948-d1230489552f
begin
	py"""
	import pickle
	
	def load_pickle(fpath):
		with open(fpath, "rb") as f:
			data = pickle.load(f)
		return data
	"""
	load_pickle(s) = py"load_pickle"(s)
end;

# ╔═╡ daacda36-1fc9-411f-b101-82944863c9f3
cube = OrderedDict(
	sp => OrderedDict(
		model_name => load_pickle(
			"$(base_dir)/WASP50_E1_$(model_id)_$(sp)/retrieval.pkl"
		)
		for (model_name, model_id) ∈ model_names
	)
	for sp ∈ species
)

# ╔═╡ 65b51ff6-0991-491f-8945-dd889ffe71dd
begin
	n_species, n_models = length(species), length(model_names)
	evidences = NamedArray(
		Matrix{Measurement{Float64}}(undef, n_models, n_species),
		(String.(keys(model_names)), species),
		("Model", "Species")
	)
	
	for (sp, models) ∈ cube
		for (model, data) in models
			evidences[model, sp] .= data["lnZ"] ± data["lnZerr"]
		end
	end
	
	ΔlnZ = evidences .- minimum(evidences)
end

# ╔═╡ 6e24e7f4-61e5-470d-8ce0-399e0fe32e90
print_results(model_names.keys, species, ΔlnZ, argmin(ΔlnZ))

# ╔═╡ 750f830c-5818-4d8f-a673-f838a9d0da46
print_results(model_names.keys, species, ΔlnZ, argmax(ΔlnZ))

# ╔═╡ 8af2ffc6-b24d-46c3-b9f5-ecc81c61cd49
plot_evidences(ΔlnZ)

# ╔═╡ 1eff1230-2423-4ac3-8e9b-f4e7bcd0121b
md"""
## Notebook setup
"""

# ╔═╡ Cell order:
# ╟─0132b4ab-0447-4546-b412-ec598b20d21d
# ╠═60dc161c-2aa2-4264-884d-6da3ead0e57b
# ╠═d7ce97c1-82f2-46f1-a5ac-73e38e032fc8
# ╠═0f65d095-09af-44d2-907b-c30e2c16b609
# ╠═093156c7-9da7-4814-9260-5173f27fa497
# ╠═daacda36-1fc9-411f-b101-82944863c9f3
# ╟─41370a85-7abc-42ac-b82e-f6d739d8b5a8
# ╠═3dfe9e7d-3d77-4b49-b25b-3e7049906d26
# ╠═65b51ff6-0991-491f-8945-dd889ffe71dd
# ╠═6e24e7f4-61e5-470d-8ce0-399e0fe32e90
# ╠═750f830c-5818-4d8f-a673-f838a9d0da46
# ╠═0064b4ce-c41d-4b4e-bba4-ef4be0430edc
# ╟─869c4e1e-ef11-4048-bb83-6710ce0b3c8e
# ╠═8af2ffc6-b24d-46c3-b9f5-ecc81c61cd49
# ╠═812210c9-e294-4d61-bdf6-a03284199188
# ╠═db524678-9ee2-4934-b1bb-6a2f13bf0fa6
# ╠═38f37547-9663-464d-9983-4fd3bdbc79b6
# ╠═b245f244-897e-4eb2-a6a6-cef7c85ca390
# ╠═66ad752e-29a2-4d1d-8985-4bda58138e31
# ╠═e801501c-a882-4f2d-bbc1-40028c1c91d8
# ╟─41a233c7-5357-453c-b7ad-36fdf9f709cb
# ╠═44b3b8cd-4b83-4b27-a948-d1230489552f
# ╟─1eff1230-2423-4ac3-8e9b-f4e7bcd0121b
# ╠═239a91a6-f68a-11eb-14fd-0ba8d08b08f9
