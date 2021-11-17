### A Pluto.jl notebook ###
# v0.17.1

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
	using DataFrameMacros
	using Dates
	using DelimitedFiles
	using Glob
	using ImageFiltering
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
	
	import CairoMakie.Makie.KernelDensity: kde
	
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
base_dir = "./data/retrievals/all_WASP50"

# ╔═╡ d7ce97c1-82f2-46f1-a5ac-73e38e032fc8
fit_R0 = "fitR0"

# ╔═╡ 093156c7-9da7-4814-9260-5173f27fa497
model_names = OrderedDict(
	"clear" => "NoHet_FitP0_NoClouds_NoHaze_$(fit_R0)",
	"clear+cloud" => "NoHet_FitP0_Clouds_NoHaze_$(fit_R0)",
	"clear+haze" => "NoHet_FitP0_NoClouds_Haze_$(fit_R0)",
	#"clear+cloud+haze" => "NoHet_FitP0_Clouds_Haze_$(fit_R0)",
	"clear+spot" => "Het_FitP0_NoClouds_NoHaze_$(fit_R0)",
	"clear+spot+cloud" => "Het_FitP0_Clouds_NoHaze_$(fit_R0)",
	"clear+spot+haze" => "Het_FitP0_NoClouds_Haze_$(fit_R0)",
	#"clear+spot+cloud+haze" => "Het_FitP0_Clouds_Haze_$(fit_R0)",
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

# ╔═╡ 0f65d095-09af-44d2-907b-c30e2c16b609
species = [
	"Na",
	"K",
	"TiO",
	"Na_TiO",
	"K_TiO",
	"Na_K_TiO",
	#"CO",
	#"H2O",
	#"NH3", # too low
	#"HCN", # too low
	#"CH4", # too low
	#"CO2", # too low
	#"FEH", # too low
]

# ╔═╡ 7b714c1e-2e3d-453f-a342-81df8283de5c
# Check if missing files
with_terminal() do
	for sp ∈ species
		for (model_name, model_id) ∈ model_names
			!isfile("$(base_dir)/WASP50_E1_$(model_id)_$(sp)/retrieval.pkl") &&
				println("WASP50_E1_$(model_id)_$(sp)")
				#println("$(base_dir)/WASP50_E1_$(model_id)_$(sp)/retrieval.pkl")
		end
	end
end

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
	
	axislegend(ax, elements, labels, nbanks=3)
	
	ylims!(ax, 0, 5.5)
	
	path = "../../ACCESS_WASP-50b/figures/detrended"
	mkpath(path)
	save("$(path)/retrieval_evidences.png", fig)
	
	fig
end

# ╔═╡ db524678-9ee2-4934-b1bb-6a2f13bf0fa6
begin
	dirpath_Na = "data/retrievals/all_WASP50/WASP50_E1_NoHet_FitP0_NoClouds_NoHaze_fitR0_Na"
	
	retr = CSV.read(
		"$(dirpath_Na)/retr_Magellan_IMACS.txt", DataFrame;
		header = [:wav, :flux, :wav_d, :wav_u, :flux_err],
		comment = "#",
	)

	retr_model_all_Na = CSV.read("$(dirpath_Na)/retr_model.txt", DataFrame;
	header = [:wav, :flux, :flux_d, :flux_u],
	comment = "#",
) 
	retr_model_sampled_Na = CSV.read(
	"$(dirpath_Na)/retr_model_sampled_Magellan_IMACS.txt", DataFrame;
	header = [:wav, :flux],
	comment = "#",
)
end

# ╔═╡ d97ec176-af5a-4f95-b891-7ceb5ae1b3e0
begin
	dirpath_Na_TiO = "data/retrievals/all_WASP50/WASP50_E1_NoHet_FitP0_NoClouds_NoHaze_fitR0_Na_TiO"
	
	retr_Na_TiO = CSV.read(
		"$(dirpath_Na_TiO)/retr_Magellan_IMACS.txt", DataFrame;
		header = [:wav, :flux, :wav_d, :wav_u, :flux_err],
		comment = "#",
	)

	retr_model_all_Na_TiO = CSV.read("$(dirpath_Na_TiO)/retr_model.txt", DataFrame;
	header = [:wav, :flux, :flux_d, :flux_u],
	comment = "#",
) 
	retr_model_sampled_Na_TiO = CSV.read(
	"$(dirpath_Na_TiO)/retr_model_sampled_Magellan_IMACS.txt", DataFrame;
	header = [:wav, :flux],
	comment = "#",
)
end

# ╔═╡ cc011a66-37bd-4543-9a58-b11e1f785e52
function retrieval!(ax, model, model_sampled; color=:blue, label="")
	lines!(ax, model.wav, model.flux, color=color, label=label)
	scatter!(ax, model_sampled.wav, model_sampled.flux;
		marker = :rect,
		markersize = 15,
		color = color,
		strokewidth = 1.5,
		strokecolor = color,
	)
	band!(ax, model.wav, model.flux_d, model.flux_u;
		color = (color, 0.25),
	)
end

# ╔═╡ e801501c-a882-4f2d-bbc1-40028c1c91d8
let
	fig = Figure(resolution=(800, 500))
	ax = Axis(fig[1, 1], xlabel="Wavelength (μm)", ylabel="Transit depth (ppm)")

	model = @subset retr_model_all_Na 0.5 < :wav < 1.0
	retrieval!(ax,  model, retr_model_sampled_Na, color=COLORS[2], label="Na")

	model = @subset retr_model_all_Na_TiO 0.5 < :wav < 1.0
	retrieval!(ax,  model, retr_model_sampled_Na_TiO, color=COLORS[1], label="Na_TiO")

	errorbars!(ax, retr.wav, retr.flux, retr.flux_err)
	scatter!(ax, retr.wav, retr.flux;
		markersize = 15,
		color = :white,
		strokewidth=1.5,
		label = "IMACS + LDSS3"
	)

	axislegend(orientation=:horizontal)

	path = "../../ACCESS_WASP-50b/figures/retrievals"
	mkpath(path)
	save("$(path)/retr.png", fig)
	
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
# ╠═093156c7-9da7-4814-9260-5173f27fa497
# ╠═daacda36-1fc9-411f-b101-82944863c9f3
# ╠═7b714c1e-2e3d-453f-a342-81df8283de5c
# ╟─41370a85-7abc-42ac-b82e-f6d739d8b5a8
# ╠═3dfe9e7d-3d77-4b49-b25b-3e7049906d26
# ╠═65b51ff6-0991-491f-8945-dd889ffe71dd
# ╠═6e24e7f4-61e5-470d-8ce0-399e0fe32e90
# ╠═750f830c-5818-4d8f-a673-f838a9d0da46
# ╠═0064b4ce-c41d-4b4e-bba4-ef4be0430edc
# ╟─869c4e1e-ef11-4048-bb83-6710ce0b3c8e
# ╠═0f65d095-09af-44d2-907b-c30e2c16b609
# ╠═8af2ffc6-b24d-46c3-b9f5-ecc81c61cd49
# ╠═812210c9-e294-4d61-bdf6-a03284199188
# ╠═db524678-9ee2-4934-b1bb-6a2f13bf0fa6
# ╠═d97ec176-af5a-4f95-b891-7ceb5ae1b3e0
# ╠═e801501c-a882-4f2d-bbc1-40028c1c91d8
# ╠═cc011a66-37bd-4543-9a58-b11e1f785e52
# ╟─41a233c7-5357-453c-b7ad-36fdf9f709cb
# ╠═44b3b8cd-4b83-4b27-a948-d1230489552f
# ╟─1eff1230-2423-4ac3-8e9b-f4e7bcd0121b
# ╠═239a91a6-f68a-11eb-14fd-0ba8d08b08f9
