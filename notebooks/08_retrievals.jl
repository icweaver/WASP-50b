### A Pluto.jl notebook ###
# v0.17.2

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
end

# ╔═╡ 0132b4ab-0447-4546-b412-ec598b20d21d
md"""
# Retrievals

!!! note "Data download"

```
rclone sync -P drive_ACCESS:papers/WASP-50b/data/retrievals data/retrievals
```

$(TableOfContents(title="📖 Table of Contents"))
"""

# ╔═╡ 60dc161c-2aa2-4264-884d-6da3ead0e57b
base_dir = "./data/retrievals/WASP50"

# ╔═╡ d7ce97c1-82f2-46f1-a5ac-73e38e032fc8
fit_R0 = "fitR0"

# ╔═╡ 093156c7-9da7-4814-9260-5173f27fa497
model_names = OrderedDict(
	"clear" => "NoHet_FitP0_NoClouds_NoHaze_$(fit_R0)",
	"clear+spot" => "Het_FitP0_NoClouds_NoHaze_$(fit_R0)",
	"clear+cloud" => "NoHet_FitP0_Clouds_NoHaze_$(fit_R0)",
	"clear+cloud+spot" => "Het_FitP0_Clouds_NoHaze_$(fit_R0)",
	"clear+haze" => "NoHet_FitP0_NoClouds_Haze_$(fit_R0)",
	#"clear+cloud+haze" => "NoHet_FitP0_Clouds_Haze_$(fit_R0)",
	"clear+haze+spot" => "Het_FitP0_NoClouds_Haze_$(fit_R0)",
	#"clear+spot+cloud+haze" => "Het_FitP0_Clouds_Haze_$(fit_R0)",
)

# ╔═╡ 0f65d095-09af-44d2-907b-c30e2c16b609
species = [
	"Na",
	"K",
	"TiO",
	"Na_K",
	"Na_TiO",
	"K_TiO",
	"Na_K_TiO",
	"CO",
	"H2O",
	# "NH3", # too low
	# "HCN", # too low
	# "CH4", # too low
	# "CO2", # too low
	# "FEH", # too low
]

# ╔═╡ 704fa634-eee0-4eef-aacf-f75f2b53f4d2
md"""
```
cube
├── Na
│   ├── clear
│   │   ├── retr
│   │   ├── retr_model
│   │   └── retr_model_sampled
│   ├── ...
│   ├── clear+spot
│   └── clear+haze+spot
├── K
├── ...
└── H2O

```
"""

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

# ╔═╡ 41370a85-7abc-42ac-b82e-f6d739d8b5a8
md"""
## Table
"""

# ╔═╡ 1c4fe72d-9872-4969-a62a-5163b5009bbb
md"""
## Retreived transmission spectra
"""

# ╔═╡ d9008477-bfae-4df3-9538-6994a639120e
retr_instr = CSV.read(
		"$(base_dir)/../retr_Magellan_IMACS.txt", DataFrame;
		header = [:wav, :flux, :wav_d, :wav_u, :flux_err],
		comment = "#",
	)

# ╔═╡ db524678-9ee2-4934-b1bb-6a2f13bf0fa6
function get_retr_model(cube, sp, model)
	retr_model = cube[sp][model]["retr_model"]
	retr_model_sampled = cube[sp][model]["retr_model_sampled"]
	return retr_model, retr_model_sampled
end

# ╔═╡ cc011a66-37bd-4543-9a58-b11e1f785e52
function retrieval!(ax, model0, model_sampled; color=:blue, label="")
	model = @rsubset model0 0.5 < :wav < 1.0
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

# ╔═╡ 00a0f9c4-cd4d-4ae2-80b7-0c044239a571
function plot_retrieval!(ax, cube, sp, model; color=:blue)
	retr_model, retr_model_sampled = get_retr_model(cube, sp, model)
	label = replace(sp, "_"=>"+") * " ($(model))"
	retrieval!(ax, retr_model, retr_model_sampled, color=color, label=label)
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

# ╔═╡ a7c68d25-a799-421b-9799-38837fa8a188
begin
	cube = OrderedDict()
	for sp ∈ species
		cube[sp] = OrderedDict()
		for (model_name, model_id) ∈ model_names
		cube[sp][model_name] = Dict()
		dirpath = "$(base_dir)/WASP50_E1_$(model_id)_$(sp)"
			cube[sp][model_name]["retr"] = load_pickle("$(dirpath)/retrieval.pkl")
			cube[sp][model_name]["retr_model"] = CSV.read(
				"$(dirpath)/retr_model.txt", DataFrame;
				header = [:wav, :flux, :flux_d, :flux_u],
				comment = "#",
			)
			cube[sp][model_name]["retr_model_sampled"] = CSV.read(
				"$(dirpath)/retr_model_sampled_Magellan_IMACS.txt", DataFrame;
				header = [:wav, :flux],
				comment = "#",
			)
		end
	end
	cube
end

# ╔═╡ 65b51ff6-0991-491f-8945-dd889ffe71dd
begin
	n_species, n_models = length(species), length(model_names)
	row_labels = model_names.keys
	col_labels = species
	
	df_evidences = DataFrame(
		Species=String[], Model=String[], lnZ=Float64[], lnZ_err=Float64[]
	)
	for (sp, model_dicts) ∈ cube
		for (model, model_dict) ∈ model_dicts
			data = model_dict["retr"]
			push!(df_evidences, (sp, model, data["lnZ"], data["lnZerr"]))
		end
	end

	@transform! df_evidences :lnZ = :lnZ .- minimum(:lnZ)
end

# ╔═╡ 1eff1230-2423-4ac3-8e9b-f4e7bcd0121b
md"""
## Notebook setup
"""

# ╔═╡ e43f1834-73dd-4859-b847-f4c552561897
begin
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
	
	COLORS = let
		pal = Makie.ColorSchemes.Paired_8 |> reverse
		[pal[7:8] ; pal[5:6] ; pal[1:2]]
	end
end

# ╔═╡ df43608e-7026-45ae-b87b-d7e0b6cea89c
let
	sort_order = sorter(model_names.keys)
	
	plt = data(df_evidences) *
		mapping(:Species => sorter(species), :lnZ => "ΔlnZ";
			dodge = :Model => sort_order,
			color = :Model => sort_order,
		) *
		visual(BarPlot)
	
	draw(plt;
		axis = (; limits=(nothing, nothing, 0, 3.5)),
		legend = (
			position = :top,
			nbanks = 2,
			orientation = :horizontal,
			titlevisible = false,
		),
		palettes = (; color=COLORS),
	)
end

# ╔═╡ 812210c9-e294-4d61-bdf6-a03284199188
function plot_evidences(nm)
	arr = nm.array
	n_subgroups, n_groups = size(arr)
	group_labels = replace.(nm.dicts[2].keys, "_"=>"+")
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
		xlabel = "Species",
		ylabel = "Δln Z",
	)
	
	barplot!(ax, tbl.x, tbl.height;
		dodge = tbl.grp,
		color = COLORS[tbl.grp]
	)
	
	labels = String.(subgroup_labels)
	elements = [PolyElement(polycolor = COLORS[i]) for i in 1:length(labels)]
	
	axislegend(ax, elements, labels, nbanks=2, orientation=:horizontal, titlevisible=true, legendtitle="hey")
	
	ylims!(ax, 0, 3.5)
	
	# path = "../../ACCESS_WASP-50b/figures/detrended"
	# mkpath(path)
	# save("$(path)/retrieval_evidences.png", fig)
	
	fig
end

# ╔═╡ e801501c-a882-4f2d-bbc1-40028c1c91d8
let
	fig = Figure(resolution=(800, 500))
	ax = Axis(fig[1, 1], xlabel="Wavelength (μm)", ylabel="Transit depth (ppm)")

	plot_retrieval!(ax, cube, "Na_TiO", "clear", color=COLORS[1])
	plot_retrieval!(ax, cube, "TiO", "clear+haze+spot", color=COLORS[6])

	# Instrument
	errorbars!(ax, retr_instr.wav, retr_instr.flux, retr_instr.flux_err)
	scatter!(ax, retr_instr.wav, retr_instr.flux;
		markersize = 15,
		color = :white,
		strokewidth=1.5,
		label = "IMACS + LDSS3"
	)

	axislegend(orientation=:horizontal)

	# path = "../../ACCESS_WASP-50b/figures/retrievals"
	# mkpath(path)
	# save("$(path)/retr.png", fig)
	
	fig
end

# ╔═╡ Cell order:
# ╠═0132b4ab-0447-4546-b412-ec598b20d21d
# ╠═60dc161c-2aa2-4264-884d-6da3ead0e57b
# ╠═d7ce97c1-82f2-46f1-a5ac-73e38e032fc8
# ╠═093156c7-9da7-4814-9260-5173f27fa497
# ╠═0f65d095-09af-44d2-907b-c30e2c16b609
# ╟─704fa634-eee0-4eef-aacf-f75f2b53f4d2
# ╠═a7c68d25-a799-421b-9799-38837fa8a188
# ╠═7b714c1e-2e3d-453f-a342-81df8283de5c
# ╟─41370a85-7abc-42ac-b82e-f6d739d8b5a8
# ╠═65b51ff6-0991-491f-8945-dd889ffe71dd
# ╠═df43608e-7026-45ae-b87b-d7e0b6cea89c
# ╟─1c4fe72d-9872-4969-a62a-5163b5009bbb
# ╠═812210c9-e294-4d61-bdf6-a03284199188
# ╠═d9008477-bfae-4df3-9538-6994a639120e
# ╠═e801501c-a882-4f2d-bbc1-40028c1c91d8
# ╠═00a0f9c4-cd4d-4ae2-80b7-0c044239a571
# ╠═db524678-9ee2-4934-b1bb-6a2f13bf0fa6
# ╠═cc011a66-37bd-4543-9a58-b11e1f785e52
# ╟─41a233c7-5357-453c-b7ad-36fdf9f709cb
# ╠═44b3b8cd-4b83-4b27-a948-d1230489552f
# ╟─1eff1230-2423-4ac3-8e9b-f4e7bcd0121b
# ╠═e43f1834-73dd-4859-b847-f4c552561897
# ╠═239a91a6-f68a-11eb-14fd-0ba8d08b08f9
