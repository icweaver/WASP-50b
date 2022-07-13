### A Pluto.jl notebook ###
# v0.19.9

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

# ╔═╡ 239a91a6-f68a-11eb-14fd-0ba8d08b08f9
begin
	import Pkg
	Pkg.activate(Base.current_project())

	using PlutoUI
	import MarkdownLiteral: @mdx
	using AlgebraOfGraphics, CairoMakie
	using CSV, DataFrames, DataFramesMeta, DelimitedFiles, Glob, OrderedCollections
	using ImageFiltering, LombScargle, Measurements, Statistics
	using Latexify, Printf
	using Dates, NaturalSort
	using CondaPkg
	using PythonCall
end

# ╔═╡ a21aad0b-5998-420e-b437-d7ad262d0fe4
begin
	const DATA_DIR = "data/retrievals"
	const FIG_DIR = "figures/retrievals"
	TableOfContents()
end

# ╔═╡ 0132b4ab-0447-4546-b412-ec598b20d21d
@mdx """
# Retrievals

!!! note "Data download"
	```
	rclone sync -P ACCESS_box:WASP-50b/$(DATA_DIR) $(DATA_DIR)
	```
	* [Direct link](https://app.box.com/s/hmqutmcs98rkip3aa4qyt68g3k1vk2du)
"""

# ╔═╡ 60dc161c-2aa2-4264-884d-6da3ead0e57b
@bind base_dir Select(glob("data/retrievals/WASP50*"))

# ╔═╡ d7ce97c1-82f2-46f1-a5ac-73e38e032fc8
@bind fit_R0 Select(["NofitR0", "fitR0"])

# ╔═╡ 589afac8-0ea5-4962-b52b-7f035e91cf44
fname_suff = let
	b = basename(base_dir)
	label = occursin("offs", b) ? "offs" : "no_offs"
	label *= "_$(fit_R0)"
	label *= "_$(b)"
end

# ╔═╡ 2c12ec4d-1184-4755-8bd8-0d7cd59fa205
md"""
## Evidences
"""

# ╔═╡ 9b43cef3-4583-4463-89ac-5f8678d1e843
help_attributes(VLines)

# ╔═╡ 093156c7-9da7-4814-9260-5173f27fa497
model_names = OrderedDict(
	"clear" => "NoHet_FitP0_NoClouds_NoHaze_$(fit_R0)",
	"spot" => "Het_FitP0_NoClouds_NoHaze_$(fit_R0)",
	"cloud" => "NoHet_FitP0_Clouds_NoHaze_$(fit_R0)",
	"cloud+spot" => "Het_FitP0_Clouds_NoHaze_$(fit_R0)",
	"haze" => "NoHet_FitP0_NoClouds_Haze_$(fit_R0)",
	#"clear+cloud+haze" => "NoHet_FitP0_Clouds_Haze_$(fit_R0)",
	"haze+spot" => "Het_FitP0_NoClouds_Haze_$(fit_R0)",
	#"clear+spot+cloud+haze" => "Het_FitP0_Clouds_Haze_$(fit_R0)",
	"flat" => "NoHet_FlatLine",
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
	# "CO",
	"H2O",
	# "NH3", # too low
	# "HCN", # too low
	# "CH4", # too low
	# "CO2", # too low
	# "FEH", # too low
]

# ╔═╡ 704fa634-eee0-4eef-aacf-f75f2b53f4d2
@mdx """
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
@with_terminal begin
	i = 0
	for sp ∈ species
		for (model_name, model_id) ∈ model_names
			if !isfile("$(base_dir)/WASP50_E1_$(model_id)_$(sp)/retrieval.pkl")
				println("WASP50_E1_$(model_id)_$(sp)")
				i += 1
				#println("$(base_dir)/WASP50_E1_$(model_id)_$(sp)/retrieval.pkl")
			end
		end
	end
	println("\n$(i) incomplete directories")
end

# ╔═╡ a0094689-a9d5-4810-baba-bd7a96c27839
dashplus(x) = replace(x, '_' => '+')

# ╔═╡ 1c4fe72d-9872-4969-a62a-5163b5009bbb
@mdx """
## Retrieved transmission spectra 🐕
"""

# ╔═╡ 5569b57c-0585-4300-927b-5d089dde0f43
function plot_instr!(ax, fpath; color=:blue, label="add label")
	retr_instr = CSV.read(
		"$(base_dir)/$(fpath)", DataFrame;
		header = [:wav, :flux, :wav_d, :wav_u, :flux_err],
		comment = "#",
	)
	wav = retr_instr.wav .* 10_000
	errorbars!(ax, wav, retr_instr.flux, retr_instr.flux_err)
	scatter!(ax, wav, retr_instr.flux;
		markersize = 15,
		color,
		strokewidth = 1.5,
		label,
	)
end

# ╔═╡ db524678-9ee2-4934-b1bb-6a2f13bf0fa6
function get_retr_model(cube, sp, model)
	retr_model = cube[sp][model]["retr_model"]
	retr_model_sampled = cube[sp][model]["retr_model_sampled"]
	return retr_model, retr_model_sampled
end

# ╔═╡ cc011a66-37bd-4543-9a58-b11e1f785e52
function retrieval!(ax, model0, model_sampled;
	color=:blue, linewidth=3, label="", show_scatter=true
)
	model = @rsubset model0 :wav < 1.0
	wav = model.wav .* 10_000
	wav_sampled = model_sampled.wav .* 10_000
	lines!(ax, wav, model.flux; color, linewidth, label)
	if show_scatter
		scatter!(ax, wav_sampled, model_sampled.flux;
			marker = :rect,
			markersize = 15,
			color,
			strokewidth = 1.5,
			strokecolor = color,
		)
	end
	band!(ax, wav, model.flux_d, model.flux_u;
		color = (color, 0.25),
	)
end

# ╔═╡ 00a0f9c4-cd4d-4ae2-80b7-0c044239a571
function plot_retrieval!(ax, cube, sp, model; color=:blue, linewidth=3, show_scatter=true)
	retr_model, retr_model_sampled = get_retr_model(cube, sp, model)
	if model == "flat"
		label = "(flat)"
	else
		label = dashplus(sp) * " ($(model))"
	end
	retrieval!(ax, retr_model, retr_model_sampled; color, linewidth, label, show_scatter)
end

# ╔═╡ 0f23e0d6-177d-4bf6-9660-f2c376b3146b
md"""
## Posterior distributions
"""

# ╔═╡ 1eff1230-2423-4ac3-8e9b-f4e7bcd0121b
@mdx """
## Notebook setup 🔧
"""

# ╔═╡ eab74923-a084-468c-9b0d-c2cc98a23913
function savefig(fig, fpath)
	mkpath(dirname(fpath))
    save(fpath, fig, pt_per_unit=1)
	@info "Saved to: $(fpath)"
end

# ╔═╡ 44b3b8cd-4b83-4b27-a948-d1230489552f
begin
	@pyexec """
	global pickle, load_pickle
	import pickle

	def load_pickle(fpath):
		with open(fpath, "rb") as f:
			data = pickle.load(f)
		return data
	"""
	load_pickle(s) = @pyeval("load_pickle")(s)
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
			"$(dirpath)/retr_model_sampled_IMACS+LDSS3C.txt", DataFrame;
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
			lnZ = pyconvert(Float64, data["lnZ"])
			lnZerr = pyconvert(Float64, data["lnZerr"])
			push!(df_evidences, (sp, model, lnZ, lnZerr))
		end
	end

	@transform! df_evidences begin
		:ΔlnZ = :lnZ .- minimum(:lnZ)
		:Species = dashplus.(:Species)
		:ΔlnZ_m = (:lnZ .± :lnZ_err) .- minimum(:lnZ .± :lnZ_err)
	end

	# df_evidences.ΔlnZ .-= df_evidences[df_evidences.Model .== "flat", :lnZ][1]

	df_evidences
end

# ╔═╡ 42e909b4-92eb-4ed7-a19c-6e54b21ae07c
df_table = unstack(df_evidences, :Model, :Species, :ΔlnZ)

# ╔═╡ 83087c2e-852c-49c0-9195-c20e787b60e7
df_table_not_model = df_table[:, Not(:Model)];

# ╔═╡ 2d375a76-f06a-43bb-ba45-a944915288ff
idx_max_model = argmax(Matrix(df_table_not_model))

# ╔═╡ d9bcdac0-fed4-4b60-a5f7-e4035814ee92
idx_min_model = argmin(Matrix(df_table_not_model))

# ╔═╡ dd81cd2c-4a2e-47a0-9251-7e59d7bc3d45
df_table.Model[idx_max_model[1]], names(df_table_not_model)[idx_max_model[2]]

# ╔═╡ ef989991-08b1-470d-86b2-016c7c952e79
df_table.Model[idx_min_model[1]], names(df_table_not_model)[idx_min_model[2]]

# ╔═╡ f834b9fc-e410-442c-b085-8cccb8e30b71
latextabular(df_table, latex=false) |> PlutoUI.Text

# ╔═╡ 930ec094-7b11-48b8-818e-15c63ed6f8a5
dists = let
	data_py = cube["Na_TiO"]["clear"]["retr"]["samples"]
	PyDict{String, Vector}(data_py)
end

# ╔═╡ 54b5c81a-835a-461c-9dfd-2d938fac3bc4
Dict(k => median(v) ± std(v) for (k, v) ∈ dists)

# ╔═╡ 95f85651-ca70-4dd8-b82c-b531a966de90
cube |> keys

# ╔═╡ 41d8cd60-aacb-46c4-9344-29a8228f43cb
models = cube["K"] |> keys

# ╔═╡ d5197c68-0045-4ca7-b5c6-f9cd7852d896
models

# ╔═╡ 9c6a178c-ac7a-4ee5-ada2-05690158fcdb
models

# ╔═╡ 6ed6b481-4b84-485a-af2f-d2ecb943cd33
samples = pyconvert(Dict{String, Vector}, cube["Na"]["clear"]["retr"]["samples"])

# ╔═╡ 8ac2b5d2-eb4b-4a4f-99f8-52830e5ef9e8
keys(samples)

# ╔═╡ e43f1834-73dd-4859-b847-f4c552561897
begin
	##############
	# PLOT CONFIGS
	##############
	const FIG_TALL = (900, 1_200)
	const FIG_WIDE = (800, 600)
	const FIG_LARGE = (1_200, 1_000)
	const COLORS_SERIES = categorical_colors(:seaborn_colorblind, 8)
	const COLORS = let
		pal = Makie.ColorSchemes.Paired_8 |> reverse
		[pal[7:8] ; pal[5:6] ; pal[1:2]; :darkgrey]
	end
	
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

# ╔═╡ 65c3647d-d786-4a2c-89c8-f44aa6aa9f80
begin
	directions = [Vec2f(1, 1), Vec2f(1, 1), Vec2f(1, -1), Vec2f(1, -1), [Vec2f(1, 1), Vec2f(1, -1)], [Vec2f(1), Vec2f(1, -1)]]
	
	patterns = [[Makie.LinePattern(direction=hatch, tilesize=(16, 16), background_color=COLORS[i], linecolor=:white)
	    for (i, hatch) in enumerate(directions)]; :darkgrey]
end

# ╔═╡ 31473b50-9d12-41a2-8cbd-a12db0cb3706
bar_theme = (; color=patterns)

# ╔═╡ df43608e-7026-45ae-b87b-d7e0b6cea89c
let
	sort_order = sorter(model_names.keys)
	
	plt = data(df_evidences) *
		mapping(:Species => sorter(dashplus.(species)), :ΔlnZ => "ΔlnZ";
			dodge = :Model => sort_order,
			color = :Model => sort_order,
		) *
		visual(BarPlot; dodge_gap=0.0, gap=0.1)
	
	fg = draw(plt;
		legend = (
			position = :top,
			nbanks = 2,
			orientation = :horizontal,
			titleposition = :left,
			titlegap = 25,
		),
		axis = (; limits = (0.5, 8.5, 0, nothing)),
		palettes = bar_theme,
		figure = (; resolution=(1000, 500)),
	)

	ax = fg.figure[1, 1]
	vlines!(ax, [i + 0.5 for i in 1:7]; color=:lightgrey)
	
	# Label(fg.figure[2, 1], fname_suff;
	# 	halign = :right,
	# 	valign = :top,
	# 	tellwidth = false,
	# 	tellheight = false,
	# 	padding = (0, 30, 0, 20),
	# 	#textsize = 26,
	# )
	
	savefig(fg.figure, "$(FIG_DIR)/evidences_$(fname_suff).pdf")

	fg
end

# ╔═╡ e801501c-a882-4f2d-bbc1-40028c1c91d8
let
	fig = Figure(resolution=(1000, 500))
	ax = Axis(
		fig[1, 1],
		xlabel = "Wavelength (Å)",
		ylabel = "Transit depth (ppm)",
		#limits = (0.3, 1.1, 17_500, 21_000),
		limits = (4_000, 10_000, 17_000, 21_500),
		#limits = (4_600, 9_800, 17_000, 21_000)
		#limits = (4_000, 13_000, 15_500, 22_500),
		xticks = LinearTicks(8),
	)

	# for (i, sp) ∈ enumerate(species)
	# 	plot_retrieval!(ax, cube, sp, "clear"; color=COLORS[mod1(i, 6)], linewidth=1)
	# end
	# for (i, model) ∈ enumerate(models)
	# 	plot_retrieval!(ax, cube, "TiO", model; color=COLORS[i], linewidth=1)
	# end
	vlines!(ax, [5892.9, 7682.0, 8189.0], color=:grey, linestyle=:dash, linewidth=0.5)
	# plot_retrieval!(ax, cube, "Na_K_TiO", "spot"; color=COLORS[2])
	# plot_retrieval!(ax, cube, "Na_K_TiO", "cloud"; color=COLORS[3])
	# plot_retrieval!(ax, cube, "Na_K_TiO", "haze"; color=COLORS[5])
	# plot_retrieval!(ax, cube, "TiO", "clear"; color=COLORS[1], linewidth=3)
	plot_retrieval!(ax, cube, "Na", "clear"; color=COLORS[1])
	plot_retrieval!(ax, cube, "Na_K_TiO", "haze+spot"; show_scatter=false, color=COLORS[6])
	plot_retrieval!(ax, cube, "Na", "flat"; show_scatter=false, color=COLORS[7])
	
	# for sp ∈ species, model ∈ models
	# 	plot_retrieval!(ax, cube, sp, model)
	# end
	
	fpath_suff = basename(base_dir)
	if occursin("offs", fpath_suff)
		plot_instr!(ax, "retr_Clay_LDSS3.txt";
			color = :black,
			label = "LDSS3C",
		)
		plot_instr!(ax, "retr_Magellan_IMACS.txt";
			color = :white,
			label = "IMACS",
		)
	else
		plot_instr!(ax, "retr_IMACS+LDSS3C.txt";
			color = :white,
			label = "IMACS + LDSS3C",
		)
	end

	# Instrument

	axislegend(orientation=:horizontal)

	savefig(fig, "$(FIG_DIR)/retr_$(fname_suff).pdf")
	
	fig
end

# ╔═╡ e3708d6f-d9a9-4e42-b25b-2d1c333fddff
html"""
<style>
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
# ╟─0132b4ab-0447-4546-b412-ec598b20d21d
# ╠═a21aad0b-5998-420e-b437-d7ad262d0fe4
# ╠═60dc161c-2aa2-4264-884d-6da3ead0e57b
# ╟─d7ce97c1-82f2-46f1-a5ac-73e38e032fc8
# ╠═589afac8-0ea5-4962-b52b-7f035e91cf44
# ╟─2c12ec4d-1184-4755-8bd8-0d7cd59fa205
# ╠═65c3647d-d786-4a2c-89c8-f44aa6aa9f80
# ╠═9b43cef3-4583-4463-89ac-5f8678d1e843
# ╠═df43608e-7026-45ae-b87b-d7e0b6cea89c
# ╠═31473b50-9d12-41a2-8cbd-a12db0cb3706
# ╠═093156c7-9da7-4814-9260-5173f27fa497
# ╠═0f65d095-09af-44d2-907b-c30e2c16b609
# ╟─704fa634-eee0-4eef-aacf-f75f2b53f4d2
# ╠═a7c68d25-a799-421b-9799-38837fa8a188
# ╟─7b714c1e-2e3d-453f-a342-81df8283de5c
# ╠═65b51ff6-0991-491f-8945-dd889ffe71dd
# ╠═a0094689-a9d5-4810-baba-bd7a96c27839
# ╠═42e909b4-92eb-4ed7-a19c-6e54b21ae07c
# ╠═83087c2e-852c-49c0-9195-c20e787b60e7
# ╠═2d375a76-f06a-43bb-ba45-a944915288ff
# ╠═d9bcdac0-fed4-4b60-a5f7-e4035814ee92
# ╠═dd81cd2c-4a2e-47a0-9251-7e59d7bc3d45
# ╠═ef989991-08b1-470d-86b2-016c7c952e79
# ╠═f834b9fc-e410-442c-b085-8cccb8e30b71
# ╠═930ec094-7b11-48b8-818e-15c63ed6f8a5
# ╠═54b5c81a-835a-461c-9dfd-2d938fac3bc4
# ╟─1c4fe72d-9872-4969-a62a-5163b5009bbb
# ╠═d5197c68-0045-4ca7-b5c6-f9cd7852d896
# ╠═e801501c-a882-4f2d-bbc1-40028c1c91d8
# ╠═00a0f9c4-cd4d-4ae2-80b7-0c044239a571
# ╠═5569b57c-0585-4300-927b-5d089dde0f43
# ╠═db524678-9ee2-4934-b1bb-6a2f13bf0fa6
# ╠═cc011a66-37bd-4543-9a58-b11e1f785e52
# ╟─0f23e0d6-177d-4bf6-9660-f2c376b3146b
# ╠═9c6a178c-ac7a-4ee5-ada2-05690158fcdb
# ╠═95f85651-ca70-4dd8-b82c-b531a966de90
# ╠═41d8cd60-aacb-46c4-9344-29a8228f43cb
# ╠═6ed6b481-4b84-485a-af2f-d2ecb943cd33
# ╠═8ac2b5d2-eb4b-4a4f-99f8-52830e5ef9e8
# ╟─1eff1230-2423-4ac3-8e9b-f4e7bcd0121b
# ╟─eab74923-a084-468c-9b0d-c2cc98a23913
# ╠═44b3b8cd-4b83-4b27-a948-d1230489552f
# ╠═e43f1834-73dd-4859-b847-f4c552561897
# ╠═239a91a6-f68a-11eb-14fd-0ba8d08b08f9
# ╟─e3708d6f-d9a9-4e42-b25b-2d1c333fddff
