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
	rclone sync -P drive_ACCESS:papers/WASP-50b/$(DATA_DIR) $(DATA_DIR)
	```
	* [Direct link](https://app.box.com/s/hmqutmcs98rkip3aa4qyt68g3k1vk2du)
"""

# ╔═╡ 60dc161c-2aa2-4264-884d-6da3ead0e57b
@bind base_dir Select(glob("./data/retrievals/WASP50*"))

# ╔═╡ d7ce97c1-82f2-46f1-a5ac-73e38e032fc8
fit_R0 = "fitR0"

# ╔═╡ 8c6e3fd8-f6cb-4250-acb8-c17c00b1b2eb
fname_suff = let
	b = basename(base_dir)
	label = occursin("offs", b) ? "offs" : "no_offs"
	label *= "_$(fit_R0)"
end

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
function retrieval!(ax, model0, model_sampled; color=:blue, linewidth=3, label="")
	model = @rsubset model0 :wav < 1.0
	wav = model.wav .* 10_000
	wav_sampled = model_sampled.wav .* 10_000
	lines!(ax, wav, model.flux; color, linewidth, label)
	scatter!(ax, wav_sampled, model_sampled.flux;
		marker = :rect,
		markersize = 15,
		color,
		strokewidth = 1.5,
		strokecolor = color,
	)
	band!(ax, wav, model.flux_d, model.flux_u;
		color = (color, 0.25),
	)
end

# ╔═╡ 00a0f9c4-cd4d-4ae2-80b7-0c044239a571
function plot_retrieval!(ax, cube, sp, model; color=:blue, linewidth=3)
	retr_model, retr_model_sampled = get_retr_model(cube, sp, model)
	label = dashplus(sp) * " ($(model))"
	retrieval!(ax, retr_model, retr_model_sampled; color, linewidth, label)
end

# ╔═╡ 1eff1230-2423-4ac3-8e9b-f4e7bcd0121b
@mdx """
## Notebook setup
"""

# ╔═╡ eab74923-a084-468c-9b0d-c2cc98a23913
function savefig(fig, fpath)
	mkpath(dirname(fpath))
    save(fpath, fig)
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
end

# ╔═╡ 42e909b4-92eb-4ed7-a19c-6e54b21ae07c
df_table = unstack(df_evidences, :Model, :Species, :ΔlnZ)

# ╔═╡ 83087c2e-852c-49c0-9195-c20e787b60e7
df_table_not_model = df_table[:, Not(:Model)];

# ╔═╡ 2d375a76-f06a-43bb-ba45-a944915288ff
idx_max_model = argmax(Matrix(df_table_not_model))

# ╔═╡ dd81cd2c-4a2e-47a0-9251-7e59d7bc3d45
df_table.Model[idx_max_model[1]], names(df_table_not_model)[idx_max_model[2]]

# ╔═╡ f834b9fc-e410-442c-b085-8cccb8e30b71
latextabular(df_table, latex=false) |> PlutoUI.Text

# ╔═╡ 930ec094-7b11-48b8-818e-15c63ed6f8a5
dists = let
	data_py = cube["Na_TiO"]["clear"]["retr"]["samples"]
	PyDict{String, Vector}(data_py)
end

# ╔═╡ 54b5c81a-835a-461c-9dfd-2d938fac3bc4
Dict(k => median(v) ± std(v) for (k, v) ∈ dists)

# ╔═╡ e43f1834-73dd-4859-b847-f4c552561897
begin
	##############
	# PLOT CONFIGS
	##############
	const FIG_TALL = (900, 1_200)
	const FIG_WIDE = (800, 600)
	const FIG_LARGE = (1_200, 1_000)
	const COLORS_SERIES = to_colormap(:seaborn_colorblind, 9)
	const COLORS = let
		pal = Makie.ColorSchemes.Paired_8 |> reverse
		[pal[7:8] ; pal[5:6] ; pal[1:2]]
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

# ╔═╡ df43608e-7026-45ae-b87b-d7e0b6cea89c
let
	sort_order = sorter(model_names.keys)
	
	plt = data(df_evidences) *
		mapping(:Species => sorter(dashplus.(species)), :ΔlnZ => "ΔlnZ";
			dodge = :Model => sort_order,
			color = :Model => sort_order,
		) *
		visual(BarPlot)
	
	fg = draw(plt;
		axis = (; limits=(nothing, nothing, 0, 3.5)),
		legend = (
			position = :top,
			nbanks = 2,
			orientation = :horizontal,
			titleposition = :left,
			titlegap = 25,
		),
		palettes = (; color=COLORS),
		figure = (; resolution=FIG_WIDE),
	)
	
	Label(fg.figure[2, 1], fname_suff;
		halign = :right,
		valign = :top,
		tellwidth = false,
		tellheight = false,
		padding = (0, 30, 0, 20),
		#textsize = 26,
	)
	
	savefig(fg.figure, "$(FIG_DIR)/evidences_$(fname_suff).png")

	fg
end

# ╔═╡ e801501c-a882-4f2d-bbc1-40028c1c91d8
let
	fig = Figure(resolution=(1_200, 400))
	ax = Axis(
		fig[1, 1],
		xlabel = "Wavelength (Å)",
		ylabel = "Transit depth (ppm)",
		#limits = (0.3, 1.1, 17_500, 21_000),
		#limits = (4_600, 9_800, 15_500, 22_500),
		limits = (4_600, 9_800, 17_500, 21_500)
	)

	plot_retrieval!(ax, cube, "Na_TiO", "spot", color=COLORS[2])
	plot_retrieval!(ax, cube, "Na", "cloud", color=COLORS[3])
	plot_retrieval!(ax, cube, "TiO", "haze", color=COLORS[5])
	plot_retrieval!(ax, cube, "Na_TiO", "clear", color=COLORS[1], linewidth=6)
	
	fpath_suff = basename(base_dir)
	if occursin("offs", fpath_suff)
		plot_instr!(ax, "retr_Clay_LDSS3.txt";
			color = :black,
			label = "LDSS3",
		)
		plot_instr!(ax, "retr_Magellan_IMACS.txt";
			color = :white,
			label = "IMACS",
		)
	else
		plot_instr!(ax, "retr_Magellan_IMACS.txt";
			color = :white,
			label = "IMACS + LDSS3",
		)
	end

	# Instrument

	axislegend(orientation=:horizontal)

	savefig(fig, "$(FIG_DIR)/retr_$(fname_suff).png")
	
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
# ╟─60dc161c-2aa2-4264-884d-6da3ead0e57b
# ╠═d7ce97c1-82f2-46f1-a5ac-73e38e032fc8
# ╠═8c6e3fd8-f6cb-4250-acb8-c17c00b1b2eb
# ╟─093156c7-9da7-4814-9260-5173f27fa497
# ╠═0f65d095-09af-44d2-907b-c30e2c16b609
# ╟─704fa634-eee0-4eef-aacf-f75f2b53f4d2
# ╠═a7c68d25-a799-421b-9799-38837fa8a188
# ╟─7b714c1e-2e3d-453f-a342-81df8283de5c
# ╠═65b51ff6-0991-491f-8945-dd889ffe71dd
# ╠═a0094689-a9d5-4810-baba-bd7a96c27839
# ╠═df43608e-7026-45ae-b87b-d7e0b6cea89c
# ╠═42e909b4-92eb-4ed7-a19c-6e54b21ae07c
# ╠═83087c2e-852c-49c0-9195-c20e787b60e7
# ╠═2d375a76-f06a-43bb-ba45-a944915288ff
# ╠═dd81cd2c-4a2e-47a0-9251-7e59d7bc3d45
# ╠═f834b9fc-e410-442c-b085-8cccb8e30b71
# ╠═930ec094-7b11-48b8-818e-15c63ed6f8a5
# ╠═54b5c81a-835a-461c-9dfd-2d938fac3bc4
# ╟─1c4fe72d-9872-4969-a62a-5163b5009bbb
# ╠═e801501c-a882-4f2d-bbc1-40028c1c91d8
# ╠═00a0f9c4-cd4d-4ae2-80b7-0c044239a571
# ╠═5569b57c-0585-4300-927b-5d089dde0f43
# ╠═db524678-9ee2-4934-b1bb-6a2f13bf0fa6
# ╠═cc011a66-37bd-4543-9a58-b11e1f785e52
# ╟─1eff1230-2423-4ac3-8e9b-f4e7bcd0121b
# ╟─eab74923-a084-468c-9b0d-c2cc98a23913
# ╠═44b3b8cd-4b83-4b27-a948-d1230489552f
# ╠═e43f1834-73dd-4859-b847-f4c552561897
# ╠═239a91a6-f68a-11eb-14fd-0ba8d08b08f9
# ╟─e3708d6f-d9a9-4e42-b25b-2d1c333fddff
