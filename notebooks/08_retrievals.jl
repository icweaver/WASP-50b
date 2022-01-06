### A Pluto.jl notebook ###
# v0.17.5

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
	using AlgebraOfGraphics
	using CairoMakie
	using CSV, DataFrames, DataFramesMeta
	using DelimitedFiles
	using Glob
	using ImageFiltering
	using Latexify
	using Measurements, Statistics
	using OrderedCollections
end

# ╔═╡ f2b2db32-f7fb-4735-849d-5bee761a5e85
begin
	ENV["PYTHON"] = ""
	Pkg.build("PyCall")
	using PyCall
	const Conda = PyCall.Conda
	Conda.add("lightkurve", :WASP50b)
end

# ╔═╡ 0132b4ab-0447-4546-b412-ec598b20d21d
md"""
# Retrievals

!!! note "Data download"

	```
	rclone sync -P drive_ACCESS:papers/WASP-50b/data/retrievals data/retrievals
	
	rsync -azRP $H:/pool/sao_access/iweaver/exoretrievals/./"WASP50_offs/*/{*.pkl,*.txt}" data/retrievals/ --exclude "*nest*"

	rsync -azRP $H:/pool/sao_access/iweaver/exoretrievals/WASP50_offs/"*/./WASP50_offs" data/retrievals/ --exclude "*nest*" 
	```

$(TableOfContents(title="📖 Table of Contents"))
"""

# ╔═╡ 60dc161c-2aa2-4264-884d-6da3ead0e57b
@bind base_dir Select(glob("./data/retrievals/WASP50*"))

# ╔═╡ 56971ef4-7512-4e85-ac41-ee446006457f
const FIG_PATH = "figures/retrievals"

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

# ╔═╡ d4356ef7-abd7-47dd-83e3-38b6a782509e
#df_wide = unstack(df_evidences, :Model, :Species, :lnZ, renamecols = dashplus)

# ╔═╡ a0094689-a9d5-4810-baba-bd7a96c27839
dashplus(x) = replace(x, '_' => '+')

# ╔═╡ 1c4fe72d-9872-4969-a62a-5163b5009bbb
md"""
## Retreived transmission spectra
"""

# ╔═╡ 5569b57c-0585-4300-927b-5d089dde0f43
function plot_instr!(ax, fpath; color=:blue, label="add label")
	retr_instr = CSV.read(
		"$(base_dir)/$(fpath)", DataFrame;
		header = [:wav, :flux, :wav_d, :wav_u, :flux_err],
		comment = "#",
	)
	errorbars!(ax, retr_instr.wav, retr_instr.flux, retr_instr.flux_err)
	scatter!(ax, retr_instr.wav, retr_instr.flux;
		markersize = 15,
		color = color,
		strokewidth=1.5,
		label = label
	)
end

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
	label = dashplus(sp) * " ($(model))"
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

	@transform! df_evidences begin
		:lnZ = :lnZ .- minimum(:lnZ)
		:Species = dashplus.(:Species)
		:ΔlnZ_m = (:lnZ .± :lnZ_err) .- minimum(:lnZ .± :lnZ_err)
	end
end

# ╔═╡ 42e909b4-92eb-4ed7-a19c-6e54b21ae07c
unstack(df_evidences, :Model, :Species, :lnZ)

# ╔═╡ 1eff1230-2423-4ac3-8e9b-f4e7bcd0121b
md"""
## Notebook setup
"""

# ╔═╡ eab74923-a084-468c-9b0d-c2cc98a23913
function savefig(fig, fpath)
	mkpath(dirname(fpath))
    save(fpath, fig)
	@info "Saved to: $(fpath)"
end

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
		mapping(:Species => sorter(dashplus.(species)), :lnZ => "ΔlnZ";
			dodge = :Model => sort_order,
			color = :Model => sort_order => "",
		) *
		visual(BarPlot)
	
	fg = draw(plt;
		axis = (; limits=(nothing, nothing, 0, 3.5)),
		legend = (
			position = :top,
			nbanks = 2,
			orientation = :horizontal,
			titlevisible = false,
		),
		palettes = (; color=COLORS),
	)

	savefig(fg.figure, "$(FIG_PATH)/evidences_$(basename(base_dir)).png")
end

# ╔═╡ e801501c-a882-4f2d-bbc1-40028c1c91d8
let
	fig = Figure(resolution=(800, 500))
	ax = Axis(
		fig[1, 1],
		xlabel = "Wavelength (μm)",
		ylabel = "Transit depth (ppm)",
		#limits = (0.3, 1.1, 17_500, 21_000),
	)

	plot_retrieval!(ax, cube, "Na_TiO", "clear", color=COLORS[1])
	plot_retrieval!(ax, cube, "TiO", "clear+haze+spot", color=COLORS[6])
	plot_instr!(ax, "retr_Magellan_IMACS.txt", color=:white, label="IMACS + LDSS3")
	#plot_instr!(ax, "retr_Clay_LDSS3.txt", color=:lightblue, label="IMACS + LDSS3")

	# Instrument

	axislegend(orientation=:horizontal)

	savefig(fig, "$(FIG_PATH)/retr_$(basename(base_dir)).png")
	
	fig
end

# ╔═╡ 01bbee9f-66bf-4e08-a91b-9870def4e62a
@with_terminal Conda.list(:WASP50b)

# ╔═╡ Cell order:
# ╟─0132b4ab-0447-4546-b412-ec598b20d21d
# ╠═60dc161c-2aa2-4264-884d-6da3ead0e57b
# ╟─56971ef4-7512-4e85-ac41-ee446006457f
# ╠═d7ce97c1-82f2-46f1-a5ac-73e38e032fc8
# ╠═093156c7-9da7-4814-9260-5173f27fa497
# ╠═0f65d095-09af-44d2-907b-c30e2c16b609
# ╟─704fa634-eee0-4eef-aacf-f75f2b53f4d2
# ╠═a7c68d25-a799-421b-9799-38837fa8a188
# ╠═7b714c1e-2e3d-453f-a342-81df8283de5c
# ╠═65b51ff6-0991-491f-8945-dd889ffe71dd
# ╟─d4356ef7-abd7-47dd-83e3-38b6a782509e
# ╠═a0094689-a9d5-4810-baba-bd7a96c27839
# ╠═df43608e-7026-45ae-b87b-d7e0b6cea89c
# ╠═42e909b4-92eb-4ed7-a19c-6e54b21ae07c
# ╟─1c4fe72d-9872-4969-a62a-5163b5009bbb
# ╠═e801501c-a882-4f2d-bbc1-40028c1c91d8
# ╠═5569b57c-0585-4300-927b-5d089dde0f43
# ╠═00a0f9c4-cd4d-4ae2-80b7-0c044239a571
# ╠═db524678-9ee2-4934-b1bb-6a2f13bf0fa6
# ╠═cc011a66-37bd-4543-9a58-b11e1f785e52
# ╟─41a233c7-5357-453c-b7ad-36fdf9f709cb
# ╠═44b3b8cd-4b83-4b27-a948-d1230489552f
# ╟─1eff1230-2423-4ac3-8e9b-f4e7bcd0121b
# ╟─eab74923-a084-468c-9b0d-c2cc98a23913
# ╠═e43f1834-73dd-4859-b847-f4c552561897
# ╠═01bbee9f-66bf-4e08-a91b-9870def4e62a
# ╠═f2b2db32-f7fb-4735-849d-5bee761a5e85
# ╠═239a91a6-f68a-11eb-14fd-0ba8d08b08f9
