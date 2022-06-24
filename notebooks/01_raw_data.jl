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

# ╔═╡ 3433ed02-c27c-4fe5-bfda-a5108a58407c
begin
	import Pkg
	Pkg.activate(Base.current_project())

	using PlutoUI
	import MarkdownLiteral: @mdx
	using AlgebraOfGraphics, CairoMakie
	using CSV, DataFrames, DataFramesMeta, DelimitedFiles, Glob, OrderedCollections
	using CCDReduction, ImageFiltering, Statistics
	using Latexify, Printf
end

# ╔═╡ e0cd40b4-c616-41b1-8d4d-b9283a46863d
using BenchmarkTools

# ╔═╡ d8019fa7-380d-4f40-9e08-420a32c34483
begin
	const DATA_DIR = "data/raw"
	const FIG_DIR = "figures/raw"
	TableOfContents()
end

# ╔═╡ fb39c593-86bd-4d4c-b9ec-e5e212a4de98
@mdx """
# Raw data

In this notebook we will view bias-subtracted sample frames from each instrument to verify the mask slit placement, as well as the positions of WASP-50 and its comparison stars on each CCD chip.

For each fits files (1 per chip), we extract the region of each chip dedicated to main data collection (`DATASEC`), and then subtract out the median bias level measured in the overscan region (`BIASSEC`). We then store the subtracted data in a data cube for each instrument.

!!! tip "Data download"
	```
	rclone sync -P ACCESS_box:WASP-50b/$(DATA_DIR) $(DATA_DIR)
	```
	* [Direct link](https://app.box.com/s/yfa96kreb67mjmvqretpsdr5s5u5xskr)

	Outline:
	```
	raw/
	├── [4.1k]  IMACS/
	│   ├── [4.1k]  ut131219/
	│   │   ├── [4.6M]  ift0026c1.fits
	│   │   ├── [4.6M]  ⋅ ⋅ ⋅
	│   │   ├── [4.6M]  ift0026c8.fits
	│   │   └── [ 510]  WASP50.coords
	│   ├── [4.1k]  ut150927/
	│   │   ├── [4.6M]  ift0543c1.fits
	│   │   ├── [4.6M]  ⋅ ⋅ ⋅
	│   │   ├── [4.6M]  ift0543c8.fits
	│   │   └── [ 413]  WASP50.coords
	│   └── [4.1k]  ut161211/
	│       ├── [4.6M]  ift0046c1.fits
	│       ├── [4.6M]  ⋅ ⋅ ⋅
	│       ├── [4.6M]  ift0046c8.fits
	│       └── [ 391]  WASP50.coords
	├── [4.1k]  LDSS3C/
	│   └── [4.1k]  ut150927/
	│       ├── [2.8M]  ccd0606c1.fits
	│       ├── [2.8M]  ccd0606c2.fits
	│       └── [1.6k]  LDSS3.150927e.txt
	├── [3.2k]  wasp50s.SMF*
	└── [4.1k]  wbins/
	    ├── [ 326]  w50_bins.dat
	    ├── [ 342]  w50_bins_LDSS3.dat
	    ├── [ 247]  w50_bins_species.dat
	    └── [ 390]  w50_bins_ut131219.dat
	```
"""

# ╔═╡ 90845d70-35d9-402d-8936-74936b069577
@mdx """
## IMACS 🎱
Starting with IMACS, let's first select the night we would like to visualize from the dropdown menu below:
"""

# ╔═╡ 8b9581db-71c6-42b6-915b-bde307755bcd
@bind DATA_DIR_IMACS Select(
	glob("$(DATA_DIR)/IMACS/ut*"),
	default = "$(DATA_DIR)/IMACS/ut161211"
)

# ╔═╡ c040905d-23bf-4484-8743-98d917db9c81
@mdx """
This figure uses the following data:
"""

# ╔═╡ fb6e6221-8136-44e2-979b-ecbbd71f740d
df_sci_IMACS = fitscollection(DATA_DIR_IMACS, abspath=false)

# ╔═╡ 0d42f6f9-d789-46a3-9e9a-381dbed2d5a5
@mdx """
We compute the bias-subtracted science frame for each of these $(nrow(df_sci_IMACS)) chips:
"""

# ╔═╡ c83c2161-95e2-4d08-9934-6d9c12c42a44
process_frames(df) = @views map(ccds(df)) do img
	img_data = img[CCDReduction.fits_indices(img.hdr["DATASEC"])...]
	img_bias = img[CCDReduction.fits_indices(img.hdr["BIASSEC"])...]
	corr = img_data .- median(img_bias)
end

# ╔═╡ 86b4ace1-a891-41e5-a29f-f7eee5f8fb17
frames_IMACS = process_frames(df_sci_IMACS)

# ╔═╡ 27b793f3-7a7c-48ad-8302-deffa2dd017b
@mdx """
And finally, we overlay the coordinates of the science mask onto the CCD array:
"""

# ╔═╡ 26feb668-4e7e-4a9d-a2b6-a5dac81e3ab7
coords_IMACS = CSV.read("$(DATA_DIR_IMACS)/WASP50.coords", DataFrame;
	delim = ' ',
	ignorerepeated = true,
	header = ["target", "chip", "x", "y"],
);

# ╔═╡ 8c5a4e21-897f-4fbc-bd4f-18adf71fa926
transits = merge(
	Dict("$(dirname(DATA_DIR_IMACS))/ut$(d)" => "Transit $(i) (IMACS)"
		for (i, d) ∈ enumerate(("131219", "150927", "161211"))
	),
)

# ╔═╡ 06a834f0-8c90-4013-af34-725166970969
@mdx """
## LDSS3C 2️⃣

We follow the same operations to visualize the two chips for LDSS3C below.
"""

# ╔═╡ 0e66d467-1098-46dc-8d06-36d488b14637
@bind DATA_DIR_LDSS3 Select(glob("$(DATA_DIR)/LDSS3C/ut*"))

# ╔═╡ 97445546-41c0-46ca-9339-a0f9fb87b2c1
f1(img) = heatmap(img)

# ╔═╡ 0b469e00-386f-401b-949d-e0715b4d9def
f2(img) = image(img)

# ╔═╡ c6205ad9-4d3c-420a-a279-81731d83603b
@mdx """
This figure uses the following data:
"""

# ╔═╡ 5c6e6f7b-70e0-49a8-b064-60dcf1440223
df_sci_LDSS3 = fitscollection(DATA_DIR_LDSS3, abspath=false)

# ╔═╡ c488270a-3126-4e38-a0c8-ee242115a3ea
frames_LDSS3 = process_frames(df_sci_LDSS3)

# ╔═╡ 83a9357d-836b-4cee-a41f-eabc8f3f12e7
coords_LDSS3 = DataFrame((
	(target="c21",     chip="c1", x=25,  y=1300),
	(target="c06",     chip="c1", x=150,  y=1200),
	(target="WASP-50", chip="c1", x=322, y=800),
	(target="c15",     chip="c2", x=198,  y=1000)
));

# ╔═╡ bf8ef5a9-0806-44b4-907d-c95d6926dabb
function plot_frame!(ax;
	img=frames_LDSS3[1], ch=1, coords=coords_LDSS3, stepsize=32, hm_kwargs=(),
)
	hm = plot!(ax, img;
		colormap = :magma,
		hm_kwargs...,
	)

	# Add labels
	chip = "c$ch"
	for obj in eachrow(@subset coords :chip .== chip)
		coord = (obj.x, obj.y) ./ stepsize
		text!(ax, obj.target;
			position = coord,
			textsize = 11,
			align = (:center, :center),
			color = :black,
			rotation = π/2,
		)
	end

	return hm
end

# ╔═╡ 8fadd0b6-6ff8-42e5-9014-4e79593e3502
@mdx """
## Wavelength bins 🌈

We show the wavelength bins used for each instrument here:
"""

# ╔═╡ 9e3d2013-a4ac-4413-912e-aa94046e2f44
@bind fpath_wbin Select(glob("$(DATA_DIR)/wbins/*.dat"))

# ╔═╡ f2ccf230-f2ac-43c2-b313-8821ef69a1e7
df_wbins = let
	df = CSV.read(fpath_wbin, DataFrame;
		header = [:wav_d, :wav_u],
		comment = "#",
		delim = ' ',
		ignorerepeated = true,
	)
	@select df begin
		:wav_cen = mean((:wav_d, :wav_u))
		:wav_Δ = :wav_u .- :wav_d
	end
end

# ╔═╡ 0d2476b1-2864-4bfc-ac37-f771aab77368
latextabular(df_wbins, latex=false) |> PlutoUI.Text

# ╔═╡ 4480ae72-3bb2-4e17-99be-28afc756332a
@mdx """
## Notebook setup 🔧
"""

# ╔═╡ b44b6591-d57b-40bd-810c-a41386412b6c
function savefig(fig, fpath)
	mkpath(dirname(fpath))
    save(fpath, fig; pt_per_unit=1)
	@info "Saved to: $(fpath)"
end

# ╔═╡ db4a4cd8-c5e8-4124-935f-0666f6e73fe2
begin
	##############
	# PLOT CONFIGS
	##############
	const FIG_TALL = 72 .* (6, 8)
	const FIG_WIDE = 72 .* (12, 6)
	const FIG_LARGE = 72 .* (12, 12)
	const COLORS_SERIES = categorical_colors(:seaborn_colorblind, 9)
	const COLORS = parse.(Makie.Colors.Colorant,
		[
			"#66C2A5",  # Green
			"#FDBF6F",  # Yellow
			"#FF7F00",  # Orange
			"#1F78B4",  # Blue
		]
	)

	set_aog_theme!()
	update_theme!(
		Theme(
			Axis = (
				xlabelsize = 18,
				ylabelsize = 18,
				topspinevisible = true,
				rightspinevisible = true,
				topspinecolor = :darkgrey,
				rightspinecolor = :darkgrey,
			),
			Label = (
				textsize = 18,
				font = AlgebraOfGraphics.firasans("Medium"),
			),
			Lines = (linewidth=3,),
			Scatter = (linewidth=10, strokewidth=0),
			Text = (font = AlgebraOfGraphics.firasans("Regular"), textsize=18),
			palette = (color=COLORS, patchcolor=[(c, 0.35) for c in COLORS]),
			figure_padding = (0, 1.5, 0, 0),
			fontsize = 18,
			rowgap = 5,
			colgap = 5,
		)
	)

	COLORS
end

# ╔═╡ 3a6ab0c0-ba08-4151-9646-c19d45749b9f
let
	fig = Figure(resolution = FIG_LARGE)
	hm = nothing
	stepsize = 1
	grid = CartesianIndices((2, 4))
	chip_order = [1, 6, 2, 5, 3, 8, 4, 7]
	for (g, ch) ∈ zip(grid, chip_order)
		# Set chip orientation based on IMACS conventions and apply chip labels
		i, j = g.I
		if ch ∈ 1:4
			xflip, yflip = x -> x, x -> x
			Label(fig[1, j], "c$(ch)", tellwidth=false)
		else
			xflip, yflip = reverse, reverse
			Label(fig[4, j], "c$(ch)", tellwidth=false)
		end
		ax = Axis(fig[i+1, j], xreversed=false, yreversed=true)
		img = @view(
			frames_IMACS[ch][xflip(begin:stepsize:end), yflip(begin:stepsize:end)]
		)'
		hm = plot_frame!(ax;
			img,
			ch,
			coords = coords_IMACS,
			stepsize,
			hm_kwargs = (highclip=:darkgrey, colorrange=(0, 200)),
		)
    end

	rowgap!(fig.layout, 0)
	colgap!(fig.layout, 0)
	
	Colorbar(fig[2:3, end+1], hm, width=20, label="Counts", ticksvisible=false)
	axs = filter(x -> x isa Axis, fig.content)
	linkaxes!(axs...)
	hidedecorations!.(axs)

	Label(fig[0, end], transits[DATA_DIR_IMACS], tellwidth=false, halign=:right)
	
    savefig(fig, "$(FIG_DIR)/sci_IMACS_$(basename(DATA_DIR_IMACS)).png")

	fig
end

# ╔═╡ 71ba9181-90e4-4d12-97c0-462b3f1df077
let
	fig = Figure(resolution = FIG_LARGE)
	stepsize = 1
	hm = nothing
	for j ∈ 1:2
		ax = Axis(fig[2, j])
		img = @view(frames_LDSS3[j][begin:stepsize:end, begin:stepsize:end])'
		hm = plot_frame!(ax;
			img,
			ch = j,
			coords = coords_LDSS3,
			stepsize,
			hm_kwargs = (highclip=:darkgrey, colorrange=(0, 200)),
		)
		Label(fig[1, j], "c$(j)", tellwidth=false)
		Label(fig[3, j], "c$(j)", tellwidth=false)
	end

	rowgap!(fig.layout, 0)
	colgap!(fig.layout, 0)

	Colorbar(fig[2, end+1], hm, width=20, ticksvisible=false, label="Counts")
	
	axs = filter(x -> x isa Axis, fig.content)
	linkaxes!(axs...)
	hidedecorations!.(axs)

	Label(fig[0, end], "Transit 2 (LDSS3C)", tellwidth=false, halign=:right)
	
    savefig(fig, "$(FIG_DIR)/sci_LDSS3_$(basename(DATA_DIR_LDSS3)).png")
	
	fig
end

# ╔═╡ Cell order:
# ╟─fb39c593-86bd-4d4c-b9ec-e5e212a4de98
# ╠═d8019fa7-380d-4f40-9e08-420a32c34483
# ╟─90845d70-35d9-402d-8936-74936b069577
# ╟─8b9581db-71c6-42b6-915b-bde307755bcd
# ╠═3a6ab0c0-ba08-4151-9646-c19d45749b9f
# ╟─c040905d-23bf-4484-8743-98d917db9c81
# ╠═fb6e6221-8136-44e2-979b-ecbbd71f740d
# ╟─0d42f6f9-d789-46a3-9e9a-381dbed2d5a5
# ╠═c83c2161-95e2-4d08-9934-6d9c12c42a44
# ╠═86b4ace1-a891-41e5-a29f-f7eee5f8fb17
# ╟─27b793f3-7a7c-48ad-8302-deffa2dd017b
# ╠═26feb668-4e7e-4a9d-a2b6-a5dac81e3ab7
# ╟─8c5a4e21-897f-4fbc-bd4f-18adf71fa926
# ╠═bf8ef5a9-0806-44b4-907d-c95d6926dabb
# ╟─06a834f0-8c90-4013-af34-725166970969
# ╟─0e66d467-1098-46dc-8d06-36d488b14637
# ╠═71ba9181-90e4-4d12-97c0-462b3f1df077
# ╠═97445546-41c0-46ca-9339-a0f9fb87b2c1
# ╠═0b469e00-386f-401b-949d-e0715b4d9def
# ╠═e0cd40b4-c616-41b1-8d4d-b9283a46863d
# ╟─c6205ad9-4d3c-420a-a279-81731d83603b
# ╠═5c6e6f7b-70e0-49a8-b064-60dcf1440223
# ╠═c488270a-3126-4e38-a0c8-ee242115a3ea
# ╠═83a9357d-836b-4cee-a41f-eabc8f3f12e7
# ╟─8fadd0b6-6ff8-42e5-9014-4e79593e3502
# ╟─9e3d2013-a4ac-4413-912e-aa94046e2f44
# ╠═f2ccf230-f2ac-43c2-b313-8821ef69a1e7
# ╠═0d2476b1-2864-4bfc-ac37-f771aab77368
# ╟─4480ae72-3bb2-4e17-99be-28afc756332a
# ╟─b44b6591-d57b-40bd-810c-a41386412b6c
# ╠═db4a4cd8-c5e8-4124-935f-0666f6e73fe2
# ╠═3433ed02-c27c-4fe5-bfda-a5108a58407c
