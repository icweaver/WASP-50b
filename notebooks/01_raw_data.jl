### A Pluto.jl notebook ###
# v0.17.3

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
	
	using AlgebraOfGraphics, CairoMakie
	using CSV, DataFrames, DataFramesMeta
	using CCDReduction
	using PlutoUI
	using Glob
	using Statistics
	
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
			Text = (; font=AlgebraOfGraphics.firasans("Medium")),
			Lines = (linewidth=3, cycle=Cycle([:color, :linestyle], covary=true)),
			Scatter = (linewidth=10,),
			fontsize = 18,
			rowgap = 0,
			colgap = 0,
		)
	)
	
	COLORS = Makie.wong_colors()
end

# ╔═╡ fb39c593-86bd-4d4c-b9ec-e5e212a4de98
md"""
# Raw Data

In this notebook we will view a bias-subtracted sample frame from each instrument. We do this to verify the mask slit placement, as well as the positions of WASP-50 and its comparison stars on each CCD chip.

For each fits files (1 for each chip), we extract just the portion on each chip specified by `DATASEC`, and then subtract out the median bias level measured in the overscan region (defined by `BIASSEC`). We then store the subtracted data in a data cube for each instrumet

$(TableOfContents())

!!! note "Data download"
	```
	 rclone sync ACCESS_box:WASP-50b/data/raw data/raw
	```
"""

# ╔═╡ 90845d70-35d9-402d-8936-74936b069577
md"""
## IMACS
Starting with IMACS, let's first select the night we would like to visualize from the dropdown menu below:
"""

# ╔═╡ 8b9581db-71c6-42b6-915b-bde307755bcd
@bind DATA_DIR_IMACS Select(glob("data/raw/IMACS/ut*"))

# ╔═╡ fb6e6221-8136-44e2-979b-ecbbd71f740d
df_sci_IMACS = fitscollection(DATA_DIR_IMACS)

# ╔═╡ b8de138b-282a-44b3-be8c-7fea2c88030d
img_overscan_corrected = @views map(ccds(df_sci_IMACS)) do img
	img_data = img[CCDReduction.fits_indices(img.hdr["DATASEC"])...]
	img_bias = img[CCDReduction.fits_indices(img.hdr["BIASSEC"])...]
	corr = img_data .- median(img_bias)
end

# ╔═╡ c83c2161-95e2-4d08-9934-6d9c12c42a44
process_frames(df) = @views map(ccds(df)) do img
	img_data = img[CCDReduction.fits_indices(img.hdr["DATASEC"])...]
	img_bias = img[CCDReduction.fits_indices(img.hdr["BIASSEC"])...]
	corr = img_data .- median(img_bias)
end

# ╔═╡ 0d42f6f9-d789-46a3-9e9a-381dbed2d5a5
md"""
We next compute the bias subtracted science frame for each of these $(nrow(df_sci_IMACS)) chips:
"""

# ╔═╡ 86b4ace1-a891-41e5-a29f-f7eee5f8fb17
frames_IMACS = process_frames(df_sci_IMACS)

# ╔═╡ fed5dd71-98e8-4212-b574-6150bccce2f7
plot(frames_IMACS[3], colormap = :magma,
			colorrange = (0, 500))

# ╔═╡ 6732965f-9ab6-4b84-9703-c212bb89d353
chip_order = [1, 6, 2, 5, 3, 8, 4, 7]

# ╔═╡ 26feb668-4e7e-4a9d-a2b6-a5dac81e3ab7
coords_IMACS = CSV.read("$(DATA_DIR_IMACS)/WASP50.coords", DataFrame;
	delim = ' ',
	ignorerepeated = true,
	header = ["target", "chip", "x", "y"],
)

# ╔═╡ be8399c7-d817-476a-bf26-becaea8b6522
@with_terminal begin
	grid = CartesianIndices((2, 4))
	chip_order = [1, 6, 2, 5, 3, 8, 4, 7]

	for (ch, g) ∈ zip(chip_order, grid)
		println(ch, g)
	end
end

# ╔═╡ d00a2334-bee8-4819-b11b-17a5c7eb84e0
2 ∈ 4:8

# ╔═╡ 3a6ab0c0-ba08-4151-9646-c19d45749b9f
let
	fig = Figure(resolution = (800, 600))
	step = 1 # For quick testing
	hms = []
	axs = []
	grid = CartesianIndices((2, 4))
	chip_order = [1, 6, 2, 5, 3, 8, 4, 7]
	for (g, ch) ∈ zip(grid, chip_order)
		# Set chip orientation based on IMACS conventions
		if ch ∈ 5:8
			xflip = reverse
			yflip = reverse
			xreversed = false
			yreversed = true
		else
			xflip = x -> x
			yflip = x -> x
			xreversed = false
			yreversed = true
		end
		ax = Axis(fig[g.I...], xreversed=xreversed, yreversed=yreversed)
		d = frames_IMACS[ch][yflip(begin:step:end), xflip(begin:step:end)]'
		nrows, ncols = size(d)
		hm = plot!(
			ax,
			d,
			colormap = :magma,
			colorrange = (0, 500),
		)

		# Label objects on chip
		chip = "c$ch"
		for obj in eachrow(@subset coords_IMACS :chip .== chip)
			coord = (obj.x, obj.y) ./ step
			text!(ax, split(obj.target, "_")[1];
				position = coord,
				textsize = 11,
				align = (:center, :center),
				color = :white,
				rotation = π/2,
			)
		end
		
		# Label chip
		text!(
			"$(chip)",
			position = (0, ncols),
			color = :yellow,
			align = (:left, :baseline),
			offset = 0.5 .* (20, 20),
		)
		push!(axs, ax)
		push!(hms, hm)
	end
	axs = reshape(axs, 2, 4)
	
	Colorbar(fig[:, end+1], hms[1], width=20, label="Counts",)
	
	linkaxes!(axs...)
	hidedecorations!.(axs)
	
	path = "../../ACCESS_WASP-50b/figures/frames"
	mkpath(path)
	save("$(path)/sci_imacs_$(basename(DATA_DIR_IMACS)).png", fig)
	
	fig #|> as_svg
end

# ╔═╡ 06a834f0-8c90-4013-af34-725166970969
md"""
## LDSS3
"""

# ╔═╡ 0e66d467-1098-46dc-8d06-36d488b14637
@bind DATA_DIR_LDSS3 Select(glob("data/raw/LDSS3/ut*"))

# ╔═╡ 5c6e6f7b-70e0-49a8-b064-60dcf1440223
df_sci_LDSS3 = fitscollection(DATA_DIR_LDSS3)

# ╔═╡ c488270a-3126-4e38-a0c8-ee242115a3ea
frames_LDSS3 = process_frames(df_sci_LDSS3)

# ╔═╡ 83a9357d-836b-4cee-a41f-eabc8f3f12e7
coords_LDSS3 = DataFrame((
	(target="c21",     chip="c1", x=25,  y=1300),
	(target="c06",     chip="c1", x=150,  y=1200),
	(target="WASP-50", chip="c1", x=322, y=800),
	(target="c15",     chip="c2", x=198,  y=1000)
))

# ╔═╡ 71ba9181-90e4-4d12-97c0-462b3f1df077
let
	fig = Figure(resolution = (800, 600))
	step = 1
	chip_idx = 1
	hms = []
	axs = []
	for j in 1:2, i in 1
		ax = Axis(fig[i, j])
		d = frames_LDSS3[chip_idx][begin:step:end, begin:step:end]'
		nrows, ncols = size(d)
		hm = plot!(
			ax,
			d,
			colormap = :magma,
			colorrange = (0, 2500),
		)

		chip = "c$j"

		# Label objects
		for obj in eachrow(@subset coords_LDSS3 :chip .== chip)
			coord = (obj.x, obj.y)
			text!(ax, obj.target;
				position = coord,
				textsize = 11,
				align = (:center, :center),
				color = :black,
				rotation = π/2,
			)
		end
		
		# Label chip
		text!(
			"$(chip)",
			position = 0.5 .* (40, 80),
			color = :yellow,
			align = (:left, :baseline),
		)
		chip_idx += 1
		push!(axs, ax)
		push!(hms, hm)
	end
	axs = reshape(axs, 1, 2)
	
	Colorbar(fig[:, end+1], hms[1], width=20, label="Counts",)
	
	linkaxes!(axs...)
	hidedecorations!.(axs)
	
	path = "../../ACCESS_WASP-50b/figures/frames"
	mkpath(path)
	save("$(path)/sci_ldss3.png", fig)
	
	fig #|> as_svg
end

# ╔═╡ 4480ae72-3bb2-4e17-99be-28afc756332a
md"""
## Notebook setup
"""

# ╔═╡ 6000db3d-0798-4f76-be31-617d43406b54
#html"""
# <style>
# #launch_binder {
# 	display: none;
# }
# body.disable_ui main {
# 		max-width : 95%;
# 	}
# @media screen and (min-width: 1081px) {
# 	body.disable_ui main {
# 		margin-left : 10px;
# 		max-width : 72%;
# 		align-self: flex-start;
# 	}
# }
# </style>
# """

# ╔═╡ Cell order:
# ╟─fb39c593-86bd-4d4c-b9ec-e5e212a4de98
# ╟─90845d70-35d9-402d-8936-74936b069577
# ╟─8b9581db-71c6-42b6-915b-bde307755bcd
# ╠═fb6e6221-8136-44e2-979b-ecbbd71f740d
# ╠═b8de138b-282a-44b3-be8c-7fea2c88030d
# ╠═c83c2161-95e2-4d08-9934-6d9c12c42a44
# ╟─0d42f6f9-d789-46a3-9e9a-381dbed2d5a5
# ╠═86b4ace1-a891-41e5-a29f-f7eee5f8fb17
# ╠═fed5dd71-98e8-4212-b574-6150bccce2f7
# ╠═6732965f-9ab6-4b84-9703-c212bb89d353
# ╠═26feb668-4e7e-4a9d-a2b6-a5dac81e3ab7
# ╠═be8399c7-d817-476a-bf26-becaea8b6522
# ╠═d00a2334-bee8-4819-b11b-17a5c7eb84e0
# ╠═3a6ab0c0-ba08-4151-9646-c19d45749b9f
# ╟─06a834f0-8c90-4013-af34-725166970969
# ╠═0e66d467-1098-46dc-8d06-36d488b14637
# ╠═5c6e6f7b-70e0-49a8-b064-60dcf1440223
# ╠═c488270a-3126-4e38-a0c8-ee242115a3ea
# ╠═83a9357d-836b-4cee-a41f-eabc8f3f12e7
# ╠═71ba9181-90e4-4d12-97c0-462b3f1df077
# ╟─4480ae72-3bb2-4e17-99be-28afc756332a
# ╠═3433ed02-c27c-4fe5-bfda-a5108a58407c
# ╠═6000db3d-0798-4f76-be31-617d43406b54
