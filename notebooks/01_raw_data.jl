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
	using CCDReduction: fits_indices
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

In this notebook we will view bias-subtracted sample frames from each instrument to verify the mask slit placement, as well as the positions of WASP-50 and its comparison stars on each CCD chip.

For each fits files (1 per chip), we extract the region of each chip dedicated to main data collection (`DATASEC`), and then subtract out the median bias level measured in the overscan region (`BIASSEC`). We then store the subtracted data in a data cube for each instrument.

$(TableOfContents())

!!! note "Data download"
	```
	 rclone sync ACCESS_box:WASP-50b/data/raw data/raw
	```
"""

# ╔═╡ 90845d70-35d9-402d-8936-74936b069577
md"""
## IMACS 🎱
Starting with IMACS, let's first select the night we would like to visualize from the dropdown menu below:
"""

# ╔═╡ 8b9581db-71c6-42b6-915b-bde307755bcd
@bind DATA_DIR_IMACS Select(glob("data/raw/IMACS/ut*"))

# ╔═╡ fb6e6221-8136-44e2-979b-ecbbd71f740d
df_sci_IMACS = fitscollection(DATA_DIR_IMACS, abspath=false)

# ╔═╡ b8de138b-282a-44b3-be8c-7fea2c88030d
img_overscan_corrected = @views map(ccds(df_sci_IMACS)) do img
	img_data = img[CCDReduction.fits_indices(img.hdr["DATASEC"])...]
	img_bias = img[CCDReduction.fits_indices(img.hdr["BIASSEC"])...]
	corr = img_data .- median(img_bias)
end

# ╔═╡ c83c2161-95e2-4d08-9934-6d9c12c42a44
process_frames(df) = @views map(ccds(df)) do img
	img_data = img[fits_indices(img.hdr["DATASEC"])...]
	img_bias = img[fits_indices(img.hdr["BIASSEC"])...]
	corr = img_data .- median(img_bias)
end

# ╔═╡ 0d42f6f9-d789-46a3-9e9a-381dbed2d5a5
md"""
We next compute the bias-subtracted science frame for each of these $(nrow(df_sci_IMACS)) chips:
"""

# ╔═╡ 86b4ace1-a891-41e5-a29f-f7eee5f8fb17
frames_IMACS = process_frames(df_sci_IMACS)

# ╔═╡ 27b793f3-7a7c-48ad-8302-deffa2dd017b
md"""
And finally, we overlay the coordinates of the science mask onto the CCD array:
"""

# ╔═╡ 26feb668-4e7e-4a9d-a2b6-a5dac81e3ab7
coords_IMACS = CSV.read("$(DATA_DIR_IMACS)/WASP50.coords", DataFrame;
	delim = ' ',
	ignorerepeated = true,
	header = ["target", "chip", "x", "y"],
);

# ╔═╡ e4a4ee16-a986-4d68-a73f-60b5c10be192
x = CartesianIndices((2, 4))

# ╔═╡ 21514113-1955-4f7c-a188-03d0440b36e2
a, b = x[1].I

# ╔═╡ e5aa8301-d755-47df-a057-214c2acc3594
let
	fig = Figure()
	ax = Axis(fig[1, 1])
	Label(fig[2, 1], "Hey";
		tellwidth=false)

	heatmap!(ax, rand(30, 4), colormap=(:viridis, 1.0))
	fig
end

# ╔═╡ 0e66d467-1098-46dc-8d06-36d488b14637
@bind DATA_DIR_LDSS3 Select(glob("data/raw/LDSS3/ut*"))

# ╔═╡ 5c6e6f7b-70e0-49a8-b064-60dcf1440223
df_sci_LDSS3 = fitscollection(DATA_DIR_LDSS3, abspath=false)

# ╔═╡ 06a834f0-8c90-4013-af34-725166970969
md"""
## LDSS3 2️⃣

We follow the same operations to visualize the $(nrow(df_sci_LDSS3)) science frames below.
"""

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
function plot_frame!(ax, img=frames_LDSS3[1], i=1, coords=coords_LDSS3;
	step = 32,
	colorrange = (0, 2500),
)
	hm = plot!(ax, img;
		colormap = :magma,
		colorrange = colorrange,
	)
	
	# Add labels
	chip = "c$i"
	for obj in eachrow(@subset coords :chip .== chip)
		coord = (obj.x, obj.y) ./ step
		text!(ax, obj.target;
			position = coord,
			textsize = 11,
			align = (:center, :center),
			color = :black,
			rotation = π/2,
		)
	end
	# text!("$(chip)";
	# 	position = 0.5 .* (40, 80),
	# 	color = :yellow,
	# 	align = (:left, :baseline),
	# )

	return hm
end

# ╔═╡ 3a6ab0c0-ba08-4151-9646-c19d45749b9f
let
	fig = Figure(resolution = (800, 600))
	hm = nothing
	step = 64
	grid = CartesianIndices((2, 4))
	chip_order = [1, 6, 2, 5, 3, 8, 4, 7]
	for (g, ch) ∈ zip(grid, chip_order)
		# Set chip orientation based on IMACS conventions
		i, j = g.I
		if ch ∈ 5:8
			xflip, yflip = reverse, reverse
			Label(fig[3, j], "c$(ch)", tellwidth=false)
		else
			xflip, yflip = x -> x, x -> x
			#Label(fig[0, j], "c$(ch)", tellwidth=false)
		end
		ax = Axis(fig[i, j], xreversed=false, yreversed=true)
		img = @view(frames_IMACS[ch][xflip(begin:step:end), yflip(begin:step:end)])'

		hm = plot_frame!(ax, img, ch, coords_IMACS;
			step = step,
			colorrange = (0, 250),
		)
		
		# hm = plot!(
		# 	ax,
		# 	d,
		# 	colormap = :magma,
		# 	colorrange = (0, 500),
		# )

		# # Label objects on chip
		# chip = "c$ch"
		# for obj in eachrow(@subset coords_IMACS :chip .== chip)
		# 	coord = (obj.x, obj.y) ./ step
		# 	text!(ax, split(obj.target, "_")[1];
		# 		position = coord,
		# 		textsize = 11,
		# 		align = (:center, :center),
		# 		color = :white,
		# 		rotation = π/2,
		# 	)
		# end
		
		# # Label chip
		# text!(
		# 	"$(chip)",
		# 	position = (0, ncols),
		# 	color = :yellow,
		# 	align = (:left, :baseline),
		# 	offset = 0.5 .* (20, 20),
		# )
		# push!(hms, hm)
	end
	
	Colorbar(fig[1:2, end+1], hm, width=20, label="Counts",)
	axs = filter(x -> x isa Axis, fig.content)
	linkaxes!(axs...)
	hidedecorations!.(axs)
	
	path = "../../ACCESS_WASP-50b/figures/frames"
	mkpath(path)
	save("$(path)/sci_imacs_$(basename(DATA_DIR_IMACS)).png", fig)
	
	fig #|> as_svg
end

# ╔═╡ 71ba9181-90e4-4d12-97c0-462b3f1df077
let
	fig = Figure(resolution = (800, 600))
	step = 1
	#hms = []

	hm = nothing
	for j ∈ 1:2
		ax = Axis(fig[1, j])
		img = @view(frames_LDSS3[j][begin:step:end, begin:step:end])'
		hm = plot_frame!(ax, img, j, coords_LDSS3, step=1)
	end

	Colorbar(fig[:, end+1], hm, width=20, label="Counts",)
	
	axs = filter(x -> x isa Axis, fig.content)
	linkaxes!(axs...)
	hidedecorations!.(axs)
	
	path = "../../ACCESS_WASP-50b/figures/frames"
	mkpath(path)
	save("$(path)/sci_ldss3.png", fig)
	
	fig #|> as_svg
end

# ╔═╡ 036fb8af-6fca-4cdb-80f3-9ecdad102868
@with_terminal begin
	x = 0
	for _ in 1:10
		x += 1
	end
	println(x)
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
# ╟─27b793f3-7a7c-48ad-8302-deffa2dd017b
# ╠═26feb668-4e7e-4a9d-a2b6-a5dac81e3ab7
# ╠═3a6ab0c0-ba08-4151-9646-c19d45749b9f
# ╠═e4a4ee16-a986-4d68-a73f-60b5c10be192
# ╠═21514113-1955-4f7c-a188-03d0440b36e2
# ╠═bf8ef5a9-0806-44b4-907d-c95d6926dabb
# ╠═e5aa8301-d755-47df-a057-214c2acc3594
# ╟─06a834f0-8c90-4013-af34-725166970969
# ╠═0e66d467-1098-46dc-8d06-36d488b14637
# ╠═5c6e6f7b-70e0-49a8-b064-60dcf1440223
# ╠═c488270a-3126-4e38-a0c8-ee242115a3ea
# ╠═83a9357d-836b-4cee-a41f-eabc8f3f12e7
# ╠═71ba9181-90e4-4d12-97c0-462b3f1df077
# ╠═036fb8af-6fca-4cdb-80f3-9ecdad102868
# ╟─4480ae72-3bb2-4e17-99be-28afc756332a
# ╠═3433ed02-c27c-4fe5-bfda-a5108a58407c
# ╠═6000db3d-0798-4f76-be31-617d43406b54
