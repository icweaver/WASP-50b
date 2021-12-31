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
end

# ╔═╡ fb39c593-86bd-4d4c-b9ec-e5e212a4de98
md"""
# Raw data

In this notebook we will view bias-subtracted sample frames from each instrument to verify the mask slit placement, as well as the positions of WASP-50 and its comparison stars on each CCD chip.

For each fits files (1 per chip), we extract the region of each chip dedicated to main data collection (`DATASEC`), and then subtract out the median bias level measured in the overscan region (`BIASSEC`). We then store the subtracted data in a data cube for each instrument.

$(TableOfContents())

!!! note "Data download"
	```
	 rclone sync -P ACCESS_box:WASP-50b/data/raw_data data/raw_data
	```
"""

# ╔═╡ 90845d70-35d9-402d-8936-74936b069577
md"""
## IMACS 🎱
Starting with IMACS, let's first select the night we would like to visualize from the dropdown menu below:
"""

# ╔═╡ 8b9581db-71c6-42b6-915b-bde307755bcd
@bind DATA_DIR_IMACS Select(glob("data/raw_data/IMACS/ut*"))

# ╔═╡ fb6e6221-8136-44e2-979b-ecbbd71f740d
df_sci_IMACS = fitscollection(DATA_DIR_IMACS, abspath=false)

# ╔═╡ 0d42f6f9-d789-46a3-9e9a-381dbed2d5a5
md"""
We next compute the bias-subtracted science frame for each of these $(nrow(df_sci_IMACS)) chips:
"""

# ╔═╡ c83c2161-95e2-4d08-9934-6d9c12c42a44
process_frames(df) = @views map(ccds(df)) do img
	img_data = img[fits_indices(img.hdr["DATASEC"])...]
	img_bias = img[fits_indices(img.hdr["BIASSEC"])...]
	corr = img_data .- median(img_bias)
end

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

# ╔═╡ 4b247772-ca20-480d-b649-143d3489ce46
f(a; f_kwargs=()) = (a, f_kwargs...)

# ╔═╡ ab0a289c-21e4-4dee-9f70-4b3a6a5d67de
f(3, f_kwargs=(b = 5,))

# ╔═╡ 0e66d467-1098-46dc-8d06-36d488b14637
@bind DATA_DIR_LDSS3 Select(glob("data/raw_data/LDSS3/ut*"))

# ╔═╡ 5c6e6f7b-70e0-49a8-b064-60dcf1440223
df_sci_LDSS3 = fitscollection(DATA_DIR_LDSS3, abspath=false)

# ╔═╡ 06a834f0-8c90-4013-af34-725166970969
md"""
## LDSS3 2️⃣

We follow the same operations to visualize the $(nrow(df_sci_LDSS3)) chips for LDSS3 below.
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
	stepsize = 32,
	hm_kwargs = (),
)
	hm = plot!(ax, img;
		colormap = :magma,
		hm_kwargs...,
	)
	
	# Add labels
	chip = "c$i"
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

# ╔═╡ 3a6ab0c0-ba08-4151-9646-c19d45749b9f
let
	fig = Figure(resolution = (800, 600))
	hm = nothing
	stepsize = 3
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

		hm = plot_frame!(ax, img, ch, coords_IMACS;
			stepsize = stepsize,
			hm_kwargs = (colorrange=(0, 2500),),
		)
	end
	
	Colorbar(fig[2:3, end+1], hm, width=20, label="Counts",)
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
	stepsize = 1
	hm = nothing
	for j ∈ 1:2
		ax = Axis(fig[1, j])
		img = @view(frames_LDSS3[j][begin:stepsize:end, begin:stepsize:end])'
		hm = plot_frame!(ax, img, j, coords_LDSS3;
			stepsize = stepsize,
			hm_kwargs = (colorrange = (0, 2500),),
		)
		Label(fig[2, j], "c$(j)", tellwidth=false)
	end

	Colorbar(fig[1, end+1], hm, width=20, label="Counts",)
	
	axs = filter(x -> x isa Axis, fig.content)
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

# ╔═╡ db4a4cd8-c5e8-4124-935f-0666f6e73fe2
begin
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

# ╔═╡ 6000db3d-0798-4f76-be31-617d43406b54
html"""
<style>
#launch_binder {
	display: none;
}
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
# ╟─fb39c593-86bd-4d4c-b9ec-e5e212a4de98
# ╟─90845d70-35d9-402d-8936-74936b069577
# ╟─8b9581db-71c6-42b6-915b-bde307755bcd
# ╠═fb6e6221-8136-44e2-979b-ecbbd71f740d
# ╟─0d42f6f9-d789-46a3-9e9a-381dbed2d5a5
# ╠═c83c2161-95e2-4d08-9934-6d9c12c42a44
# ╠═86b4ace1-a891-41e5-a29f-f7eee5f8fb17
# ╟─27b793f3-7a7c-48ad-8302-deffa2dd017b
# ╠═26feb668-4e7e-4a9d-a2b6-a5dac81e3ab7
# ╠═3a6ab0c0-ba08-4151-9646-c19d45749b9f
# ╠═4b247772-ca20-480d-b649-143d3489ce46
# ╠═ab0a289c-21e4-4dee-9f70-4b3a6a5d67de
# ╠═bf8ef5a9-0806-44b4-907d-c95d6926dabb
# ╟─06a834f0-8c90-4013-af34-725166970969
# ╟─0e66d467-1098-46dc-8d06-36d488b14637
# ╠═5c6e6f7b-70e0-49a8-b064-60dcf1440223
# ╠═c488270a-3126-4e38-a0c8-ee242115a3ea
# ╠═83a9357d-836b-4cee-a41f-eabc8f3f12e7
# ╠═71ba9181-90e4-4d12-97c0-462b3f1df077
# ╟─4480ae72-3bb2-4e17-99be-28afc756332a
# ╠═db4a4cd8-c5e8-4124-935f-0666f6e73fe2
# ╠═3433ed02-c27c-4fe5-bfda-a5108a58407c
# ╟─6000db3d-0798-4f76-be31-617d43406b54
