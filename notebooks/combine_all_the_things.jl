### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 45c6d6f6-86c0-11ec-017f-95b00d819374
begin
	import Pkg
	Pkg.activate(Base.current_project())
	
	using JLD2, FileIO, CairoMakie, AlgebraOfGraphics, Statistics
end

# ╔═╡ 3e3b35da-8aff-492a-b3db-6d07fe418577
LC_spectra = load("data/reduced_data/cubes/LC_spectra.jld2");

# ╔═╡ 63ed64d0-8b5f-4873-a098-ef572c1ef602
LC_spectra["LDSS3C_ut150927_flat"]

# ╔═╡ 452cd5af-c17c-4eab-b07a-c60939e4cdeb
x = [1, 2, 3 ,4]

# ╔═╡ 7db2747a-4616-4733-a169-d851772625a7
function spec_plot!(ax, wav, μ, σ; color=:blue, norm=1.0, label="")
	band!(ax, wav, μ .- σ, μ .+ σ, color=(color, 0.25))
	lines!(ax, wav, μ; color, label) #cycle = Cycle(:linestyle))
end

# ╔═╡ 008ee0a0-9a98-4442-a43c-bb30a3960057
med_std(A; dims=1) = (median(A; dims), std(A; dims)) .|> vec

# ╔═╡ ccbfeab9-a4cf-460f-88c5-5298824f3523
begin
	##############
	# PLOT CONFIGS
	##############
	const FIG_TALL = 72 .* (6, 8)
	const FIG_WIDE = 72 .* (12, 6)
	const FIG_LARGE = 72 .* (12, 12)
	const COLORS_SERIES = to_colormap(:seaborn_colorblind, 9)
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
			Scatter = (linewidth=10,),
			Text = (font = AlgebraOfGraphics.firasans("Regular"), textsize=18),
			palette = (color=COLORS, patchcolor=[(c, 0.35) for c in COLORS]),
			figure_padding = 1.5,
			fontsize = 18,
			rowgap = 5,
			colgap = 5,
		)
	)

	COLORS
end

# ╔═╡ 77d5b506-6887-4bb0-aed9-61cc9af88630
function plot_specs!(ax, LC_spectra; norm=1.0, label="")
	wav = LC_spectra["wav"]
	for (i, (label, f)) in enumerate(sort(LC_spectra["spec"]))
		μ, σ = f
		spec_plot!(ax, wav, μ, σ; color=COLORS_SERIES[i], norm, label)
	end
end

# ╔═╡ a0d9d9b3-427d-4903-bbcf-9877c4a9d4ea
let
	fig = Figure(resolution=FIG_LARGE)
	grid = CartesianIndices((2, 2))
	axs = []
	for (i, (transit, cube)) ∈ enumerate(LC_spectra)
		g_idxs = grid[i]
		ax = Axis(fig[g_idxs.I...])
		push!(axs, ax)
		plot_specs!(ax, cube)
	end
	
	axs = reshape(axs, (2, 2))
	linkaxes!.(axs...)
	hidexdecorations!.(axs[1, :])
	hideydecorations!.(axs[:, 2])
	
	fig
end

# ╔═╡ Cell order:
# ╠═45c6d6f6-86c0-11ec-017f-95b00d819374
# ╠═3e3b35da-8aff-492a-b3db-6d07fe418577
# ╠═63ed64d0-8b5f-4873-a098-ef572c1ef602
# ╠═a0d9d9b3-427d-4903-bbcf-9877c4a9d4ea
# ╠═452cd5af-c17c-4eab-b07a-c60939e4cdeb
# ╠═77d5b506-6887-4bb0-aed9-61cc9af88630
# ╠═7db2747a-4616-4733-a169-d851772625a7
# ╠═008ee0a0-9a98-4442-a43c-bb30a3960057
# ╠═ccbfeab9-a4cf-460f-88c5-5298824f3523
