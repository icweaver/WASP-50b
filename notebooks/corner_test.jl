### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 693a7390-d525-11eb-17dd-5d4965c0a82a
using AlgebraOfGraphics, CairoMakie, BenchmarkTools, PlutoUI

# ╔═╡ 2a13eb6f-2ad0-4de2-8707-245592e7b85d
function create_corner(;n_params=3, hide_decorations=false)	
	# Create empty corner plot grid
	fig = Figure(resolution=(1_400, 1_400))
	
	for j in 1:n_params, i in 1:n_params
		# Create subplot apply global settings
		ax = Axis(fig[i, j];
			aspect = 1,
			xticklabelrotation = π/4,
			xticks = LinearTicks(3),
			yticks = LinearTicks(3),
			limits = ((-14, 14), (-14, 14)),
		)
		
		if hide_decorations
			# Hide upper triangle
			j > i && (hidedecorations!(ax); hidespines!(ax))
			# Hide y ticks on diagonals
			j == i && (hideydecorations!(ax))
			# Hide x ticks on all diagonal elements except the bottom one
			j == i && i != n_params && (
				hidexdecorations!(ax, grid=false);
			)
			# Hide ticks on interior lower triangle
			j < i && i != n_params && j != 1 && (
				hideydecorations!(ax, grid=false);
				hidexdecorations!(ax, grid=false);
			)
			# Hide remaining xyticks
			j < i && j == 1 && i != n_params && (
				hidexdecorations!(ax, grid=false);
			)
			j < i && i == n_params && j != 1 && (
				hideydecorations!(ax, grid=false);
			)
		end
	end
	
	return fig
end

# ╔═╡ 10ca5e76-9080-4e3a-87ae-9c95cae2deae
with_terminal() do
	@btime create_corner(n_params=3, hide_decorations=false)
	@btime create_corner(n_params=3, hide_decorations=true)
end

# ╔═╡ d31092c2-736e-466b-a5f8-5dedc34381ad
create_corner(n_params=3, hide_decorations=false)

# ╔═╡ 5ef5e798-3432-4de7-aeaa-428c53b532a6
create_corner(n_params=3, hide_decorations=true)

# ╔═╡ 623cba18-37d8-4d3c-ba29-731098dca33d
begin
	set_aog_theme!()
	update_theme!(
		Theme(
			Axis = (xlabelsize=18, ylabelsize=18,),
			Label = (textsize=18,),
			Lines = (linewidth=3,),
			Scatter = (linewidth=10,),
			fontsize = 18,
			rowgap = 0,
			colgap = 0,
		)
	)
end

# ╔═╡ Cell order:
# ╠═10ca5e76-9080-4e3a-87ae-9c95cae2deae
# ╠═d31092c2-736e-466b-a5f8-5dedc34381ad
# ╠═5ef5e798-3432-4de7-aeaa-428c53b532a6
# ╠═2a13eb6f-2ad0-4de2-8707-245592e7b85d
# ╠═623cba18-37d8-4d3c-ba29-731098dca33d
# ╠═693a7390-d525-11eb-17dd-5d4965c0a82a
