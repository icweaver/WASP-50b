### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ d51d16f0-9c03-11ec-3d03-c15ca2178c81
begin
	import Pkg
	Pkg.activate(Base.current_project())
	
	using PlutoUI
	using MarkdownLiteral: @mdx
end

# ╔═╡ 12ab0900-2bca-45da-a4b6-a317e67660fa
md"""
## References
"""

# ╔═╡ 007a54e5-508f-4106-a87b-7c703f3c80a7
begin
	fsed = "``f_\\mathrm{sed}``"

	Kzz = "``K_\\mathrm{zz}``"

	g = "``g``"

	Teq = "``T_\\mathrm{eq}``"
	
	AM_2001 = "[Ackerman & Marley 2001](https://ui.adsabs.harvard.edu/abs/2001ApJ...556..872A/abstract)"
	
	R_1978 = "[Rossow (1978)](https://ui.adsabs.harvard.edu/abs/1978Icar...36....1R/abstract)"

	OO_2017 = "[Ohno & Okuzumi (2017)](https://ui.adsabs.harvard.edu/abs/2017ApJ...835..261O/abstract)"

	GMA_2018 = "[Gao, Marley, & Ackerman (2018)](https://ui.adsabs.harvard.edu/abs/2018ApJ...855...86G/exportcitation)"
end;

# ╔═╡ 3670f3a7-e01a-4a4f-bf96-4c5f5b3589ef
@mdx """
## $(OO_2017)

* Present a cloud model that takes the microphysics of condensation and coalescence into account, unlike the previous approach by $(AM_2001), which wraps these two processes into $(fsed) as a single free parameter.

* They find that coalescence is driven by differential settling under gravity, which is essential for understanding the later stages of cloud formation, where the cloud particles grow to larger than 10 μm in radius $(R_1978).
"""

# ╔═╡ 9ea4de5d-c63b-4157-807c-b15d038b52fa
@mdx """
## $(GMA_2018)

* $(fsed) study for BDs/sub-stellar objects, taking into account microphysical processes such as cloud particle nucletion rate and condensation nuclei flux


* find that $(fsed) is strongly dependent on (in descending order of impact):

1) material properties and formation pathways of clouds

2) vertical mixing, as parameterized by $(Kzz) (eddy diffusivity)

3) only slightly dependent on g given constant $(Kzz), but a more realistic representation of $(Kzz) that scales with $(g) (e.g., $(AM_2001)) may tell a different story

* Predict that $(fsed) ∝ depth of cloud base (i.e., decreasing $(Teq)) due to higher abundances of condenstate vapor at depth → aiding growth of larger particles
"""

# ╔═╡ 4e112c5e-6f49-4f04-9e85-8365a344e8b3
TableOfContents()

# ╔═╡ Cell order:
# ╟─3670f3a7-e01a-4a4f-bf96-4c5f5b3589ef
# ╟─9ea4de5d-c63b-4157-807c-b15d038b52fa
# ╠═12ab0900-2bca-45da-a4b6-a317e67660fa
# ╠═007a54e5-508f-4106-a87b-7c703f3c80a7
# ╟─4e112c5e-6f49-4f04-9e85-8365a344e8b3
# ╟─d51d16f0-9c03-11ec-3d03-c15ca2178c81
