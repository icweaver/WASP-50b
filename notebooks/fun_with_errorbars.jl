### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 752a6443-4553-45ef-83f1-5de71aa5ff12
using PlutoUI # For widgets and such

# ╔═╡ 55c53d34-247a-4ea7-8649-978f5460b702
using PyCall # For calling Python code in Julia

# ╔═╡ 018625fc-9732-49af-a45b-9903278a17d1
using Measurements

# ╔═╡ 4ccb34df-4b68-4849-bab9-6b4864c2cbf0
using Plots

# ╔═╡ af390593-7842-4823-beb9-874aa0d41377
begin
	using Unitful # The base package with unit support
	using UnitfulAstro # Extends the Unitful package to add astro unit support
	using UnitfulRecipes # Extends the Plots package to handle units
end

# ╔═╡ f6bfe41e-fb59-4dbf-bc03-27f8741d0670
md"""
# Fun with errorbars

$(TableOfContents())
"""

# ╔═╡ 3a77788c-ef4e-4d5b-a388-94ed7f83e597
md"""
## Python function 🐍
"""

# ╔═╡ 7cededae-1a31-4d7e-8f4a-7c20ed68a935
md"""
Let's wrap Victoria's hand-rolled Python code in a Julia function to make comparing results a bit easier:
"""

# ╔═╡ f8b62dee-747f-4303-a085-d8e7a5d60686
begin
	py"""
	import numpy as np
	
	def weighted_avg(values, errors):
		avg = np.sum(values/(errors**2)) / np.sum(1/(errors**2))
		avg_unc = np.sqrt(1/(np.sum(errors**-2)))
		return(avg, avg_unc)
	"""
	weighted_avg(values, errors) = py"weighted_avg"(values, errors)
end

# ╔═╡ 17435d61-d03a-4a74-8a40-e4dc46efa1ea
md"""
## Julia function 🍉
"""

# ╔═╡ cccb5aec-9d90-4c49-99e5-a92183df5329
md"""
Cool, now let's use a pure Julia function to check our work:
"""

# ╔═╡ 83a5deea-69c0-42b0-95e4-ef7d32c47dcc
md"""
That's it 🌈
"""

# ╔═╡ cfcab75c-187a-4cf3-b6d1-0d63256c36ab
md"""
## Comparison 📐
"""

# ╔═╡ 320f2b31-e958-4c47-b9cc-1421bbcb1d05
md"""
First let's generate some fake data:
"""

# ╔═╡ 0dddde61-d668-4dbf-8c25-0c3bdbad27b2
vals = rand(100)

# ╔═╡ 97d71c47-7f79-4537-9463-8901f1d0ed6c
errs = rand(100)

# ╔═╡ 703e6970-4d3f-4f12-8c2f-233c8b9ee9a4
md"""
### `Measurements` aside
"""

# ╔═╡ 2f7d4bbd-23ea-420f-938d-0ab197ae0f74
md"""
In python, we pass these as separate arguments to `weighted_avg`. In Julia, we can do it in one shot by combining the `vals` and `errs` into their own array of `val` ± `err` entries. For example, this can be written verbatim as:

```julia
3 ± 4 # 3 \pm <tab> 4 on the keyboard
```

and Julia will know that we are not dealing with a regular ol' integer or float anymore, but a bonafide new type called a `Measurement`, which we can do all sorts of cool things with, thanks to Julia's extensible and composable ecosystem.
"""

# ╔═╡ d1495894-034f-4afc-b2e8-e4b0c157abb0
md"""
!!! note
	Feel free to try it out in the cell below!
"""

# ╔═╡ 31597114-4a30-4a2c-ab8a-3f33bddbb996


# ╔═╡ 4f073fdd-7d33-46e3-9d2b-39f1b513777b
md"""
Creating an array of `Measurements` from our separate `vals` and `errs` is then as easy as:
"""

# ╔═╡ e187a1af-7008-4015-a2a0-7cf2ff3c0e7f
data = vals .± errs

# ╔═╡ 92d29c9f-e267-4a08-9f90-424c90e84286
md"""
Where the `.` before the `±` operator tells Julia to broadcast things element-wise.

One thing we can do is call the function `weightedmean` that was exported from the `Measurements` package above to perform our weighted averaging. Both the manual Python implementation above and the baked-in Julia implementation use [inverse-variance weighting](https://en.wikipedia.org/wiki/Inverse-variance_weighting), so let's give them both a spin:
"""

# ╔═╡ b00261b8-9d42-4828-81a9-566e9cd31614
md"""
### Results
"""

# ╔═╡ a1164ef1-19b9-4c61-b960-3e31647d6d20
md"""
First we run the python function:
"""

# ╔═╡ 4fa326f9-8d90-4647-886e-19eb1a16b24e
data_avg_python = weighted_avg(vals, errs)

# ╔═╡ 322e4a83-bfcf-429b-aafc-e851501fa7ea
md"""
And next the Julia function:
"""

# ╔═╡ 56d600dc-d9fe-4cf1-9884-c8f4ce7ff8ef
data_avg_julia = weightedmean(data)

# ╔═╡ 86d4b449-638a-4025-a895-4d2513e1d479
md"""
Looks pretty close, but the auto-rounding in the `Measurements` display above makes it a bit difficult to be sure. We can display all decimal places by calling the relevant fields in the `Measurement` object `data_avg_julia`:
"""

# ╔═╡ 93bf5aae-6040-42f3-b041-359b35919621
(data_avg_julia.val, data_avg_julia.err)

# ╔═╡ 8da27f5f-926a-48b3-b8ca-d94fc3a41bbe
md"""
Nice! Julia even has a handly little `≈` operator (typed `\approx<TAB>`) that we can use to confirm equality (up to floating point precision):
"""

# ╔═╡ 48f1e77d-7b55-46d1-89ef-35efc254111d
(data_avg_julia.val, data_avg_julia.err) .≈ data_avg_python

# ╔═╡ 5154f8ad-937a-4b15-8339-d9d129d4b526
md"""
Solid ✔
"""

# ╔═╡ 82c22e46-2298-4596-b20d-032fc8c9dcb9
md"""
## Plots 📊
"""

# ╔═╡ 6ae30b4b-065c-4015-9249-d5bd9b829aa2
md"""
Thanks to Julia's multiple dispatch paradigm, we can use the general `plot` method from the `Plots` package to plot our `Measurements` type data with errorbars automatically included 🌈
"""

# ╔═╡ ac3709fb-fde5-4c6c-a7e0-b9244c655c78
plotly() # nice backend for interactive plots

# ╔═╡ cd6ce0e9-a219-48d7-963d-c3c8a7003d5b
plot(data, xlabel="x data", ylabel="ydata", markerstrokecolor=:auto)

# ╔═╡ 5b49800e-445e-4fde-b74e-1c398e925858
md"""
Speaking of composability, we can even combine our `Measurements` data with units and have those automatically displayed as well!
"""

# ╔═╡ 356b49dc-f631-4566-b05d-c105e56e2933
md"""
We can then make objects that have uncertainties *and* units by doing:
"""

# ╔═╡ 320f2d1c-a1e1-4e9b-b1c8-17f2a818c2b6
x = (3 ± 4)u"Rsun"

# ╔═╡ 4d245c55-0606-43eb-a214-056c5fe6c02b
md"""
Putting it all together:
"""

# ╔═╡ 8f02bd06-12ce-4308-b474-857173abf245
data_with_units = data * u"Rsun"

# ╔═╡ 6ea52901-06df-46ae-9bf9-35c0414a05b6
plot(data_with_units, xlabel="x data", ylabel="y data", markerstrokecolor=:auto)

# ╔═╡ Cell order:
# ╟─f6bfe41e-fb59-4dbf-bc03-27f8741d0670
# ╠═752a6443-4553-45ef-83f1-5de71aa5ff12
# ╟─3a77788c-ef4e-4d5b-a388-94ed7f83e597
# ╟─7cededae-1a31-4d7e-8f4a-7c20ed68a935
# ╠═55c53d34-247a-4ea7-8649-978f5460b702
# ╠═f8b62dee-747f-4303-a085-d8e7a5d60686
# ╟─17435d61-d03a-4a74-8a40-e4dc46efa1ea
# ╟─cccb5aec-9d90-4c49-99e5-a92183df5329
# ╠═018625fc-9732-49af-a45b-9903278a17d1
# ╟─83a5deea-69c0-42b0-95e4-ef7d32c47dcc
# ╟─cfcab75c-187a-4cf3-b6d1-0d63256c36ab
# ╟─320f2b31-e958-4c47-b9cc-1421bbcb1d05
# ╠═0dddde61-d668-4dbf-8c25-0c3bdbad27b2
# ╠═97d71c47-7f79-4537-9463-8901f1d0ed6c
# ╟─703e6970-4d3f-4f12-8c2f-233c8b9ee9a4
# ╟─2f7d4bbd-23ea-420f-938d-0ab197ae0f74
# ╟─d1495894-034f-4afc-b2e8-e4b0c157abb0
# ╠═31597114-4a30-4a2c-ab8a-3f33bddbb996
# ╟─4f073fdd-7d33-46e3-9d2b-39f1b513777b
# ╠═e187a1af-7008-4015-a2a0-7cf2ff3c0e7f
# ╟─92d29c9f-e267-4a08-9f90-424c90e84286
# ╟─b00261b8-9d42-4828-81a9-566e9cd31614
# ╟─a1164ef1-19b9-4c61-b960-3e31647d6d20
# ╠═4fa326f9-8d90-4647-886e-19eb1a16b24e
# ╟─322e4a83-bfcf-429b-aafc-e851501fa7ea
# ╠═56d600dc-d9fe-4cf1-9884-c8f4ce7ff8ef
# ╟─86d4b449-638a-4025-a895-4d2513e1d479
# ╠═93bf5aae-6040-42f3-b041-359b35919621
# ╟─8da27f5f-926a-48b3-b8ca-d94fc3a41bbe
# ╠═48f1e77d-7b55-46d1-89ef-35efc254111d
# ╟─5154f8ad-937a-4b15-8339-d9d129d4b526
# ╟─82c22e46-2298-4596-b20d-032fc8c9dcb9
# ╟─6ae30b4b-065c-4015-9249-d5bd9b829aa2
# ╠═4ccb34df-4b68-4849-bab9-6b4864c2cbf0
# ╠═ac3709fb-fde5-4c6c-a7e0-b9244c655c78
# ╠═cd6ce0e9-a219-48d7-963d-c3c8a7003d5b
# ╟─5b49800e-445e-4fde-b74e-1c398e925858
# ╠═af390593-7842-4823-beb9-874aa0d41377
# ╟─356b49dc-f631-4566-b05d-c105e56e2933
# ╠═320f2d1c-a1e1-4e9b-b1c8-17f2a818c2b6
# ╟─4d245c55-0606-43eb-a214-056c5fe6c02b
# ╠═8f02bd06-12ce-4308-b474-857173abf245
# ╠═6ea52901-06df-46ae-9bf9-35c0414a05b6
