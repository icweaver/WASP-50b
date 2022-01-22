### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ 22b67cb0-6687-11ec-1868-bb216a9703f4
begin
	import Pkg
	Pkg.activate(Base.current_project())
	
	using PythonCall, CondaPkg, PlutoUI
end

# ╔═╡ f150770d-cbcd-4eda-a0a5-b759a21f5b9b
md"""
# Fun with 🐍

In this notebook we will be showing some useage examples for the new package `PythonCall.jl`. It behaves similarly to the older package `PyCall.jl`, but with a few key differences that make interacting with Python quite nice:

* automatic, project-specific environments for easy reproducibility
* package management via [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) for blazing fast Python package management 🔥
* easy plot display support

Let's try some things out!

$(TableOfContents())
"""

# ╔═╡ 30181f53-97e2-4d21-83d2-9394ef64aba2
md"""
## `pyeval`
To see how we can start using this package, let's start by creating some simple Python objects:
"""

# ╔═╡ b1f1706d-da3b-48ae-bed1-4f5fb9e95726
py_list = @pyeval "[1, 2, 3]"

# ╔═╡ fa184670-718c-4e68-8e23-fe6e72a651c5
py_dict = @pyeval "{'one': 1, 'two': 2, 'three': 3}"

# ╔═╡ 53b89008-8ac0-4f20-955d-1d875530f121
py_range = @pyeval "range(0, 10, 2)"

# ╔═╡ 81a0a61d-0ff8-4610-915b-3f32a2742f5d
md"""
The `@pyeval` macro behaves like Python's [`eval`](https://docs.python.org/3/library/functions.html#eval) and evaluates whatever Python expression (written as a string literal or `cmd`) is passed to it. The returned objects can participate in the usual Python methods now:
"""

# ╔═╡ 15b2d0ab-bc1a-435a-9557-99537d117d91
py_list.append(4); py_list

# ╔═╡ ffa07e57-7797-4000-8be4-a166895bcc54
py_dict["three"]

# ╔═╡ 5e0919d8-0a02-46fc-92fe-0a5e40d1708c
py_range[0] # Python's zero-based indexing is automatically understood

# ╔═╡ 0b3322b5-cd15-49d9-a446-35224a4b7637
md"""
!!! tip
	It's recommended to pass strings ("") instead of cmds (\`\`) to play nice with Pluto's `ExpressionExplorer`
 
	Related issue?: <https://github.com/fonsp/Pluto.jl/issues/961>
"""

# ╔═╡ 3359adff-8eb5-4911-a914-047ecf3663fe
md"""
`PythonCall.jl` also provides a function version of `@pyeval` if we would like to do things like string interpolation first:
"""

# ╔═╡ 6be1822a-d905-4de9-8345-4e160e53ea64
let
	four = 4
	pyeval("[1, 2, 3, $four]", Main)
end

# ╔═╡ ed066bda-07b4-4c52-9302-15c66cb4e1d8
md"""
More info about the usage for each version is available in the Live docs.
"""

# ╔═╡ a7e24a65-cfa8-4e03-9664-f8f424324917
md"""
To a certain extent, Julia's functions can also operate on these objects automatically:
"""

# ╔═╡ 1e94e185-af8c-4796-9358-6a2c15c7fd43
py_range |> sum

# ╔═╡ b810a02b-8329-4f51-8cfc-c6920ad1cf5d
py_range |> collect

# ╔═╡ 40cb6860-7ac6-40ca-8b1f-30360fd6aef7
md"""
As an alternative, `PythonCall` makes it very easy to convert these objects to their native Julia counterparts with the `@pyconvert`/`pyconvert` macro/function:
"""

# ╔═╡ 57686014-99e8-4733-8f1d-ba80306b9da9
(@pyconvert StepRange py_range) |> collect

# ╔═╡ 1360029e-76f3-468a-9b2f-d27482bbc525
pyconvert(Dict, py_dict)

# ╔═╡ fde4db08-e5e1-4e92-b1c7-9a41edb04476
md"""
Next, we will take a look at running Python statements.
"""

# ╔═╡ 71ba808a-ebde-43ad-b264-ccf8e676891a
md"""
## `pyexec`

We saw how to run simple Python expressions with `pyeval`. Now we will turn to executing Python statements, which can be used to perform more complex tasks like storing values in variables and defining functions. This is accomplished with `PythonCalls.jl's` `@pyexec`/`pyexec` macro/function, which behave's like Python's [`exec`](https://docs.python.org/3/library/functions.html#exec):
"""

# ╔═╡ 193e65ae-c344-4dda-b92b-e0b136eb7581
@pyexec (x=1, y=2) => "ans = x + y" => ans

# ╔═╡ 0c8e14ef-782d-420b-8906-3a139bacf331
pyexec(@NamedTuple{ans}, "ans = 1 + 2", Main)

# ╔═╡ 1eccfc11-4ea1-4b79-a688-2644d6d1d0fd
md"""
## Combining `pyexec` and `pyeval`

We can now compose these ideas to start interacting with full blocks of Python code:
"""

# ╔═╡ 5b597d2a-2483-4b80-860b-839fa3ddeaec
begin
	@pyexec """
	global greeting
	def greeting(name):
		return f"Hi {name} 👋"
	"""
	greeting(name) = @pyeval("greeting")(name)
end

# ╔═╡ 5aa1ba97-d55c-4794-a839-e90effb84bbe
greeting("Pluto citizen")

# ╔═╡ e9e6c4e0-6834-4b6b-ac20-ff722f9a5cd9
md"""
## Using packages
`PythonCall.jl` has a nice companion package named `CondaPkg.jl`, which we can use to easily install Python packages into an environment in the same directory as this notebook:
"""

# ╔═╡ 0bf45621-7776-403e-b3da-5311a5c30e20
begin
	CondaPkg.add.(("matplotlib", "numpy"))
	CondaPkg.resolve()
end

# ╔═╡ f5d2228c-b45e-4dce-99fb-668f6c50df30
md"""
We now use the `@py` macro to import and interact with these packages:
"""

# ╔═╡ a159c36c-83c8-460e-a7e2-e85c7df8d9da
@py begin
	import matplotlib.pyplot as plt
	import numpy as np
end

# ╔═╡ c83fdb79-b0d9-4830-b6cf-74d2a669ceed
let
	xs = np.random.rand(10, 4, 4)
	
	fig, axes = plt.subplots(2, 2, sharex=true, sharey=true)
	
	for (ax, x) ∈ zip(axes.flat, xs)
		ax.plot(x)
	end

	fig.tight_layout()
	
	plt.gcf()
end

# ╔═╡ 74653590-ff90-4c83-97cd-68a07b61ecd5
md"""
!!! warning "TODO"
	Update to `#main` branch to show inline plotting above
"""

# ╔═╡ 681ece7a-6800-4779-8118-28c3179cd43a
md"""
## Environment details

Where is all this stuff being downloaded/run? Let's see!
"""

# ╔═╡ 655674de-56c1-4386-8fda-aa5c95b6271f
@with_terminal CondaPkg.status()

# ╔═╡ 9f68d823-64ac-48cc-b35e-b131a0bd5c50
@with_terminal run(CondaPkg.MicroMamba.cmd(`list`))

# ╔═╡ 782a5542-860b-4862-9bb9-dc4b1c2f48bb
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
# ╟─f150770d-cbcd-4eda-a0a5-b759a21f5b9b
# ╟─30181f53-97e2-4d21-83d2-9394ef64aba2
# ╠═b1f1706d-da3b-48ae-bed1-4f5fb9e95726
# ╠═fa184670-718c-4e68-8e23-fe6e72a651c5
# ╠═53b89008-8ac0-4f20-955d-1d875530f121
# ╟─81a0a61d-0ff8-4610-915b-3f32a2742f5d
# ╠═15b2d0ab-bc1a-435a-9557-99537d117d91
# ╠═ffa07e57-7797-4000-8be4-a166895bcc54
# ╠═5e0919d8-0a02-46fc-92fe-0a5e40d1708c
# ╟─0b3322b5-cd15-49d9-a446-35224a4b7637
# ╟─3359adff-8eb5-4911-a914-047ecf3663fe
# ╠═6be1822a-d905-4de9-8345-4e160e53ea64
# ╟─ed066bda-07b4-4c52-9302-15c66cb4e1d8
# ╟─a7e24a65-cfa8-4e03-9664-f8f424324917
# ╠═1e94e185-af8c-4796-9358-6a2c15c7fd43
# ╠═b810a02b-8329-4f51-8cfc-c6920ad1cf5d
# ╟─40cb6860-7ac6-40ca-8b1f-30360fd6aef7
# ╠═57686014-99e8-4733-8f1d-ba80306b9da9
# ╠═1360029e-76f3-468a-9b2f-d27482bbc525
# ╟─fde4db08-e5e1-4e92-b1c7-9a41edb04476
# ╟─71ba808a-ebde-43ad-b264-ccf8e676891a
# ╠═193e65ae-c344-4dda-b92b-e0b136eb7581
# ╠═0c8e14ef-782d-420b-8906-3a139bacf331
# ╠═1eccfc11-4ea1-4b79-a688-2644d6d1d0fd
# ╠═5b597d2a-2483-4b80-860b-839fa3ddeaec
# ╠═5aa1ba97-d55c-4794-a839-e90effb84bbe
# ╟─e9e6c4e0-6834-4b6b-ac20-ff722f9a5cd9
# ╠═0bf45621-7776-403e-b3da-5311a5c30e20
# ╟─f5d2228c-b45e-4dce-99fb-668f6c50df30
# ╠═a159c36c-83c8-460e-a7e2-e85c7df8d9da
# ╠═c83fdb79-b0d9-4830-b6cf-74d2a669ceed
# ╟─74653590-ff90-4c83-97cd-68a07b61ecd5
# ╟─681ece7a-6800-4779-8118-28c3179cd43a
# ╠═655674de-56c1-4386-8fda-aa5c95b6271f
# ╠═9f68d823-64ac-48cc-b35e-b131a0bd5c50
# ╠═22b67cb0-6687-11ec-1868-bb216a9703f4
# ╟─782a5542-860b-4862-9bb9-dc4b1c2f48bb
