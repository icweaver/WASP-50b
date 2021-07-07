### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ ab57c20d-542c-4079-8d7d-21e523f7bdaa
begin
	using PyCall, Conda
	Conda.add(["batman-package"]; channel="conda-forge")
end

# ╔═╡ 7c4c556c-7399-4585-a9c4-428734b9a9b2
md"""
# Fun with 🐍

Interfacing with Python is very seamless, thanks to the wonderful package *PyCall.jl*. This defines the string macro `py""`, which we can use to wrap our native Python code and then call it from Julia. There is [A BUNCH](https://github.com/JuliaPy/PyCall.jl#readme) that we can configure here, but for now let's just go with the default setup, which is to install everything into a self contained environment inside of our `~/.julia` folder:
"""

# ╔═╡ ca74eea8-26de-4ecf-9e08-75658a3ae56a
md"""!!! note
	We have also used the [`Conda.jl`](https://github.com/JuliaPy/Conda.jl) package here to install `batman` from the conda-forge channel. If there is an issue running the above cell, try un-plugging it and plugging it back in by running the following in a fresh REPL (which can be copy-and-pasted at once):

	```julia
	using Pkg
	Pkg.add("Conda")
	ENV["PYTHON"] = ""
	Pkg.add("PyCall")
	Pkg.build("PyCall")
	using Conda
	Conda.add(["batman-package"]; channel="conda-forge")
	```
	Similarly, `PyCall` can be pointed to an exisiting conda environment as well:


	```julia	
	using Pkg
	ENV["PYTHON"] = "/home/mango/miniconda3/envs/batman/bin/python"
	Pkg.add("PyCall")
	Pkg.build("PyCall")
	using PyCall
	```
"""

# ╔═╡ 793a51ba-f589-42f4-b562-bfcc0291243d
md"""
Now let's see it in action by using numpy to define a function `neg_norm` that takes the inputs ``(x, y)`` and returns the negative of its vector norm to return:

```math
	-\sqrt{x^2 + y^2} \quad.
```
"""

# ╔═╡ 2be60363-fed7-4146-ada6-bb287ab94678
begin
	py"""
	import numpy as np
	import juliet
	
	def neg_norm(x, y):
		return -np.linalg.norm([x, y])
	"""
	neg_norm(x, y) = py"neg_norm"(x, y)
end

# ╔═╡ bb7aaa9e-1b22-4706-a582-dcb090c85b99
neg_norm(3, 4)

# ╔═╡ 8d512b4c-a0ad-453f-8b51-c09cc63c6a49
md"""
Looks good! Now let's try a module from `batman`:
"""

# ╔═╡ 97b525a0-cbca-4ca7-a415-23a8090eb7ac
begin
	py"""
	import numpy as np
	from batman import _rsky
		
	# Compute distance between centers
	def r_batman(t_0, period, aR_star, incl, ecc, omega):
		return _rsky._rsky(
		np.linspace(0, period, 500),
		t_0,
		period,
		aR_star,
		incl,
		ecc,
		omega,
		1,
		1,
	)
	"""
	r_batman(;t_0, period, aR_star, incl, ecc, omega) = py"r_batman"(
		t_0=t_0, period=period, aR_star=aR_star, incl=incl, ecc=ecc, omega=omega
	)
end

# ╔═╡ abd5ec08-7a8c-45c8-89f6-65d5b34a42cf
md"""
!!! note
	The `f(;arg1, arg2)` syntax is used to explicitly specify keyword arguments in Julia. Similarly, `f(args...)` is Julia's version of `f(**kwargs)`
"""

# ╔═╡ 68649e6a-d26d-405f-9754-e7902407f866
transit_params = (
	t_0 = 0.0,
	period = 2.0,
	aR_star = 7.0,
	incl = 1.5,
	ecc = 0.0,
	omega = 0.0,
)

# ╔═╡ d212844e-b382-4c89-a0f0-809313d0b8db
r = r_batman(;transit_params...)

# ╔═╡ 4cb0f5a5-d448-4d69-a92a-8f044c07d1a2
md"""
Not bad 🍉
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Conda = "8f4d0f93-b110-5947-807f-2305c1781a2d"
PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"

[compat]
Conda = "~1.5.2"
PyCall = "~1.92.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Conda]]
deps = ["JSON", "VersionParsing"]
git-tree-sha1 = "299304989a5e6473d985212c28928899c74e9421"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.5.2"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "81690084b6198a2e1da36fcfda16eeca9f9f24e4"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.1"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "6a8a2a625ab0dea913aba95c11370589e0239ff0"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.6"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "c8abc88faa3f7a3950832ac5d6e690881590d6dc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "1.1.0"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "169bb8ea6b1b143c5cf57df6d34d022a7b60c6db"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.92.3"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[VersionParsing]]
git-tree-sha1 = "80229be1f670524750d905f8fc8148e5a8c4537f"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.2.0"
"""

# ╔═╡ Cell order:
# ╟─7c4c556c-7399-4585-a9c4-428734b9a9b2
# ╠═ab57c20d-542c-4079-8d7d-21e523f7bdaa
# ╟─ca74eea8-26de-4ecf-9e08-75658a3ae56a
# ╠═793a51ba-f589-42f4-b562-bfcc0291243d
# ╠═2be60363-fed7-4146-ada6-bb287ab94678
# ╠═bb7aaa9e-1b22-4706-a582-dcb090c85b99
# ╟─8d512b4c-a0ad-453f-8b51-c09cc63c6a49
# ╠═97b525a0-cbca-4ca7-a415-23a8090eb7ac
# ╟─abd5ec08-7a8c-45c8-89f6-65d5b34a42cf
# ╠═68649e6a-d26d-405f-9754-e7902407f866
# ╠═d212844e-b382-4c89-a0f0-809313d0b8db
# ╟─4cb0f5a5-d448-4d69-a92a-8f044c07d1a2
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
