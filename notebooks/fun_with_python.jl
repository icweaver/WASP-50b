### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ 22b67cb0-6687-11ec-1868-bb216a9703f4
begin
	using PythonCall, CondaPkg, PlutoUI
	import MarkdownLiteral: @mdx
end

# ╔═╡ f150770d-cbcd-4eda-a0a5-b759a21f5b9b
@mdx """
# Fun with 🐍

In this notebook we will be showing some useage examples for the new package `PythonCall.jl`. It behaves similarly to the older package `PyCall.jl`, but with a few key differences that make interacting with Python quite nice:

* automatic, project-specific environments for easy reproducibility
* package management via [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) for blazing fast Python package management 🔥
* easy plot display support

Let's try some things out!
"""

# ╔═╡ acf2a1eb-e99f-442f-86de-5ebe5a2f601c
TableOfContents()

# ╔═╡ 30181f53-97e2-4d21-83d2-9394ef64aba2
@mdx """
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
@mdx """
The `@pyeval` macro behaves like Python's [`eval`](https://docs.python.org/3/library/functions.html#eval) and evaluates whatever Python expression (written as a string literal or `cmd`) is passed to it. The returned objects can participate in the usual Python methods now:
"""

# ╔═╡ 15b2d0ab-bc1a-435a-9557-99537d117d91
py_list.append(4); py_list

# ╔═╡ ffa07e57-7797-4000-8be4-a166895bcc54
py_dict["three"]

# ╔═╡ 5e0919d8-0a02-46fc-92fe-0a5e40d1708c
py_range[0] # Python's zero-based indexing is automatically understood

# ╔═╡ 3359adff-8eb5-4911-a914-047ecf3663fe
@mdx """
`PythonCall.jl` also provides a function version of `@pyeval` if we would like to do things like string interpolation first:
"""

# ╔═╡ 6be1822a-d905-4de9-8345-4e160e53ea64
let
	four = 4
	pyeval("[1, 2, 3, $four]", Main)
end

# ╔═╡ ed066bda-07b4-4c52-9302-15c66cb4e1d8
@mdx """
More info about the usage for each version is available in the Live docs.
"""

# ╔═╡ a7e24a65-cfa8-4e03-9664-f8f424324917
@mdx """
To a certain extent, Julia's functions can also operate on these objects automatically:
"""

# ╔═╡ 1e94e185-af8c-4796-9358-6a2c15c7fd43
py_range |> sum

# ╔═╡ b810a02b-8329-4f51-8cfc-c6920ad1cf5d
py_range |> collect

# ╔═╡ 40cb6860-7ac6-40ca-8b1f-30360fd6aef7
@mdx """
As an alternative, `PythonCall` makes it very easy to convert these objects to their native Julia counterparts with the `@pyconvert`/`pyconvert` macro/function:
"""

# ╔═╡ 57686014-99e8-4733-8f1d-ba80306b9da9
(@pyconvert StepRange py_range) |> collect

# ╔═╡ 1360029e-76f3-468a-9b2f-d27482bbc525
pyconvert(Dict, py_dict)

# ╔═╡ fde4db08-e5e1-4e92-b1c7-9a41edb04476
@mdx """
Next, we will take a look at running Python statements.
"""

# ╔═╡ 71ba808a-ebde-43ad-b264-ccf8e676891a
@mdx """
## `pyexec`

We saw how to run simple Python expressions with `pyeval`. Now we will turn to executing Python statements, which can be used to perform more complex tasks like storing values in variables and defining functions. This is accomplished with `PythonCalls.jl's` `@pyexec`/`pyexec` macro/function, which behave's like Python's [`exec`](https://docs.python.org/3/library/functions.html#exec):
"""

# ╔═╡ 193e65ae-c344-4dda-b92b-e0b136eb7581
@pyexec (x=1, y=2) => "ans = x + y" => ans

# ╔═╡ 0c8e14ef-782d-420b-8906-3a139bacf331
pyexec(@NamedTuple{ans}, "ans = 1 + 2", Main)

# ╔═╡ 1eccfc11-4ea1-4b79-a688-2644d6d1d0fd
@mdx """
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
@mdx """
## Using packages
`PythonCall.jl` has a nice companion package named `CondaPkg.jl`, which we can use to easily install Python packages into an environment in the same directory as this notebook:
"""

# ╔═╡ 0bf45621-7776-403e-b3da-5311a5c30e20
begin
	CondaPkg.add.(("matplotlib", "numpy"))
	CondaPkg.resolve()
end

# ╔═╡ f5d2228c-b45e-4dce-99fb-668f6c50df30
@mdx """
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

# ╔═╡ 681ece7a-6800-4779-8118-28c3179cd43a
@mdx """
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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CondaPkg = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
MarkdownLiteral = "736d6165-7244-6769-4267-6b50796e6954"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"

[compat]
CondaPkg = "~0.2.4"
MarkdownLiteral = "~0.1.1"
PlutoUI = "~0.7.32"
PythonCall = "~0.5.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CodecBzip2]]
deps = ["Bzip2_jll", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "2e62a725210ce3c3c2e1a3080190e7ca491f18d7"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.7.2"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.CommonMark]]
deps = ["Crayons", "JSON", "URIs"]
git-tree-sha1 = "4aff51293dbdbd268df314827b7f409ea57f5b70"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.5"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.CondaPkg]]
deps = ["MicroMamba", "Pkg", "TOML"]
git-tree-sha1 = "b257b31acb6b41e95759bfd0bfd2729252209ea6"
uuid = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
version = "0.2.4"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "22df5b96feef82434b07327e2d3c770a9b21e023"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MarkdownLiteral]]
deps = ["CommonMark", "HypertextLiteral"]
git-tree-sha1 = "0d3fa2dd374934b62ee16a4721fe68c418b92899"
uuid = "736d6165-7244-6769-4267-6b50796e6954"
version = "0.1.1"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.MicroMamba]]
deps = ["CodecBzip2", "Downloads", "Scratch", "Tar"]
git-tree-sha1 = "5b0b23d102c53e83fb06d934b847ae242180e648"
uuid = "0b3b1443-0f03-428d-bdfb-f27f9c1191ea"
version = "0.1.3"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "92f91ba9e5941fc781fecf5494ac1da87bdac775"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "ae6145ca68947569058866e443df69587acc1806"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.32"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "2cf929d64681236a2e074ffafb8d568733d2e6af"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PythonCall]]
deps = ["CondaPkg", "Dates", "Libdl", "MacroTools", "Markdown", "Pkg", "Requires", "Serialization", "Tables", "UnsafePointers"]
git-tree-sha1 = "402a4d004a9e7c9656758f479f456205f7f6a3ce"
uuid = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
version = "0.5.1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnsafePointers]]
git-tree-sha1 = "c81331b3b2e60a982be57c046ec91f599ede674a"
uuid = "e17b2a0c-0bdf-430a-bd0c-3a23cae4ff39"
version = "1.0.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─f150770d-cbcd-4eda-a0a5-b759a21f5b9b
# ╠═acf2a1eb-e99f-442f-86de-5ebe5a2f601c
# ╟─30181f53-97e2-4d21-83d2-9394ef64aba2
# ╠═b1f1706d-da3b-48ae-bed1-4f5fb9e95726
# ╠═fa184670-718c-4e68-8e23-fe6e72a651c5
# ╠═53b89008-8ac0-4f20-955d-1d875530f121
# ╟─81a0a61d-0ff8-4610-915b-3f32a2742f5d
# ╠═15b2d0ab-bc1a-435a-9557-99537d117d91
# ╠═ffa07e57-7797-4000-8be4-a166895bcc54
# ╠═5e0919d8-0a02-46fc-92fe-0a5e40d1708c
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
# ╟─681ece7a-6800-4779-8118-28c3179cd43a
# ╠═655674de-56c1-4386-8fda-aa5c95b6271f
# ╠═9f68d823-64ac-48cc-b35e-b131a0bd5c50
# ╠═22b67cb0-6687-11ec-1868-bb216a9703f4
# ╟─782a5542-860b-4862-9bb9-dc4b1c2f48bb
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
