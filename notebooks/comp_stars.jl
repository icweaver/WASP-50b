### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ e2d82196-7945-11ec-2d31-ad911cabbfba
using CSV, DataFrames

# ╔═╡ 183ac353-4850-4cc6-be57-ba849fd50739
using DataFramesMeta

# ╔═╡ 403b0807-419c-4804-86b8-aca28751a1a6
using AstroAngles

# ╔═╡ 1be7d092-b655-48b3-87e1-8ef006aedd28
using PythonCall

# ╔═╡ b258ea27-bd58-41c1-b5e6-b9c6ccdf69d9
using CondaPkg

# ╔═╡ 8293e156-e76e-41bc-a3cd-07877eb23663
using LinearAlgebra # Defines `norm` for us

# ╔═╡ b3df4b34-9287-4de0-bb26-e717ddad067d
using Latexify

# ╔═╡ 27c0f880-b775-458f-9092-a0d4d743d9ab
using PlutoUI # Handy for making widgets and such, like a table of contents

# ╔═╡ 80664df1-a135-44f9-b83f-41029a503c77
md"""
# Observations at LCO

In this notebook, we are going to take an *ugly* SMF file and convert it into a gorgeous ``\LaTeX`` table, complete with star coordinates pulled down from VizieR via [astroquery](https://astroquery.readthedocs.io/en/latest/) ✨

$(TableOfContents())
"""

# ╔═╡ 46bee12a-81ca-43ea-932f-8cfa22ad1e9a
md"""
## Introduction 🤚

Las Campanas Observatory ([LCO](http://www.lco.cl/)) is home to several world-class telescopes that are used in a wide range of astronomy fields. In particular, The 6.5 meter twin [Magellan telescopes](https://obs.carnegiescience.edu/Magellan) provide essential facilities for observing exoplanet atmopsheres with their [IMACS and LDSS3C spectrographs](http://www.lco.cl/instruments/).

**[More background intro for general public about multi-object slit spectroscopy and need for SMF files]**

SMF files are essentially just text files with a list of target names and positions that are sent to the telescope to let the observatory know how to cut the mask for observations. They also include additonal information for the slits, but this will not be the focus for this notebook. Now this is all well and good, except for that these SMF files tend to look something like this:
"""

# ╔═╡ 4a079964-5ac0-47a0-9ff0-c294df62772e
"NAME wasp50s
OBSERVER Cookie Monster
TITLE wasp50 science mask
MADE 2013-08-09  MaskGen  2.27.92 
TELESCOPE 71161.1 0.80906 3 Magellan
INSTRUMENT IMACS_sc 
DISPERSER  IMACS_grism_300  1  17.50000
POSITION 0 02:55:00.000 -10:57:00.00 2000.0  0.0 0.0 0.0
!Constellation (Eri) Eridanus
!# No rotator warnings above horizon.
!.OC wasp50s 02:55:00.000 -10:57:00.00 2000.0 0.0 0.0 -136.15 OFF  02:55:11.3 -10:41:29 2000.0 02:54:41.9 -11:12:17 2000.0 
WAVELENGTH 6749.3
TEMPERATURE 12.0
EPOCH 2013.70211
WLIMIT 4000.1 7900.1
DREF  0
HANGLE    0
DATE  56548
SLIT WASP50    02:54:45.137 -10:53:53.03   75.570   64.529  3.450  3.795  3.795    0.00 # Vmag= 11.721
SLIT c06       02:54:49.392 -10:51:54.94   53.950  105.319  3.450  3.795  3.795    0.00 # Vmag= 11.585
SLIT c13       02:54:23.502 -11:02:28.19  185.851 -113.558  3.450  3.795  3.795    0.00 # Vmag= 13.075
SLIT c15       02:54:32.141 -10:54:36.94  141.733   49.389  3.450  3.795  3.795    0.00 # Vmag= 13.252
SLIT c18       02:54:56.565 -11:02:27.36   17.461 -113.018  3.450  3.795  3.795    0.00 # Vmag= 13.568
SLIT c20       02:54:14.107 -11:03:06.82  233.973 -127.107  3.450  3.795  3.795    0.00 # Vmag= 13.662
SLIT c21       02:54:52.547 -10:51:26.96   37.908  114.984  3.450  3.795  3.795    0.00 # Vmag= 13.668
SLIT c23       02:55:16.318 -10:53:50.69  -82.969   65.341  3.450  3.795  3.795    0.00 # Vmag= 13.708
SLIT c28       02:55:35.977 -10:56:20.97 -183.127   13.446  3.450  3.795  3.795    0.00 # Vmag= 13.885
HOLE ref015    02:54:53.542 -11:08:45.55   32.898 -244.204  1.725 1  0.862  0.862    0.00 # Vmag= 14.897
HOLE ref031    02:55:46.106 -10:51:22.87 -235.193  116.663  1.725 1  0.862  0.862    0.00 # Vmag= 15.430
HOLE ref036    02:54:05.105 -10:56:02.52  280.104   19.807  1.725 1  0.862  0.862    0.00 # Vmag= 15.490
HOLE ref042    02:54:40.376 -11:01:26.73   99.779  -92.124  1.725 1  0.862  0.862    0.00 # Vmag= 15.714
HOLE ref044    02:54:06.394 -10:50:53.83  273.783  126.843  1.725 1  0.862  0.862    0.00 # Vmag= 15.728
HOLE ref046    02:54:53.749 -10:55:28.22   31.764   31.665  1.725 1  0.862  0.862    0.00 # Vmag= 15.762
HOLE ref053    02:54:56.291 -10:45:54.96   18.908  230.088  1.725 1  0.862  0.862    0.00 # Vmag= 15.865
HOLE ref054    02:55:19.594 -10:49:31.69  -99.775  154.931  1.725 1  0.862  0.862    0.00 # Vmag= 15.894
HOLE ref058    02:55:44.049 -11:00:33.94 -224.424  -74.096  1.725 1  0.862  0.862    0.00 # Vmag= 15.997
HOLE ref068    02:55:08.761 -11:01:26.28  -44.527  -91.920  1.725 1  0.862  0.862    0.00 # Vmag= 16.266
HOLE ref071    02:54:25.846 -11:09:07.73  174.283 -252.376  1.725 1  0.862  0.862    0.00 # Vmag= 16.322
HOLE ref072    02:54:42.547 -10:52:41.78   88.772   89.150  1.725 1  0.862  0.862    0.00 # Vmag= 16.351
HOLE ref076    02:54:18.818 -10:46:32.58  210.387  217.464  1.725 1  0.862  0.862    0.00 # Vmag= 16.407
HOLE ref080    02:55:46.477 -10:59:38.07 -236.851  -54.783  1.725 1  0.862  0.862    0.00 # Vmag= 16.567
HOLE ref085    02:54:19.177 -10:58:22.60  207.882  -28.623  1.725 1  0.862  0.862    0.00 # Vmag= 17.101
<and then some random long string>
" |> Text

# ╔═╡ 3dcf92a9-21c2-4cfd-9089-e6fb501c7caf
md"""
!!! note
	More info about what each column means can be found [here](https://code.obs.carnegiescience.edu/maskgen/v214/smdfile.txt/view).
"""

# ╔═╡ 3941e471-6378-4f6b-a117-81aa08eff0c7
md"""
We need a way to ingest this highly unstructured file and extract some useful information from it. For this notebook, the final product will be this summary report that we can paste into a ``\LaTeX`` document:
"""

# ╔═╡ af553292-ae3e-4966-9746-79de4b86a7e1
md"""
## Game plan 🏈
To generate the above, we are going to do the following:

1. Extract only the rows from the SMF file that refer to the slits (the rest just refer to other objects that we don't care about for this analysis)

2. Use the RA and Dec coordinates to find the corresponding star in each row using `astroquery`

3. Save each star's ID and a few other important quantities in a final report table

Alright, let's get started.
"""

# ╔═╡ 769c4b03-4ac0-4a25-9e7d-21825261c3dd
md"""
## 1. Data ingestion 📖

We are first going to use [`CSV.jl`](https://github.com/JuliaData/CSV.jl) to read the SMF file into a `DataFrame` (provided by [`DataFrames.jl`](https://github.com/JuliaData/DataFrames.jl)):
"""

# ╔═╡ 89136b56-99f6-49bf-a823-6447fee208fc
df_SMF = CSV.read(download("https://github.com/icweaver/Pluto_sample_notebooks/raw/main/data/wasp50s.SMF"), DataFrame;
	header = false, # Wish there was one
	skipto = 19, # This is the row where the good stuff starts
	footerskip = 1, # Skip that random string at the bottom
	delim = ' ', # Spaces, spaces everywhere
	ignorerepeated = true, # Treat repeated delimiters as a single one
	comment = "HOLE", # Ignore all lines starting with this
	types = String, # Treat every column as a `String`
)

# ╔═╡ 16d610c6-6b7e-46a2-bd54-487392248302
md"""
Next let's just select the columns that we will need for the report. We use the `@select` macro from the [`DataFramesMeta.jl`](https://github.com/JuliaData/DataFramesMeta.jl) package to accomplish this in a convenient way: 
"""

# ╔═╡ 8ceb1ae4-eecd-4e05-94cf-6e84749d6727
df = @select df_SMF begin
	:name = :Column2
	:_RAJ2000 = @. string(hms2deg(:Column3)) * "d" # astroquery needs columns
	:_DEJ2000 = @. string(dms2deg(:Column4)) * "d" # in this format
end

# ╔═╡ ed944c35-09a6-4b69-8ccc-20e238fa488d
md"""
!!! tip
	There is A LOT more that we can do with these kinds of macros. A nice overview from the creator of `DataFrames.jl` can be found [here](https://bkamins.github.io/julialang/2021/11/19/dfm.html).
"""

# ╔═╡ 07d51156-9b82-4b9b-b6a6-a8bf09966880
md"""
We also used `hms2deg` and `dms2deg` from [`AstroAngles.jl`](https://github.com/JuliaAstro/AstroAngles.jl) to convert the RA and Dec into a format that can automatically be read by `astroquery`, which we will do next.
"""

# ╔═╡ 21e927a2-ae8b-4623-ba17-e98ca3c04f60
md"""
## 2. Pull down some data 📚

Here comes the fun part. We are going to feed this table into the Python package `astroquery`, which will query an all-sky star catalog for us and return the corresponding table of results. Let's start by setting up a minimal Python environment and download `astroquery` into it: 
"""

# ╔═╡ 75a654c0-f15d-4cf0-b8bd-8f7a10de0e59
md"""
We use [`PythonCall.jl`](https://github.com/cjdoris/PythonCall.jl) to install a light-weight Python environment (using [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)) into a temporary folder via [`CondaPkg.jl`](https://github.com/cjdoris/CondaPkg.jl). `PythonCall.jl` will also enable us to seamlessly interact with Python objects from Julia. Finally, we download `astroquery` and a few other Python packages that will be needed for compatibility reasons later:
"""

# ╔═╡ 5337ca99-20c5-4418-8fff-76191d718522
CondaPkg.add.(("astroquery", "astropy", "pandas")); CondaPkg.resolve()

# ╔═╡ 137bb020-1aae-4eb3-9ad7-976163f9f685
md"""
Now that that's all set up, let's import some things!
"""

# ╔═╡ 5926658d-bcb1-477a-bebc-9414eaa3813c
@py begin
	import astroquery.vizier as vz
	import astropy.table as ap
end

# ╔═╡ 9a4f6266-f761-4eca-8eba-dd8f611a56bf
md"""
For SOME reason, the only flavor of table that `astroquery` accepts is an astropy `Table`, so we are going to do the following roundabout conversions to make the 🐍 happy:

```
DataFrame (julia) ---> DataFrame (pandas) ---> Table (astropy)
```

We accomplish the first step with `PythonCall.jl`'s `pytable` function.
"""

# ╔═╡ f9d21956-347b-48fc-abe2-04d733c92668
results_astropy = vz.Vizier.query_region(
	df |> pytable |> ap.Table.from_pandas,
	radius = "1 arcsec",
	catalog = "UCAC4",
	frame = "icrs", # I think this is default
)[0]

# ╔═╡ 9c710f1f-4029-4678-8ed8-654ffe107a59
md"""
!!! note
	We also could have done:
	
	```julia
	ap.Table.from_pandas(pytable(df))
	```

	for the conversion, but pipes `|>` can make things look clearer sometimes.
"""

# ╔═╡ 905127e7-16ae-4d8d-8445-7197cc5dee2a
md"""
Alright, we have a table of results now! These are the list of objects (from the [UCAC4 catalog](https://vizier.u-strasbg.fr/viz-bin/VizieR?-source=I/322)) that, row-by-row from `df`, are within 1 arcsecond of the RA and Decs specified. We set this radius to be small enough so that only one target from the catalog (ideally the star we are searching for) matches the coords given.
"""

# ╔═╡ 1f2afbc8-5e04-4f18-b047-90775443d572
md"""
## 3. Generate the report ✏️

Now that we have our list of targets, let's extract a few key columns and compute a color metric, for example. First let's convert the returned astropy `Table` back to a native Julia `DataFrame` to make things a bit easier to work with:
"""

# ╔═╡ 15eef4f3-53ac-45f0-953c-0e46e14f3e46
results = results_astropy.to_pandas() |> PyTable

# ╔═╡ ea756fff-853a-41b8-902b-a495cee1fd43
md"""
We generated a `PyPandasDataFrame` via `PythonCall.jl`'s `PyTable` function, which is basically a generic table type in Julia. We'd like to access all of the cool features from `DataFrames.jl`, so we will construct one directly from `results`:
"""

# ╔═╡ 37d2bc2e-ee75-42a9-ba4f-df15c99caf48
_df_paper = DataFrame(
	Star = df.name,
	UCAC4_ID = results.UCAC4,
	RAJ2000 = results.RAJ2000,
	DecJ200 = results.DEJ2000,
	B = results.Bmag,
	V = results.Vmag,
	J = results.Jmag,
	K = results.Kmag,
)

# ╔═╡ d6f7f9eb-5776-4419-98e0-f18cd4cee12f
md"""
Now that we have our color magnitudes for each target, we compute a color-space metric ``D``:

```math
D \equiv \sqrt{\left[(B-V)_{c}-(B-V)_{t}\right]^{2}+\left[(J-K)_{c}-(J-K)_{t}\right]^{2}}
```

relative to our main star, "WASP-50":
"""

# ╔═╡ 3a2891dc-36bf-4e1f-913a-4ba53b47eac6
"""
	compute_D(B, V, J, K, [ΔBVₜ=0.0, ΔJKₜ=0.0])

* `B, V, J, K`: comparison star magnitudes:

```math
	B_c, V_c, J_c, K_c
```

* `ΔBVₜ`, `ΔJKₜ`: target star delta mags (defaults to zero for each):

```math
(B-V)_{t}, (J-K)_{t}
```
"""
compute_D(B, V, J, K; ΔBVₜ=0.0, ΔJKₜ=0.0) = norm(((B-V)-ΔBVₜ, (J-K)-ΔJKₜ))

# ╔═╡ 96035871-6bcb-4703-9507-3365f92a93d9
md"""
and add it to our table:
"""

# ╔═╡ f3326ab8-21da-403b-b290-7694cefb63fe
df_paper = @chain _df_paper begin
	@aside ΔBVₜ, ΔJKₜ = _.B[1] - _.V[1], _.J[1] - _.K[1]
	@transform :D = compute_D.(:B, :V, :J, :K; ΔBVₜ=ΔBVₜ, ΔJKₜ=ΔJKₜ)
end

# ╔═╡ 8fbc74d9-c993-4757-865f-b7ec3d84495c
t = latexify(df_paper, latex=false, fmt="%.2f")

# ╔═╡ c9b37015-057e-4838-b510-5762819b2462
md"""
!!! note
	We used the `aside` macro from `DataFramesMeta.jl` to compute intermediate values needed for ``D`` before adding it to our table, which was all done in a chain of actions with `@chain` from the [`Chain.jl`](https://github.com/jkrumbiegel/Chain.jl) package that is automatically included (exported) with `DataFramesMeta.jl`.

	We could also have used bangs (!) to mutate tables in-place and save memory, e.g.:

	```julia
	@chain df_paper begin
		@aside ΔBVₜ, ΔJKₜ = _.B[1] - _.V[1], _.J[1] - _.K[1]
		@transform! :D = compute_D.(:B, :V, :J, :K; ΔBVₜ=ΔBVₜ, ΔJKₜ=ΔJKₜ)
	end
	```

	but for clarity, we are just creating a new table each time in this notebook.
"""

# ╔═╡ 06d8b523-4a05-43ec-973b-8fe55fa613a1
md"""
Ok, now we can display `df_paper` in a ``LaTeX`` context with [`Latexify.jl`](https://github.com/korsbo/Latexify.jl):
"""

# ╔═╡ 4f3afae4-7b3c-42ed-81a0-936e1348700a
latexify(df_paper, env=:tabular, latex=false, fmt="%.2f") |> Text

# ╔═╡ 52100a8e-3adf-4337-9985-645a10d9e8cd
md"""
!!! note
	We used `latex=false` to remove the dollar signs around each number, because this will be used in a [deluxetable](https://journals.aas.org/aastexguide/#preamble_deluxetable) environment that can handle automatic math mode columns.

	To display the nicely rendered ``\LaTeX`` table at the beginning of the document, we just used the default Markdown backend:

	```julia
	latexify(df_paper, latex=false, fmt="%.2f")
	```
"""

# ╔═╡ 34d1417b-2d56-4de0-800e-e65d5f9817fb
md"""
We now have a table that we can copy-and-paste into a paper that is definitely not being ignored right now.
"""

# ╔═╡ 4188732d-1d88-4e98-8510-f74224641919
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
AstroAngles = "5c4adb95-c1fc-4c53-b4ea-2a94080c53d2"
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CondaPkg = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DataFramesMeta = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"

[compat]
AstroAngles = "~0.1.3"
CSV = "~0.10.1"
CondaPkg = "~0.2.4"
DataFrames = "~1.3.1"
DataFramesMeta = "~0.10.0"
Latexify = "~0.15.9"
PlutoUI = "~0.7.30"
PythonCall = "~0.5.0"
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

[[deps.AstroAngles]]
git-tree-sha1 = "41621fa5ed5f7614b75eea8e0b3cfd967b284c87"
uuid = "5c4adb95-c1fc-4c53-b4ea-2a94080c53d2"
version = "0.1.3"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "fbee070c56e0096dac13067eca8181ec148468e1"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.1"

[[deps.Chain]]
git-tree-sha1 = "339237319ef4712e6e5df7758d0bccddf5c237d9"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.4.10"

[[deps.CodecBzip2]]
deps = ["Bzip2_jll", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "2e62a725210ce3c3c2e1a3080190e7ca491f18d7"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.7.2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.CondaPkg]]
deps = ["MicroMamba", "Pkg", "TOML"]
git-tree-sha1 = "b257b31acb6b41e95759bfd0bfd2729252209ea6"
uuid = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
version = "0.2.4"

[[deps.Crayons]]
git-tree-sha1 = "b618084b49e78985ffa8422f32b9838e397b9fc2"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.0"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "cfdfef912b7f93e4b848e80b9befdf9e331bc05a"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.1"

[[deps.DataFramesMeta]]
deps = ["Chain", "DataFrames", "MacroTools", "OrderedCollections", "Reexport"]
git-tree-sha1 = "ab4768d2cc6ab000cd0cec78e8e1ea6b03c7c3e2"
uuid = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
version = "0.10.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "04d13bfa8ef11720c24e4d840c0033d145537df7"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.17"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

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

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "8d70835a3759cdd75881426fced1508bb7b7e1b6"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

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

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

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

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.MicroMamba]]
deps = ["CodecBzip2", "Downloads", "Scratch", "Tar"]
git-tree-sha1 = "5b0b23d102c53e83fb06d934b847ae242180e648"
uuid = "0b3b1443-0f03-428d-bdfb-f27f9c1191ea"
version = "0.1.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

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
git-tree-sha1 = "5c0eb9099596090bb3215260ceca687b888a1575"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.30"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "db3a23166af8aebf4db5ef87ac5b00d36eb771e2"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "2cf929d64681236a2e074ffafb8d568733d2e6af"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.3"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PythonCall]]
deps = ["CondaPkg", "Dates", "Libdl", "MacroTools", "Markdown", "Pkg", "Requires", "Serialization", "Tables", "UnsafePointers"]
git-tree-sha1 = "4069e80d13c3b33a2b5680a1464baf264c36e5a0"
uuid = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
version = "0.5.0"

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

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "15dfe6b103c2a993be24404124b8791a09460983"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.11"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

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

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnsafePointers]]
git-tree-sha1 = "c81331b3b2e60a982be57c046ec91f599ede674a"
uuid = "e17b2a0c-0bdf-430a-bd0c-3a23cae4ff39"
version = "1.0.0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "c69f9da3ff2f4f02e811c3323c22e5dfcb584cfa"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.1"

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
# ╟─80664df1-a135-44f9-b83f-41029a503c77
# ╟─46bee12a-81ca-43ea-932f-8cfa22ad1e9a
# ╟─4a079964-5ac0-47a0-9ff0-c294df62772e
# ╟─3dcf92a9-21c2-4cfd-9089-e6fb501c7caf
# ╟─3941e471-6378-4f6b-a117-81aa08eff0c7
# ╟─8fbc74d9-c993-4757-865f-b7ec3d84495c
# ╟─af553292-ae3e-4966-9746-79de4b86a7e1
# ╟─769c4b03-4ac0-4a25-9e7d-21825261c3dd
# ╠═e2d82196-7945-11ec-2d31-ad911cabbfba
# ╠═89136b56-99f6-49bf-a823-6447fee208fc
# ╟─16d610c6-6b7e-46a2-bd54-487392248302
# ╠═183ac353-4850-4cc6-be57-ba849fd50739
# ╠═403b0807-419c-4804-86b8-aca28751a1a6
# ╠═8ceb1ae4-eecd-4e05-94cf-6e84749d6727
# ╟─ed944c35-09a6-4b69-8ccc-20e238fa488d
# ╟─07d51156-9b82-4b9b-b6a6-a8bf09966880
# ╟─21e927a2-ae8b-4623-ba17-e98ca3c04f60
# ╠═1be7d092-b655-48b3-87e1-8ef006aedd28
# ╟─75a654c0-f15d-4cf0-b8bd-8f7a10de0e59
# ╠═b258ea27-bd58-41c1-b5e6-b9c6ccdf69d9
# ╠═5337ca99-20c5-4418-8fff-76191d718522
# ╟─137bb020-1aae-4eb3-9ad7-976163f9f685
# ╠═5926658d-bcb1-477a-bebc-9414eaa3813c
# ╟─9a4f6266-f761-4eca-8eba-dd8f611a56bf
# ╠═f9d21956-347b-48fc-abe2-04d733c92668
# ╟─9c710f1f-4029-4678-8ed8-654ffe107a59
# ╟─905127e7-16ae-4d8d-8445-7197cc5dee2a
# ╟─1f2afbc8-5e04-4f18-b047-90775443d572
# ╠═15eef4f3-53ac-45f0-953c-0e46e14f3e46
# ╟─ea756fff-853a-41b8-902b-a495cee1fd43
# ╠═37d2bc2e-ee75-42a9-ba4f-df15c99caf48
# ╟─d6f7f9eb-5776-4419-98e0-f18cd4cee12f
# ╠═8293e156-e76e-41bc-a3cd-07877eb23663
# ╠═3a2891dc-36bf-4e1f-913a-4ba53b47eac6
# ╟─96035871-6bcb-4703-9507-3365f92a93d9
# ╠═f3326ab8-21da-403b-b290-7694cefb63fe
# ╟─c9b37015-057e-4838-b510-5762819b2462
# ╟─06d8b523-4a05-43ec-973b-8fe55fa613a1
# ╠═b3df4b34-9287-4de0-bb26-e717ddad067d
# ╠═4f3afae4-7b3c-42ed-81a0-936e1348700a
# ╟─52100a8e-3adf-4337-9985-645a10d9e8cd
# ╟─34d1417b-2d56-4de0-800e-e65d5f9817fb
# ╟─27c0f880-b775-458f-9092-a0d4d743d9ab
# ╟─4188732d-1d88-4e98-8510-f74224641919
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
