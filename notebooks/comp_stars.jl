### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ 27c0f880-b775-458f-9092-a0d4d743d9ab
begin
	import Pkg
	Pkg.activate(Base.current_project())
	using PlutoUI # Handy for making widgets and such, like a table of contents
end

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

SMF files are essentially just text files with a list of target names and positions that are sent to the telescope to let the observatory know how to cut the mask for observations. They also include additonal information for the slits, but this will not be the focus for this notebook. Now this is all well and good, except for the fact that these SMF files tend to look something like this:
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
## 1. Read the SMF 🔎

We are first going to use [`CSV.jl`](https://github.com/JuliaData/CSV.jl) to read the SMF file into a `DataFrame` (provided by [`DataFrames.jl`](https://github.com/JuliaData/DataFrames.jl)):
"""

# ╔═╡ 89136b56-99f6-49bf-a823-6447fee208fc
df_SMF = CSV.read(download("https://github.com/icweaver/WASP-50b/raw/main/notebooks/data/raw_data/wasp50s.SMF"), DataFrame;
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
## 2. Query the star catalog 📖

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

	for the conversion, but pipes `|>` can help make multi-step workflows more transparent.
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
We generated a [`PyPandasDataFrame`](https://cjdoris.github.io/PythonCall.jl/stable/compat/#Tabular-data-and-Pandas) via `PythonCall.jl`'s `PyTable` function, which is basically a generic table type in Julia. We'd like to access the full feature set of `DataFrames.jl`, so we will construct one directly from `results`:
"""

# ╔═╡ 37d2bc2e-ee75-42a9-ba4f-df15c99caf48
_df_paper = DataFrame(
	Star = df.name,
	UCAC4_ID = results.UCAC4,
	# RAJ2000 = results.RAJ2000,
	# DecJ200 = results.DEJ2000,
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
	@transform :D = compute_D.(:B, :V, :J, :K; ΔBVₜ, ΔJKₜ)
end

# ╔═╡ 8fbc74d9-c993-4757-865f-b7ec3d84495c
latexify(df_paper, latex=false, fmt="%.2f")

# ╔═╡ c9b37015-057e-4838-b510-5762819b2462
md"""
!!! note
	We used the `aside` macro from `DataFramesMeta.jl` to compute intermediate values needed for ``D`` before adding it to our table, which was all done in a chain of actions with `@chain` from the [`Chain.jl`](https://github.com/jkrumbiegel/Chain.jl) package that is automatically included (exported) with `DataFramesMeta.jl`.

	We could also have used bangs (!) to mutate tables in-place and save memory, e.g.:

	```julia
	@chain df_paper begin
		@aside ΔBVₜ, ΔJKₜ = _.B[1] - _.V[1], _.J[1] - _.K[1]
		@transform! :D = compute_D.(:B, :V, :J, :K; ΔBVₜ, ΔJKₜ)
	end
	```

	but for clarity, we are just creating a new table each time in this notebook.
"""

# ╔═╡ 06d8b523-4a05-43ec-973b-8fe55fa613a1
md"""
Ok, now we can display `df_paper` in a ``LaTeX`` context with [`Latexify.jl`](https://github.com/korsbo/Latexify.jl):
"""

# ╔═╡ 4f3afae4-7b3c-42ed-81a0-936e1348700a
t = latexify(df_paper, env=:tabular, latex=false, fmt="%.2f") |> Text

# ╔═╡ 875ca00d-0c0c-403e-bbf1-0d2b70a95107
t

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
We now have a copy-and-pasteable table ✔
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

# ╔═╡ Cell order:
# ╟─80664df1-a135-44f9-b83f-41029a503c77
# ╟─46bee12a-81ca-43ea-932f-8cfa22ad1e9a
# ╟─4a079964-5ac0-47a0-9ff0-c294df62772e
# ╟─3dcf92a9-21c2-4cfd-9089-e6fb501c7caf
# ╟─3941e471-6378-4f6b-a117-81aa08eff0c7
# ╟─8fbc74d9-c993-4757-865f-b7ec3d84495c
# ╟─875ca00d-0c0c-403e-bbf1-0d2b70a95107
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
# ╠═27c0f880-b775-458f-9092-a0d4d743d9ab
# ╟─4188732d-1d88-4e98-8510-f74224641919
