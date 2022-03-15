## WASP-50b

Collection of [Pluto.jl](https://github.com/fonsp/Pluto.jl) notebooks used to create figures for this project. To
interact with Python, we use the [PythonCall.jl](https://github.com/cjdoris/PythonCall.jl) package:

* [Quick tutorial](https://icweaver.github.io/WASP-50b/html/fun_with_python.html)
* [Example](https://icweaver.github.io/WASP-50b/html/comp_stars.html)

HTML pages for all figure notebooks are available below and include directions for installing and running each one.

1. [Raw data](https://icweaver.github.io/WASP-50b/html/01_raw_data.html)
1. [Reduced data -- IMACS](https://icweaver.github.io/WASP-50b/html/02_reduced_data_IMACS.html)
1. [Reduced data -- LDSS3](https://icweaver.github.io/WASP-50b/html/03_reduced_data_LDSS3.html)
1. [Detrended white-light curves](https://icweaver.github.io/WASP-50b/html/04_detrended_wlcs.html)
1. [Detrended binned light curves](https://icweaver.github.io/WASP-50b/html/05_detrended_blcs.html)
1. [Photometric monitoring](https://icweaver.github.io/WASP-50b/html/06_photometric_monitoring.html)
1. [Transmission spectra](https://icweaver.github.io/WASP-50b/html/07_transmission_spectra.html)
1. [Retrievals](https://icweaver.github.io/WASP-50b/html/08_retrievals.html)
1. [HGHJ Population](https://icweaver.github.io/WASP-50b/html/09_pop.html)

The data repository can be accessed [here](https://app.box.com/s/fwohk8q6dp9wgufa3gv14b492xain11t).

### General organization and navigation
The top of each notebook points to a link to download the specific dataset for each notebook. Each figure in the paper
is given its own section. Within each section, the overal goal is first stated and then the figure is displayed. The
analysis code to produce each figure is then shown directly underneath. The plotting-specific code is hidden for
clarity, but can be revealed by clicking on the "eyeball" icon next to each figure in the downloaded notebook. All
dependent cells in the notebook will run and save the figures in the same directory. The path to each figure will also
be displayed in the notebook surrounded by a blue box.
