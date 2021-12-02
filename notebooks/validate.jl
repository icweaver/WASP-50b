using Glob

function get_species(dirpath)
    tokens_dir = split(dirpath, "_")
    return tokens_dir[findfirst(x -> occursin("fitR0", x), tokens_dir)+1:end]
end

function create_dict_opts(dirpath)
    dict_opts = Dict()
    for line in eachline(open("$(dirpath)/opts.py", "r"))
        if (length(line) ≥ 1) && (line[1] != '#')
            tokens = strip.(split(line, "="))
            if tokens[1] ∈ flags_opts
                if tokens[1] == "molecules"
                    token = filter(x -> any(isletter, x), split(tokens[2], "\""))
                elseif tokens[2] == "False"
                    token = false
                else
                    token = true
                end
                dict_opts[tokens[1]] = token
            end
        end
    end
    return dict_opts
end

function check(dict_opts, name_opts, val_dir)
    if dict_opts[name_opts] == val_dir
        return true
    else
        return false
    end
end

function print_results(dict_opts, name_opts, val_dir)
    if dict_opts[name_opts] != val_dir
        #println("$(name_opts) passes and is set to $(val_dir)")
    #else
        println("$(name_opts) fails: dir=$(val_dir) but opts=$(dict_opts[name_opts])")
    end
end

const flags_opts = ["heterogeneity", "clouds", "hazes", "fit_R0", "molecules"]

for dirpath in glob("bak_with_offs/all_WASP50/*")
    het_dir, clouds_dir, haze_dir, fitR0_dir = occursin.(
    ("_Het", "_Clouds", "_Haze", "_fitR0"), dirpath)

    species_dir = get_species(dirpath)

    dict_opts = create_dict_opts(dirpath)

    if !all((
        check(dict_opts, "heterogeneity", het_dir),
        check(dict_opts, "clouds", clouds_dir),
        check(dict_opts, "hazes", haze_dir),
        check(dict_opts, "fit_R0", fitR0_dir),
        check(dict_opts, "molecules", species_dir),
    ))
        printstyled("failed: $(dirpath)\n", color=:red)
        print_results(dict_opts, "heterogeneity", het_dir),
        print_results(dict_opts, "clouds", clouds_dir),
        print_results(dict_opts, "hazes", haze_dir),
        print_results(dict_opts, "fit_R0", fitR0_dir),
        print_results(dict_opts, "molecules", species_dir),
        println()
    end
end
