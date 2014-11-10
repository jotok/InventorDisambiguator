module Disambig

include("load.jl")

function main(directory=".")
    X = load_attribute_matrix(directory)
    attribute_metric = load_attribute_metric(directory)
    words = load_attribute_dictionary(directory)
    disambiguator_input = load_disambiguator_input(directory)

    patno = attribute_metric["patno"]
    inventor = attribute_metric["inventor"]
    assignee = attribute_metric["assignee"]
    class = attribute_metric["class"]
    lastname = attribute_metric["lastname"]
    inventor_id = attribute_metric["inventor_id"]
    key = disambiguator_input["key"]
    namex = disambiguator_input["namex"]

    XX = X # X is already a boolean array
    XXt = X'

    tic()

    # matrix and sub-matrix transpose to expedite memory access

    keep = map(w -> beginswith(w, "LN"), words)
    XXLN = XX[:, find(keep)]
    XXLNt = XXLN'

    keep = map(w -> beginswith(w, "FN"), words)
    XXFN = XX[:, find(keep)]
    XXFNt = XXFN'

    keep = map(w -> beginswith(w, "NM"), words)
    XXNM = XX[:, find(keep)]
    XXNMt = XXNM'

    keep = map(w -> beginswith(w, "AS"), words)
    XXAS = XX[:, find(keep)]
    XXASt = XXAS'

    keep = map(w -> beginswith(w, "CT"), words)
    XXCT = XX[:, find(keep)]
    XXCTt = XXCT'

    keep = map(w -> beginswith(w, "CL"), words)
    XXCL = XX[:, find(keep)]
    XXCLt = XXCL'

    patno_1 = collect(1:endof(patno))

    # step = ones(Int, length(patno))
    step = sparsevec(ones(Int, length(patno)))

    inventor_id_new = copy(inventor_id)
    namey = copy(namex)

    for index in 1:length(patno)
        if step[index] == 0
            continue
        end

        target = XXLNt[:, index]
        C_lastname = XXLN[:, find(target)]

        target = XXFNt[:, index]
        C_firstname = XXFN[:, find(target)]

        target = XXNMt[:, index]
        C_name = XXNM[:, find(target)]

        target = XXASt[:, index]
        C_assignee = XXAS[:, find(target)]

        target = XXCTt[:, index]
        C_city = XXCT[:, find(target)]

        target = XXCLt[:, index]
        C_class = XXCL[:, find(target)]

        lump_index_2 = find(step .* (C_name .* (C_assignee + C_city + C_class)))
        lump_index_1 = find(step .* (C_firstname .* C_lastname .* C_assignee .* C_city .* C_class))
        lump_patno_2 = patno[lump_index_2]
        lump_patno_1 = patno[lump_index_1]
        lump_index_1_ = Int[]

        # if exists in patno_1 and in patno_2, then take out of patno_1 (full match)
        for (ii, ix) in enumerate(lump_index_1)
            if !any([pn == lump_patno_1[ii] for pn in lump_patno_2])
                push!(lump_index_1_, ix)
            end
        end

        lump_index = union(lump_index_1_, lump_index_2)

        # find repeated patent No.s adopted

        lump_initial = copy(lump_index)
        lump_after = copy(lump_initial)

        for indexy in lump_initial
            target = XXLNt[:, indexy]
            C_lastname = XXLN[:, find(target)]

            target = XXFNt[:, indexy]
            C_firstname = XXFN[:, find(target)]

            target = XXNMt[:, indexy]
            C_name = XXNM[:, find(target)]

            target = XXASt[:, indexy]
            C_assignee = XXAS[:, find(target)]

            target = XXCTt[:, indexy]
            C_city = XXCT[:, find(target)]

            target = XXCLt[:, indexy]
            C_class = XXCL[:, find(target)]

            lump_index_2 = find(step .* (C_name .* (C_assignee + C_city + C_class)))
            lump_index_1 = find(step .* (C_firstname .* C_lastname .* C_assignee .* C_city .* C_class))
            lump_patno_2 = patno[lump_index_2]
            lump_patno_1 = patno[lump_index_1]
            lump_index_1_ = Int[]

            # if exists in patno_1 and in patno_2, then take out patno_1 (full match)
            for (ii, ix) in enumerate(lump_index_1)
                if !any([pn == lump_patno_1[ii] for pn in lump_patno_2])
                    push!(lump_index_1_, ix) 
                end
            end

            lump_index = union(lump_index_1_, lump_index_2)
            lump_after_2 = copy(lump_index)

            lump_after = union(lump_after, lump_after_2)

        end

        lump = union(lump_initial, lump_after)

        if length(lump) == 0
            step[index] = 0
            continue
        end

        for i in length(lump):-1:1
            step[lump[i]] = 0
            inventor_id_new[lump[i]] = inventor_id_new[index]
            namey[lump[i]] = namey[index]
        end

    end

    open(joinpath(directory, "_disambiguator_output.tsv"), "w") do f
        for (k, iid) in zip(key, inventor_id_new)
            @printf(f, "%s\t%s\n", k, iid)
        end
    end

    toc()
end

if !isinteractive()
    main(ARGS[1])
end

end
