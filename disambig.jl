module Disambig

include("load.jl")

typealias SparseMatrix SparseMatrixCSC{Int, Int}

function make_lump_index_1(step, C_firstname, C_lastname, C_assignee, C_city, C_class)

    nnzS = max(nnz(step), nnz(C_firstname), nnz(C_lastname), nnz(C_assignee), nnz(C_city), nnz(C_class))
    rowvalS = Array(Int, nnzS)
    ptrS = 1

    step_rowval = step.rowval
    fn_rowval = C_firstname.rowval
    ln_rowval = C_lastname.rowval
    as_rowval = C_assignee.rowval
    ct_rowval = C_city.rowval
    cl_rowval = C_class.rowval

    if isempty(fn_rowval) || isempty(ln_rowval) || isempty(as_rowval) || isempty(ct_rowval) || isempty(cl_rowval) || isempty(step_rowval)
        return Array(Int, 0)
    end

    step_i = fn_i = ln_i = as_i = ct_i = cl_i = 1
    r = max(fn_rowval[1], ln_rowval[1], as_rowval[1], ct_rowval[1], cl_rowval[1])

    while fn_i <= endof(fn_rowval) && ln_i <= endof(ln_rowval) && as_i <= endof(as_rowval) && ct_i <= endof(ct_rowval) && cl_i <= endof(cl_rowval) && step_i <= endof(step_rowval)

        fn_r = fn_rowval[fn_i]
        while fn_r < r && fn_i < endof(fn_rowval)
            fn_i += 1
            fn_r = fn_rowval[fn_i]
        end

        ln_r = ln_rowval[ln_i]
        while ln_r < r && ln_i < endof(ln_rowval)
            ln_i += 1
            ln_r = ln_rowval[ln_i]
        end

        as_r = as_rowval[as_i]
        if as_r < r && as_i < endof(as_rowval)
            as_i += 1
            as_r = as_rowval[as_i]
        end

        ct_r = ct_rowval[ct_i]
        if ct_r < r && ct_i < endof(ct_rowval)
            ct_i += 1
            ct_r = ct_rowval[ct_i]
            continue
        end

        cl_r = cl_rowval[cl_i]
        if cl_r < r && cl_i < endof(cl_rowval)
            cl_i += 1
            cl_r = cl_rowval[cl_i]
            continue
        end

        step_r = step_rowval[step_i]
        if step_r < r && step_r < endof(step_rowval)
            step_i += 1
            step_r = step_rowval[step_i]
            continue
        end

        if fn_r == ln_r == as_r == ct_r == cl_r == step_r
            rowvalS[ptrS] = fn_r
            ptrS += 1
        end

        fn_i += 1
        ln_i += 1
        as_i += 1
        ct_i += 1
        cl_i += 1
        step_i += 1
    end

    deleteat!(rowvalS, ptrS:length(rowvalS))
    return rowvalS
end

function make_lump_index_2(step, C_name::SparseMatrix, C_assignee::SparseMatrix, C_city::SparseMatrix, C_class::SparseMatrix)

    nnzS = nnz(step) + nnz(C_name)
    rowvalS = Array(Int, nnzS)
    ptrS = 1

    step_rowval = step.rowval
    nm_rowval = C_name.rowval
    as_rowval = C_assignee.rowval
    ct_rowval = C_city.rowval
    cl_rowval = C_class.rowval

    if isempty(nm_rowval) || isempty(step_rowval)
        return Array(Int, 0)
    end

    step_i = nm_i = as_i = ct_i = cl_i = 1
    while step_i <= endof(step_rowval) && nm_i <= endof(nm_rowval)
        step_r = step_rowval[step_i]
        nm_r = nm_rowval[nm_i]

        if step_r < nm_r
            step_i += 1
            continue
        end

        if step_r > nm_r
            nm_i += 1
            continue
        end

        # step_r == nm_r

        while as_i <= endof(as_rowval) && as_rowval[as_i] < step_r
            as_i += 1
        end

        while ct_i <= endof(ct_rowval) && ct_rowval[ct_i] < step_r
            ct_i += 1
        end

        while cl_i <= endof(cl_rowval) && cl_rowval[cl_i] < step_r
            cl_i += 1
        end

        if as_i > endof(as_rowval) || ct_i > endof(ct_rowval) || cl_i > endof(cl_rowval)
            break
        end

        if as_rowval[as_i] == step_r || ct_rowval[ct_i] == step_r || cl_rowval[cl_i] == step_r
            rowvalS[ptrS] = step_r
            ptrS += 1
        end

        step_i += 1
        nm_i += 1
    end

    deleteat!(rowvalS, ptrS:length(rowvalS))
    return rowvalS
end

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

        lump_index_2 = make_lump_index_2(step, C_name, C_assignee, C_city, C_class)
        lump_index_1 = make_lump_index_1(step, C_firstname, C_lastname, C_assignee, C_city, C_class)
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

            lump_index_2 = make_lump_index_2(step, C_name, C_assignee, C_city, C_class)
            lump_index_1 = make_lump_index_1(step, C_firstname, C_lastname, C_assignee, C_city, C_class)
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
