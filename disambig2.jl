module Disambig

include("load.jl")

typealias AttributeDict Dict{Int, Set{Int}}

type AttributeMatrix
    LN::AttributeDict
    FN::AttributeDict
    NM::AttributeDict
    AS::AttributeDict
    CT::AttributeDict
    CL::AttributeDict
end

immutable Inventor
    lastname::Int
    firstname::Int
    name::Int
    assignee::Int
    city::Int
    class::Int
end

function getattribute(code, X, words)
    keep = [beginswith(w, code) for w in words]
    Xcode = X[:, find(keep)]
    Xcodet = Xcode'

    code_inventor = AttributeDict()
    inventor = Array(Int, X.m)

    for j in 1:Xcode.n
        code_inventor[j] = Set(find(Xcode[:, j]))
    end

    for i in 1:Xcode.m
        nz = find(Xcodet[:, i])
        if length(nz) != 1
            throw(BoundsError("Found $(length(x)) $code attributes for inventor $i"))
        end

        inventor[i] = nz[1]
    end

    return code_inventor, inventor
end

function prepattributes(X, words)
    LN, iLN = getattribute("LN", X, words)
    FN, iFN = getattribute("FN", X, words)
    NM, iNM = getattribute("NM", X, words)
    AS, iAS = getattribute("AS", X, words)
    CT, iCT = getattribute("CT", X, words)
    CL, iCL = getattribute("CL", X, words)

    amat = AttributeMatrix(LN, FN, NM, AS, CT, CL)

    ilist = Array(Inventor, X.m)
    for i in 1:X.m
        ilist[i] = Inventor(iLN[i], iFN[i], iNM[i], iAS[i], iCT[i], iCL[i])
    end

    amat, ilist
end

function lumpit(amat, inventor, step)
    C_lastname = amat.LN[inventor.lastname]
    C_firstname = amat.FN[inventor.firstname]
    C_name = amat.NM[inventor.name]
    C_assignee = amat.AS[inventor.assignee]
    C_city = amat.CT[inventor.city]
    C_class = amat.CL[inventor.class]

    lump_index_2 = intersect(step, intersect(C_name, union(C_assignee, C_city, C_class)))
    lump_index_1 = intersect(step, intersect(C_firstname, C_lastname, C_assignee, C_city, C_class))

    union(lump_index_1, lump_index_2)
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

    tic()

    amat, ilist = prepattributes(X, words)

    patno_1 = collect(1:endof(patno))
    step = Set(patno_1)

    inventor_id_new = copy(inventor_id)
    namey = copy(namex)

    for index in 1:length(patno)
        if !in(index, step)
            continue
        end

        inventor = ilist[index]
        lump_index = lumpit(amat, inventor, step)
        lump_initial = copy(lump_index)

        # find repeated patent No.s adopted

        for indexy in lump_initial
            inventor_ = ilist[indexy]
            lump_index_ = lumpit(amat, inventor_, step)
            union!(lump_index, lump_index_)
        end

        if length(lump_index) == 0
            setdiff!(step, [index])
            continue
        end

        for i in setdiff(lump_index, [index])
            inventor_id_new[i] = inventor_id_new[index]
            namey[i] = namey[index]
        end

        setdiff!(step, lump_index)
    end

    open(joinpath(directory, "_disambiguator2_output.tsv"), "w") do f
        for (k, iid) in zip(key, inventor_id_new)
            print(f, k, '\t', iid)
        end
    end

    toc()
end

if !isinteractive()
    main(ARGS[1])
end

end
