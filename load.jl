function load_attribute_matrix(directory=".")
    cols = readdlm(joinpath(directory, "_attribute_matrix"), ',', Int)
    sparse(cols[:, 1], cols[:, 2], cols[:, 3])
end

function load_attribute_metric(directory=".")
    cols = readdlm(joinpath(directory, "_attribute_metric"), '\t', String, quotes=false)

    result = Dict{String, Vector{String}}()
    result["patno"] = cols[:, 1]
    result["inventor"] = cols[:, 3]
    result["assignee"] = cols[:, 4]
    result["class"] = cols[:, 5]
    result["lastname"] = cols[:, 6]
    result["inventor_id"] = cols[:, 7]

    result
end

function load_attribute_dictionary(directory=".")
    cols = readdlm(joinpath(directory, "_attribute_dictionary"), '\t', quotes=false)
    cols[:, 2]
end

function load_disambiguator_input(directory=".")
    key = String[]
    namex = String[]

    open(joinpath(directory, "_disambiguator_input.csv")) do f
        for line in eachline(f)
            fields = [strip(s) for s in split(line, "\t")]
            push!(key, fields[1])
            push!(namex, join(fields[2:4], " "))
        end
    end

    map!(t -> strip(replace(t, r"\s+", " ")), namex)
    ["key" => key, "namex" => namex]
end
