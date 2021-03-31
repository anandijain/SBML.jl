# this is used to download the data from biomd
using CSV, DataFrames, JSON3

function curl_sbml_metadata()
    offsets = 0:100:2200
    urls = "https://www.ebi.ac.uk/biomodels/search?query=sbml&offset=" .* string.(offsets) .* "&numResults=100&format=json"
    for i in 1:length(urls) 
        run(`curl $(urls[i]) -o "data/sbml_$(i).json"`)
    end
end

function jsonfn_to_df(fn)
    json = read(fn, String);
    json = JSON3.read(json)
    DataFrame(jsontable(json.models))
end

function biomd_sbml_metadata()
    curl_sbml_metadata()
    fns = readdir("data"; join=true)
    dfs = jsonfn_to_df.(fns)
    df = vcat(dfs...)
    CSV.write("sbml_biomodels.csv", df)
    df
end

function sbml_zip_urls(df)
    base = "https://www.ebi.ac.uk/biomodels/search/download?models="
    ids = df.id
    N = 100 
    chunks = [ids[i:i+99] for i in 1:100:2200]
    append!(chunks, [ids[2201:end]])
    qs = join.(chunks, ",") 
    base .* qs  
end

"given the metadata dataframe from `biomd_sbml_metadata()`"
function curl_sbml_zips(df)
    urls = sbml_zip_urls(df)
    for i in 1:length(urls)
        @show i
        run(`curl -X GET "$(urls[i])" -H "accept: application/zip" -o data/yo_$i.zip`)
    end
end

# TODO ZipFile.jl function for extracting

# TODO SBML_TEST_SUITE grabber
