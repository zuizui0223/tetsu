using Circuitscape
using DelimitedFiles
if length(ARGS) != 2
    error("Expected configuration INI followed by output TSV path")
end
result = Circuitscape.compute(ARGS[1])
open(ARGS[2],"w") do stream
    writedlm(stream,result,'\t')
end
println("ACTUAL JULIA Circuitscape network pairwise matrix size ",size(result))
