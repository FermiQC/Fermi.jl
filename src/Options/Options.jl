# default options go here
# if user specifies option they will be overwritten
CurrentOptions = Dict{String,Any}(
                                  "reference" => "rhf",
                                  "quiet" => "true",
                                  "mp2_type" => "conv",
                                  "cc_type" => "conv",
                                  "e_conv" => 10,
                                  "d_conv" => 8,
                                  "basis" => "sto-3g"
                                 )