# default options go here
# if user specifies option they will be overwritten
CurrentOptions = Dict{String,Any}(
                                  "molstring" => """
                                  O        1.2091536548      1.7664118189     -0.0171613972
                                  H        2.1984800075      1.7977100627      0.0121161719
                                  H        0.9197881882      2.4580185570      0.6297938832
                                  """,
                                  "basis" => "sto-3g",
                                  "charge" => 0,
                                  "multiplicity" => 1,
                                  "reference" => "rhf",
                                  "quiet" => "true",
                                  "mp2_type" => "conv",
                                  "precision" => "double",
                                  "cc_type" => "conv",
                                  "e_conv" => 10,
                                  "d_conv" => 8,
                                  "cc_max_iter" => 50,
                                  "cc_max_rms" => 10^-10,
                                  "cc_e_conv" => 10^-10,
                                  "diis" => false,
                                  "num_frozen" => 0,
                                  "cas_frozen" => 0,
                                  "cas_active" => "all",
                                  "cas_nroot" => 1,
                                  "min_matrix_elem" => 10^-9
                                 )
