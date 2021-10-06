for (( r = 1; r < 4; r++))
do
    julia create_inputs.jl
    for (( i = 1; i < 23; i++))
    do
        cd S$i
        nohup psi4
        cd ..
    done
    mkdir R$r
    mv S* R$r/
done

