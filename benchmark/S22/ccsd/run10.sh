cd fermi
for i in 1 2 3 4 5 6 7 8 9 10
do
    nohup julia input.jl
    mkdir R$i
    mv output.dat R$i
    mv nohup.dat R$i

    cd ../psi4
    nohup psi4 -n 12
    mkdir R$i
    mv output.dat R$i
    mv timer.dat R$i

    cd ../fermi
done
