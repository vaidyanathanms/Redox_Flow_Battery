molecule='TMED'
charge=('0' 'p1' 'p2' 'm1' 'm2')
curr_dir=$PWD

for qval in "${charge[@]}"; do
    cp gauss_var.sh $molecule/DME/q_$qval/
    cd $molecule/DME/q_$qval
    echo $PWD
    rm optim_nodiffuse_$molecule.com *.chk Default.Route Gau* std* *.log gauss.sh
    wait
    sed -e "s/6-31+G\**/6-31G(d,p)/g" optim_$molecule.com > optim_nodiffuse_$molecule.com
    wait
    sed -e "s/molname/$molecule/g" -e "s/charge/q$qval/g" gauss_var.sh > gauss.sh
    wait
    rm gauss_var.sh
    sbatch gauss.sh
    cd $curr_dir
done
