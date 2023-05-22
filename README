

Intel fortran (with BLAS, MKL interfacing with LAPACK) can be installed here:
# https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html


At the beginning of a session, you must summon Intel fortran with:

  On a desktop:
# source /opt/intel/oneapi/setvars.sh 

  On an example HPC server:
# . ~/.bashrc
# module purge
# module load toolchain/intel/2018.5.274
# module load numlib/ScaLAPACK/2.1.0-gompi-2021b-fb




Make a NN with 2 hidden layers with architecture X-Y-Z-1 (e.g. X=528,Y=20,Z=29):
# echo "528 20 29" | awk '{print 2; print $1; print $2; print $3; Nw=$1*($2)+$2*($3+1)+$3*(1+1)+1; print Nw; for(i=1;i<=$1+1;i++){printf "%.8f  %.8f\n", 0.0, 1.0}; i1=1;0;for(i=1;i<=$2;i++){i2+=$1+1;printf "%6d %6d    %6d %6d\n", 0+1, $1+1, i1, i2;i1+=$1+1};for(i=1;i<=$3;i++){i2+=$2+1;printf "%6d %6d    %6d %6d\n", $1+1+1,$1+1+$2+1, i1, i2;i1+=$2+1};printf "%6d %6d    %6d %6d\n", $1+1+$2+1+1, $1+1+$2+1+$3+1, i1, i1+$3; for(i=1;i<=Nw;i++){printf "%.8f\n", 0.50}}'




(The command below assumes the PIP/FI fortran file is SIMPLIFIED (e.g. s2(3)*s3(3) -> s5(3))
To convert an exiting FI fortran file to a dFI fortran file (the derivatives), do something like (for orders up to 9):
# awk '{flag=1} /&/ {flag=0;line=line" "$0} {if(flag==1){print line" "$0;line=""}}' obj/BrCH5.f90 | grep 'p([0-9]*)=' | while read line; do p=$(echo "$line" | grep -oh 'p([0-9]*)' | grep -oh '[0-9]*'); seq 1 21 | while read i; do echo "$line" | sed -e 's/[ &]*//g' -e 's/^p.*=//' | sed 's/+/ /g' | xargs -n1 | while read term; do if echo "$term" | grep -q 'r('"$i"')'; then echo "$term" | sed 's/r('"$i"')/1/' | sed -e 's/\*1//' -e 's/1\*//'; fi; seq 2 9 | while read n; do if echo "$term" | grep -q 's'"$n"'('"$i"')'; then let "nm1 = $n - 1"; echo "$term" | sed 's/s'"$n"'('"$i"')/'"$n"'*s'"$nm1"'('"$i"')/'; fi; done | sed 's/s1/r/'; done | xargs | sed 's/ /+/g' | awk '{if (NF==0) {print "0.0d0"} else {print}}' | sed 's/^/dp('"$p"','"$i"')=/' | sed 's/=1$/=1.0d0/'; done; done > obj/BrCH5-derivatives.f90
(The fortran header, footer, and variable declarations must still be added after doing this command)

