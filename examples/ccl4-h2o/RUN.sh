#!/bin/bash

GROMPP=`which grompp`

GDENS1=`which gitim`
GDENS2=`which ../../gitim`

if [ -z  $GDENS1 ] ; then
	if [ -z $GDENS2  ] ; then 
		echo "No working version of gitim  was found neither in the path, not in the parent folder, cannot continue."
                exit
        else
		GDENS=$GDENS2
        fi
else 
        GDENS=$GDENS1
fi

if [ -z $GROMPP ] ; then 
	echo "could not find grompp in the path, cannot continue." 
	exit
fi

echo "gitim needs a .tpr file, let's create it now...."
grompp -f grompp.mdp -p topol.top -c ccl4-h2o.gro  -maxwarn 10  -o ccl4-h2o.tpr> grompp.log 2>&1 
echo "0 4 5" > masscom.dat
echo "Using $GDENS ..."
if [ -a ccl4-h2o.tpr ] ; then 
        name=dens_H2O_atomic_m ; echo "Intrinsic mass density profile w.r.t. water -> $name.xvg"
	echo "SOL SOL CCl4"  | $GDENS  -intrinsic -dens mass -f ccl4-h2o.gro -ng 3 -center -sl 200 -s ccl4-h2o.tpr  -o $name.xvg  > $name.log 2>&1 
        name=dens_CCl4_atomic_m ; echo "Intrinsic mass density profile w.r.t. ccl4-> $name.xvg"
	echo "CCl4 SOL CCl4" | $GDENS  -intrinsic -dens mass -f ccl4-h2o.gro -ng 3 -center -sl 200 -s ccl4-h2o.tpr  -o $name.xvg  > $name.log   2>&1 
        name=dens_H2O_atomic_n ; echo "Intrinsic number density profile w.r.t. water -> $name.xvg"
	echo "SOL SOL CCl4"  | $GDENS  -intrinsic -dens number -f ccl4-h2o.gro -ng 3 -center -sl 200 -s ccl4-h2o.tpr    -o $name.xvg  > $name.log   2>&1 
        name=dens_CCl4_atomic_n ; echo "Intrinsic number density profile w.r.t. ccl4 -> $name.xvg"
	echo "CCl4 SOL CCl4" | $GDENS  -intrinsic -dens number -f ccl4-h2o.gro -ng 3 -center -sl 200 -s ccl4-h2o.tpr    -o $name.xvg  > $name.log  2>&1 
        name=dens_H2O_atomicMC_m ; echo "Intrinsic mass density profile w.r.t. water with Monte Carlo normalization-> $name.xvg"
	echo "SOL SOL CCl4"  | $GDENS  -intrinsic -MCnorm -dens mass -f ccl4-h2o.gro -ng 3 -center -sl 200 -s ccl4-h2o.tpr  -o $name.xvg  > $name.log 2>&1 
        name=dens_H2O_molecularMC_m ; echo "Intrinsic molecular mass density profile w.r.t. water with Monte Carlo normalization-> $name.xvg"
	echo "SOL SOL CCl4"  | $GDENS  -intrinsic -MCnorm -dens mass -com -f ccl4-h2o.gro -ng 3 -center -sl 200 -s ccl4-h2o.tpr  -o $name.xvg  > $name.log 2>&1 
else 
	echo "Some errors probably occurred running grompp. Please have a look at grompp.log"

fi  



