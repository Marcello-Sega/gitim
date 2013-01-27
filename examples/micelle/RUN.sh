#!/bin/bash

GROMPP=`which grompp`

GDENS1=`which gitim`
GDENS2=`which ../../gitim`

if [ -z $GDENS1 ] ; then
	if [  -z $GDENS2  ] ; then 
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
grompp -f grompp.mdp -p micelle.top -c micelle.gro  -o micelle.tpr -maxwarn 20 > grompp.log 2>&1 
echo "0 4 5" > masscom.dat
echo "Using $GDENS ..."
if [ -a micelle.tpr ] ; then 
        name=dens_DPC_atomic_m ; echo "Intrinsic mass density profile w.r.t. DPC -> $name.xvg"
	echo "DPC SOL DPC"  | $GDENS  -intrinsic -dens mass -f micelle.gro -ng 3 -center -sl 200 -s micelle.tpr -MCnorm -alpha 0.2 -geometry sphere -o $name.xvg  > $name.log 2>&1 
else 
	echo "Some errors probably occurred running grompp. Please have a look at grompp.log"

fi  



