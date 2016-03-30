#!/bin/bash
if [ -e "$TRAVIS_BUILD_DIR/build/CMakeFiles/CMakeOutput.log" ];  then 
	cat $TRAVIS_BUILD_DIR/build/CMakeFiles/CMakeOutput.log; 
else 
	echo "No CMakeOutput.log found..."; 
fi
if [ -e "$TRAVIS_BUILD_DIR/build/CMakeFiles/CMakeError.log" ];  then 
	cat $TRAVIS_BUILD_DIR/build/CMakeFiles/CMakeError.log; 
else 
	echo "No CmakeError.log found..."; 
fi

