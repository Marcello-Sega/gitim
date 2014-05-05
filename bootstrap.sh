find ./ -name "CMakeCache.txt" -or -name "CMakeFiles" -exec rm -rf {} \;
cmake .
