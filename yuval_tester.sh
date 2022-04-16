#!/bin/bash



testers_path=$1
output_file="/tmp/output.txt"


function individual_test() {
	# the first argument shall be the interface being tested: c/py
	# the second argument shall be the goal being tested
	# the third argument shall be the input file being used 

	# the function returns the output of the 'diff' operation upon the output testing file and the actual output

	if [[ "${1}" == "py" ]]; then # if we are testing the python interface
		python3 spkmeans.py 0 
	elif [[ "${1}" == "c" ]]; then # if we are testing the C interface
		./spkmeans $2 $testers_path/$3 > $output_file		
	else
		return -1
	fi

	diff_result=$(diff $output_file $testers_path/outputs/$1/$2/$3)
	return diff_result
}



function test_goal() {
	# the first argument shall be the interface being tested c/py
	# the second argument shall be the goal being tested

	
}



function test_interface() {
	# the first argument shall be the interface being tested
	interface=$1

	
}






echo "Testing the interface for: \e[4;33m\e[1;33mC\e[0m"
echo -e "\n\e[4;34m\e[1;34mRESULTS\e[0m"
bash comp.sh &> /dev/null # compiling



echo "Testing the interface for: \e[4;33m\e[1;33mPython\e[0m"
echo -e "\n\e[4;34m\e[1;34mRESULTS\e[0m"
python3 setup.py build_ext --inplace &> /dev/null # building



echo -e "\e[1;31mFAILED\e[0m" # print out a 'failed' message
echo -e '\033[1;32mSUCCESS\e[0m' # print out a 'success' message
