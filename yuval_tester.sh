#!/bin/bash




# INSTRUCTIONS: Read, else get a bozo message.
#	1. You shall place this shell script within the directory that contains all of the files
#	That you need to assign.
#	
#	2. You shouldn't run this from Nova/Beast/rack-mad/any tau-related server.
#	If you still wish to do so, you would need to replace the path of $output_file to somewhere
#	that isn't write-protected (tau-related servers have /tmp/ write-protected)
#
#	3. Be sure to remove any build/dist/egg files from the working directory. Could
#	allegdly lead to undefined behaviors of the test script.
#
#	4. The most important and crucial instruction: you shall supply as the first and only
#	cmd-line argument, the relative (or absolute) path to the directory of tests that yuval created.
#	You'll get a bozo message if you don't do so.
#
#	5. Yuval might increase the amount of test output files, so beware. You would probably
#	need to change some arguments here and there.



# global variables
testers_path=$1
output_file="/tmp/output.txt"



# auxiliary functions
function verdict_diff() {
	# the first argument shall be the length (in characters) of the result of diff
	if [[ $1 -eq 0 ]]; then
		echo -e '\033[1;32mSUCCESS\e[0m' # print out a 'success' message
	else
		echo -e "\e[1;31mFAILED\e[0m" # print out a 'failed' message
	fi
}



function individual_test() {
	# the first argument shall be the interface being tested: c/py
	# the second argument shall be the goal being tested
	# the third argument shall be the input file being used

	# the function returns the output of the 'diff' operation upon the output testing file and the actual output

	if [[ "${1}" == "py" ]]; then # if we are testing the python interface
		python3 spkmeans.py 0 $2 $testers_path/$3 &> $output_file
	elif [[ "${1}" == "c" ]]; then # if we are testing the C interface
		./spkmeans $2 $testers_path/$3 > $output_file
	else
		echo "Individual test function failed: Invalid interface"
		return -1
	fi

	diff_result=$(diff $output_file $testers_path/outputs/$1/$2/$3)
	echo -e "DIFF RESULT FOR: ${1}: ${2}: ${3}:\n${diff_result}\n\n" >> test_transcript_$1.txt

	return ${#diff_result}
}



function test_goal() {
	# the first argument shall be the interface being tested c/py
	# the second argument shall be the goal being tested

	if [[ "${2}" == "jacobi" ]]; then
		for i in {0..19}; do
			individual_test $1 $2 jacobi_$i.txt
			diff_res_len=$?
			echo -n "${1^^}: ${2^^}: ${testers_path}/jacobi_${i}.txt: "
			verdict_diff $diff_res_len
		done
	else
		for i in {0..9}; do
			individual_test $1 $2 spk_$i.txt
			diff_res_len=$?
			echo -n "${1^^}: ${2^^}: ${testers_path}/spk_${i}.txt: "
			verdict_diff $diff_res_len
		done
	fi
}



function test_interface() {
	# the first argument shall be the interface being tested
	interface=$1

	test_goal $interface wam
	test_goal $interface ddg
	test_goal $interface lnorm
	test_goal $interface jacobi

	if [[ $interface == "py" ]]; then
		test_goal $interface spk
	fi
}




function run_comprehensive_test() {
	# testing the C interface
	rm test_transcript_c.txt &> /dev/null
	touch test_transcript_c.txt
	echo -e "Testing the interface for: \e[4;33m\e[1;33mC\e[0m"
	echo -e "\n\e[4;34m\e[1;34mRESULTS\e[0m"
	bash comp.sh &> /dev/null # compiling
	test_interface c

	# testing the CPython interface
	rm test_transcript_py.txt &> /dev/null
	touch test_transcript_py.txt
	echo -e "\n\nTesting the interface for: \e[4;33m\e[1;33mPython\e[0m"
	echo -e "\n\e[4;34m\e[1;34mRESULTS\e[0m"
	python3 setup.py build_ext --inplace &> /dev/null # building
	test_interface py

	# Summary 
	echo -e "\n\n\e[1;31mNOTICE:\e[0m More detailed results are in: \e[4;34m\e[1;34mtest_transcript_c.txt\e[0m and \e[4;34m\e[1;34mtest_transcript_py.txt\e[0m"
}


if [[ $# -ne 1 ]]; then
	echo -e "
 _               ____      _  _____ ___ ___            ____   ___ ________  
| |        _    |  _ \    / \|_   _|_ _/ _ \     _    | __ ) / _ \__  / _ \ 
| |      _| |_  | |_) |  / _ \ | |  | | | | |  _| |_  |  _ \| | | |/ / | | |
| |___  |_   _| |  _ <  / ___ \| |  | | |_| | |_   _| | |_) | |_| / /| |_| |
|_____|   |_|   |_| \_\/_/   \_\_| |___\___/    |_|   |____/ \___/____\___/ 
                                                                            
          ____ ___  ____  _____           ____  _  _____ _     _     
   _     / ___/ _ \|  _ \| ____|    _    / ___|| |/ /_ _| |   | |    
 _| |_  | |  | | | | |_) |  _|    _| |_  \___ \| ' / | || |   | |    
|_   _| | |__| |_| |  __/| |___  |_   _|  ___) | . \ | || |___| |___ 
  |_|    \____\___/|_|   |_____|   |_|   |____/|_|\_\___|_____|_____|
                                                                     
 ___ ____ ____  _   _ _____ 
|_ _/ ___/ ___|| | | | ____|
 | |\___ \___ \| | | |  _|  
 | | ___) |__) | |_| | |___ 
|___|____/____/ \___/|_____|"
else
	run_comprehensive_test $1
fi
