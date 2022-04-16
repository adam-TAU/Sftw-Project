import sys
from typing import List
import numpy as np
import pandas as pd
import spkmeans




def main(K: int, goal: str, infile: str) -> None:
    
    # In case that we desire a normalized spectral clustering:
	if goal == "spk":
		with open(infile) as tmpFile:
			num_data = len(tmpFile.readlines()) 
    
        # We want K to be less or equal to num_data, and not equal to 1
		assert_valid_input(K <= num_data)
		assert_valid_input(K != 1)

        # Fetch the matrix of points produced from the eigen vectors of the normalized graph laplacian matrix of the given vectors
		T_points = spkmeans.goal(K, goal, infile)

        # Initialize centroids picked by the Kmeans++ algorithm
		initial_centroids_indices = initialize_centroids(len(T_points[0]), T_points)

        # Perform the kmeans clustering algorith on the
		spkmeans.kmeans_fit(T_points, len(T_points), len(T_points[0]), initial_centroids_indices, len(initial_centroids_indices)) 

	else: # In any other case that isn't a normalized spectral clustering - just perform the desired operation corresponding to the "goal" parameter
		spkmeans.goal(K, goal, infile)



def initialize_centroids(K: int, datapoints: List[List[float]]) -> List[List[float]]:
    """
    Given a list of datapoints and an amount of centroids we wish to output 'K',
    weightedly pick the initial centroids for the kmeans process we'll perform later.
    This function will return their indices inside the given set of datapoints.
    """
    np.random.seed(0)
    index_range = list(range(len(datapoints)))
    centroids_indices = list()
    centroids_indices.append(int(np.random.choice(index_range)))

    for i in range(1, K):
        distance_list = [min([calc_distance(dp, datapoints[prev_obs]) for prev_obs in centroids_indices]) for dp in datapoints]
        sum_dist = sum(distance_list)
        probability_list = [dist/sum_dist for dist in distance_list]
        centroids_indices.append(int(np.random.choice(index_range, p=probability_list)))

    return centroids_indices


def calc_distance(u1: List[float], u2: List[float]) -> float:
    """
    Given two same length vectors, calculate the euclidean norm of their distance,
    and return it as the value
    """
    assert(len(u1) == len(u2))
    return sum([(u1[ind] - u2[ind]) ** 2 for ind in range(len(u1))])


def assert_valid_input(cond: bool):
    """
    If the provided condition is not satisfied, exits with the message 'Invalid Input!'
    """
    if not cond:
        print("Invalid Input!")
        sys.exit(1)


def assert_generic(cond: bool):
    """
    If the provided condition is not satisfied, exits with the message 'An Error Has Occurred'
    """
    if not cond:
        print("An Error Has Occurred")
        sys.exit(1)



def check_positive_numstr(string: str):
    """
    Checks if the input string represents a valid positive integer,
    and returns both the integer and whether it is valid.
    The value of the integer is 0 if it could not be parsed correctly.
    """
    if not string.isdecimal():
        return 0, False
    x=int(string)
    return x, (x > 0)



# Only run the program if it's actually called from the command line.
# Argument validation also happens here.
if __name__ == "__main__":
	num_args=len(sys.argv)
	assert_valid_input(num_args == 4)

    # We want K to be both numeric and positive.
	K, valid=check_positive_numstr(sys.argv[1])
	assert_valid_input(K == 0 or valid)
	
	# verify the correctness of our goal name
	goal = sys.argv[2]
	assert_valid_input(goal in ['wam', 'ddg', 'lnorm', 'jacobi', 'spk'])

	# The document specified that filenames must end with .txt or .csv, so we verify this here.
	infile = sys.argv[3]
	assert_valid_input( infile.endswith((".txt", ".csv")) )
	   
	# call the main mechanism
	main(K, goal, infile)




