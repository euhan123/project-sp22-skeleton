"""Solves an instance.

Modify this file to implement your own solvers.

For usage, run `python3 solve.py --help`.
"""

import argparse
from audioop import add
from locale import currency
from pathlib import Path
from typing import Callable, Dict

from point import Point
from instance import Instance
from solution import Solution
from file_wrappers import StdinFileWrapper, StdoutFileWrapper


def solve_naive(instance: Instance) -> Solution:
    return Solution(
        instance=instance,
        towers=instance.cities,
    )

def solve_greedy(instance: Instance) -> Solution:
    length = instance.D
    coverR = instance.R_s
    penaltyR = instance.R_p
    cities = instance.cities
    N = instance.N
    towers = []
    towerDict = {}
    citiesCovered = {}

    def findPenalty(point):
        penalty = 170
        for i in towers:
            if (Point.distance_sq(point, i) <= penaltyR * penaltyR):
                penalty *= 2

        return penalty

    def findCities(point):
        count = 0
        for i in cities:
            string = str(i.x) + " " + str(i.y)
            if (Point.distance_sq(point, i) <= coverR * coverR) and string not in citiesCovered:
                count += 1

        return count
    
    def checkCities(point):
        for i in cities:
            if (Point.distance_sq(point, i) <= coverR * coverR):
                xcoord = i.x
                ycoord = i.y
                string = str(xcoord) + " "  + str(ycoord)
                citiesCovered[string] = 0


    greedySolution = Solution(instance = instance, towers = towers)
    while (not greedySolution.valid()):
        minimumPenalty = [[], 0]
        maxCity = 0
        minX = 0
        minY = 0
        for i in range(length):
            for j in range(length):
                string = str(i) + " " + str(j)
                if (string not in towerDict):
                    tower = Point(i, j)
                    currCities = findCities(tower)                        
                    currPenalty = findPenalty(tower)
                    if currCities > maxCity:
                        minX = i
                        minY = j
                        maxCity = currCities
                        minimumPenalty[1] = currPenalty
                    elif currCities == maxCity:
                        if currPenalty < minimumPenalty[1]:
                            minX = i
                            minY = j
                            minimumPenalty[1] = currPenalty
        
        string = str(minX) + " " + str(minY)
        towerDict[string] = 0
        addTower = Point(minX, minY)
        checkCities(addTower)
        towers.append(addTower)
        greedySolution = Solution(instance = instance, towers = towers)

    return greedySolution

def solve_greedy2(instance: Instance) -> Solution:
    length = instance.D
    coverR = instance.R_s
    penaltyR = instance.R_p
    cities = instance.cities
    N = instance.N
    towers = []
    citiesCovered = []

    def findPenalty(point):
        penalty = 170
        for i in towers:
            if (Point.distance_sq(point, i) <= penaltyR * penaltyR and point != i):
                penalty *= 2
        return penalty

    def findCities(point):
        count = 0
        for i in cities:
            if (Point.distance_sq(point, i) <= coverR * coverR) and i not in citiesCovered:
                count += 1
        return count
    
    def checkCities(citiesCovered):
        citiesCovered = []
        for i in towers:
            for j in cities:
                if(Point.distance_sq(i, j) <= coverR*coverR):
                    citiesCovered.append(j)
                
    def sum_penalty():
        sum = 0
        for i in towers:
            sum += findPenalty(i)
        return sum
    
    def better_tower(point):
        total = sum_penalty()
        maximum = findPenalty(point)
        tow = point
        for i in towers:
            if(maximum < findPenalty(i)):
                maximum = findPenalty(i)
                tow = i
        curr = 0
        for i in range(length):
            for j in range(length):
                tower = Point(i, j)
                if(tower not in towers):
                    temp = findCities(tower)
                    towers.remove(tow)
                    towers.append(tower)
                    total_temp = sum_penalty()
                    if((total > total_temp and temp >= curr) or (total >= total_temp and temp > curr) and tower != tow ):
                        curr = temp
                        tow = tower
                        total = total_temp
                    else:
                        towers.remove(tower)
                        towers.append(tow)        
        return tow
    greedySolution = Solution(instance = instance, towers = towers)
    
    while (not greedySolution.valid()):
        minimumPenalty = 0
        maxCity = 0
        minX = 0
        minY = 0
        
        for i in range(length):
            for j in range(length):
                tower = Point(i, j)
                if(tower not in towers):
                    currCities = findCities(tower)           
                    currPenalty = findPenalty(tower)
                    if currCities > maxCity:
                        minX = i
                        minY = j
                        maxCity = currCities
                        minimumPenalty = currPenalty
                    elif currCities == maxCity:
                        if currPenalty < minimumPenalty:
                            minX = i
                            minY = j
                            minimumPenalty = currPenalty
        addTower = Point(minX, minY)
        towers.append(addTower)
        temp = towers
        temp_tower = better_tower(addTower)
        while(temp != towers):
            temp = towers
            checkCities(citiesCovered)
            temp_tower  = better_tower(temp_tower)
        checkCities(citiesCovered)
    greedySolution = Solution(instance = instance, towers = towers)    
 
    return greedySolution


SOLVERS: Dict[str, Callable[[Instance], Solution]] = {
    "naive": solve_naive,
    "greedy": solve_greedy
}


# You shouldn't need to modify anything below this line.
def infile(args):
    if args.input == "-":
        return StdinFileWrapper()

    return Path(args.input).open("r")


def outfile(args):
    if args.output == "-":
        return StdoutFileWrapper()

    return Path(args.output).open("w")


def main(args):
    with infile(args) as f:
        instance = Instance.parse(f.readlines())
        solver = SOLVERS[args.solver]
        solution = solver(instance)
        assert solution.valid()
        with outfile(args) as g:
            print("# Penalty: ", solution.penalty(), file=g)
            solution.serialize(g)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Solve a problem instance.")
    parser.add_argument("input", type=str, help="The input instance file to "
                        "read an instance from. Use - for stdin.")
    parser.add_argument("--solver", required=True, type=str,
                        help="The solver type.", choices=SOLVERS.keys())
    parser.add_argument("output", type=str,
                        help="The output file. Use - for stdout.",
                        default="-")
    main(parser.parse_args())
