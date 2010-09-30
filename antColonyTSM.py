import random
import math

from Tkinter import *

class Matrix(object):
    def __init__(self, rows, columns, default=0):
        self.m = []
        for i in range(rows):
            self.m.append([default for j in range(columns)])

    def __getitem__(self, index):
        return self.m[index]

def paintTour(cities, tour):
    scale = 4
    root = Tk()
    root.title('Minmal Tour')

    canvas = Canvas(root, width=100*scale, height=100*scale, bg = 'white')
    canvas.pack()
    Button(root, text='Quit', command=root.quit).pack()

    #paint cities
    for city in cities:
        canvas.create_oval(city[0]*scale - 2,city[1]*scale - 2,city[0]*scale + 2,city[1]*scale + 2, width=2)

    #paint tour
    tourLength = len(tour)
    for point in range(1,tourLength):
        canvas.create_line(cities[tour[point]][0]*scale,cities[tour[point]][1]*scale,cities[tour[point-1]][0]*scale,cities[tour[point-1]][1]*scale,width=2)
    canvas.create_line(cities[tour[tourLength-1]][0]*scale,cities[tour[tourLength-1]][1]*scale,cities[tour[0]][0]*scale,cities[tour[0]][1]*scale,width=2)
    root.mainloop()   

def calculateTourLength(tour, distances):
    numberofCities = len(tour)
    tourLength = 0
    for i in range(1, numberofCities):
        tourLength = tourLength + distances[tour[i-1]][tour[i]]
    tourLength = tourLength + distances[tour[0]][tour[numberofCities-1]]
    return tourLength

def adjustPheromoneForAntTours(tau, tours, tourLength):
    numberOfAnts = len(tours);
    numberOfCities = len(tours[0])
    for a in range(numberOfAnts):
        for i in range(1,numberOfCities):
            tau[tours[a][i-1]][tours[a][i]] = tau[tours[a][i-1]][tours[a][i]] + 1/tourLength[a]
            tau[tours[a][i]][tours[a][i-1]] = tau[tours[a][i]][tours[a][i-1]] + 1/tourLength[a]
        tau[tours[a][0]][tours[a][numberOfCities-1]] = tau[tours[a][0]][tours[a][numberOfCities-1]] + 1/tourLength[a]
        tau[tours[a][numberOfCities-1]][tours[a][0]] = tau[tours[a][numberOfCities-1]][tours[a][0]] + 1/tourLength[a]
    return tau

def evaporatePheromone(tau, numberOfCities):
    for i in range(numberOfCities):
        for j in range(numberOfCities):        
            tau[i][j] = tau[i][j] * (1-rho)
    return tau

def recalculateWeightedTau(tau, numberOfCities):
    for i in range(numberOfCities):
        for j in range(numberOfCities):        
            weightedTau[i][j] =  pow(tau[i][j],alpha)    
    return weightedTau

def constructAntTour(numberOfCities, weightedTau, weightedDistances):
    # initialize tour
    visited = {0:True}
    tour = {0:0}
    strength = {}
    for i in range(1, numberOfCities):
        visited[i] = False

    # construct tour
    for i in range(1,numberOfCities):

        # calculate probability for each possible next city
        addedProb = 0
        for j in range(numberOfCities):
            if not visited[j]: 
                addedProb = addedProb + (weightedTau[tour[i-1]][j] / weightedDistances[tour[i-1]][j])
            strength[j] = addedProb

        # choose next city in tour
        randProb = addedProb * random.random()
        for j in range(numberOfCities-1):
            if (strength[j+1] >= randProb) and (strength[j] <= randProb):
               # next city found
               tour[i] = j+1;
               visited[j+1] = True
               break;

    return tour


cities = [
	[37, 52], [49, 49], [52, 64], [20, 26], [40, 30], [21, 47],
	[17, 63], [31, 62], [52, 33], [51, 21], [42, 41], [31, 32],
	[ 5, 25], [12, 42], [36, 16], [52, 41], [27, 23], [17, 33],
	[13, 13], [57, 58], [62, 42], [42, 57], [16, 57], [ 8, 52],
	[ 7, 38], [27, 68], [30, 48], [43, 67], [58, 48], [58, 27],
	[37, 69], [38, 46], [46, 10], [61, 33], [62, 63], [63, 69],
	[32, 22], [45, 35], [59, 15], [ 5,  6], [10, 17], [21, 10],
	[ 5, 64], [30, 15], [39, 10], [32, 39], [25, 32], [25, 55], 
	[48, 28], [56, 37], [30, 40]
]


numberOfCities = len(cities)
numberOfAnts = 500
numberofGenerations = 50
alpha = 1
beta = 2
rho = 0.5
minTourLength = 1000000
distances = Matrix(numberOfCities,numberOfCities)
weightedDistances = Matrix(numberOfCities,numberOfCities)
tau = Matrix(numberOfCities,numberOfCities)
weightedTau = Matrix(numberOfCities,numberOfCities)

#calculate distancess and initialize pheromone
for i in range(numberOfCities):
    for j in range(numberOfCities):
        tau[i][j] = 0.1 * random.random()
        weightedTau[i][j] = pow(tau[i][j],alpha)
        if i <> j:
            distances[i][j] = math.sqrt(pow(cities[i][0] - cities[j][0],2) + pow(cities[i][1] - cities[j][1],2))
            weightedDistances[i][j] = pow(distances[i][j],beta) 
        
for g in range(numberofGenerations):
    tours = {}
    tourLength = {}

    for a in range(numberOfAnts):
        tours[a] = constructAntTour(numberOfCities, weightedTau, weightedDistances)
        tourLength[a] = calculateTourLength(tours[a], distances)
        if tourLength[a] < minTourLength:
            minTourLength = tourLength[a]
            minTour = tours[a]
            print 'Generation:%d minimalTourLength:%4.2f' % (g, minTourLength)

    tau = adjustPheromoneForAntTours(tau, tours, tourLength)
    tau = evaporatePheromone(tau, numberOfCities)
    weightedTau = recalculateWeightedTau(tau, numberOfCities)

print "Finished"
print 'minimalTourLength:%4.2f' % (minTourLength)
print minTour
paintTour(cities,minTour)






