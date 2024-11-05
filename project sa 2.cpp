#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <ctime>

struct City {
    int id;
    double x, y;
};

double euclideanDistance(const City &a, const City &b) {
    return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

double calculateTotalDistance(const std::vector<City> &cities, const std::vector<int> &tour) {
    double totalDistance = 0.0;
    for (size_t i = 0; i < tour.size() - 1; i++) {
        totalDistance += euclideanDistance(cities[tour[i]], cities[tour[i + 1]]);
    }
    totalDistance += euclideanDistance(cities[tour.back()], cities[tour[0]]);  // Return to the starting point
    return totalDistance;
}

std::vector<int> generateInitialTour(int numCities) {
    std::vector<int> tour(numCities);
    std::iota(tour.begin(), tour.end(), 0);  // Fill with 0, 1, 2, ..., numCities - 1
    std::shuffle(tour.begin(), tour.end(), std::mt19937{std::random_device{}()});
    return tour;
}

std::vector<int> generateNeighbor(const std::vector<int> &tour) {
    std::vector<int> newTour = tour;
    int i = rand() % tour.size();
    int j = rand() % tour.size();
    std::swap(newTour[i], newTour[j]);  // Swap two random cities
    return newTour;
}

bool acceptNewSolution(double newCost, double oldCost, double temperature) {
    if (newCost < oldCost) {
        return true;
    }
    double acceptanceProbability = std::exp((oldCost - newCost) / temperature);
    return acceptanceProbability > static_cast<double>(rand()) / RAND_MAX;
}

std::vector<City> readCitiesFromFile(const std::string &filename) {
    std::ifstream file(filename);
    std::vector<City> cities;
    if (file.is_open()) {
        int id;
        double x, y;
        while (file >> id >> x >> y) {
            cities.push_back({id, x, y});
        }
        file.close();
    } else {
        std::cerr << "Error opening file " << filename << std::endl;
    }
    return cities;
}

int main() {
    srand(static_cast<unsigned>(time(0)));

    std::string filename = "xqf131.tsp";  // Change this as needed
    std::vector<City> cities = readCitiesFromFile(filename);
    int numCities = cities.size();

    // Simulated Annealing parameters
    double initialTemperature = 10000.0;
    double coolingRate = 0.999;
    int maxIterations = 10000;

    // Generate an initial random tour
    std::vector<int> currentTour = generateInitialTour(numCities);
    double currentCost = calculateTotalDistance(cities, currentTour);

    std::vector<int> bestTour = currentTour;
    double bestCost = currentCost;

    // Simulated Annealing loop
    double temperature = initialTemperature;
    for (int iteration = 0; iteration < maxIterations; ++iteration) {
        std::vector<int> newTour = generateNeighbor(currentTour);
        double newCost = calculateTotalDistance(cities, newTour);

        if (acceptNewSolution(newCost, currentCost, temperature)) {
            currentTour = newTour;
            currentCost = newCost;

            if (currentCost < bestCost) {
                bestTour = currentTour;
                bestCost = currentCost;
            }
        }

        temperature *= coolingRate;
    }

    // Output results
    std::cout << "Final cost: " << bestCost << std::endl;
    std::cout << "Tour sequence: ";
    for (int city : bestTour) {
        std::cout << city << " ";
    }
    std::cout << std::endl;

    return 0;
}
Place the TSP data file (e.g., xqf131.tsp) in the same directory as the code.
Compile and run the code:
g++ tsp_sa.cpp -o tsp_sa
./tsp_sa