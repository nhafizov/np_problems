#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <ctime>
#include <random>
#include <cmath>

// точка
struct Point {
    double x = 0;
    double y = 0;

    Point() {}

    Point(double _x, double _y) : x(_x), y(_y) {}
};

// считаем расстояние между двумя точками
double get_euclid_distance(const Point &first, const Point &second) {
    return std::sqrt((first.x - second.x) * (first.x - second.x) +
                     (first.y - second.y) * (first.y - second.y));
}

// общее расстояние пути
double get_tour_Distance(const std::vector<Point> &listOfPoints, const std::vector<size_t> &tour) {
    double sum_dist = 0;
    for (auto i = 0; i != tour.size() - 1; ++i) {
        sum_dist += get_euclid_distance(listOfPoints[tour[i]], listOfPoints[tour[i + 1]]);
    }
    sum_dist += get_euclid_distance(listOfPoints[tour[0]], listOfPoints[tour[tour.size() - 1]]);
    return sum_dist;
}

// генератор случайных чисел из равномерного распределения
size_t get_integer_random(size_t N) {
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_int_distribution<size_t> distr(0, N - 1); // define the range
    return distr(eng);
}

double get_double_random() {
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_real_distribution<double> distr(0, 1); // define the range
    return distr(eng);
};

double get_probability(double prev_score, double new_score, double temperature) {
    return std::exp(std::abs(new_score - prev_score) / temperature);
}

// первое жадное решение
std::vector<size_t> get_greedy_sol(const std::vector<Point> &listOfPoints) {
    // жадное решение
    std::vector<size_t> greedy_sol(listOfPoints.size(), 0);
    // использованные точки
    std::vector<bool> used(listOfPoints.size(), false);
    greedy_sol[0] = 0;
    // вектор, в котором отмечаем посещенные точки
    used[0] = true;
    // кол-во точек в решении
    size_t cur_size = 1;
    while (cur_size != listOfPoints.size()) {
        // расстояние до ближайшей точки
        double cur_nearest_dist = DBL_MAX;
        // индекс ближайшей точки
        size_t cur_nearest_city = SIZE_MAX;
        // последний элемент в массиве
        size_t last_elem = greedy_sol.back();
        // цикл по всем точкам
        for (size_t i = 0; i != listOfPoints.size(); ++i) {
            // пропускаем посещенные
            if (used[i]) continue;
            // считаем расстояние от последней точки в пути до нового
            double tmp_dist = get_euclid_distance(listOfPoints[last_elem], listOfPoints[i]);
            // если новый город ближе, обновляем
            if (tmp_dist < cur_nearest_dist) {
                cur_nearest_city = i;
                cur_nearest_dist = tmp_dist;
            }
        }
        // добавляем новую точку в путь
        greedy_sol[cur_size] = cur_nearest_city;
        // отмечаем, что он пройден
        used[cur_nearest_city] = true;
        // обновляем кол-во точек в пути
        ++cur_size;
    }
    return greedy_sol;
}

std::vector<size_t> new_solution(const std::vector<Point> &listOfPoints, const std::vector<size_t> &prev_sol) {
    std::vector<size_t> new_solution = prev_sol;
    auto city1 = get_integer_random(listOfPoints.size());
    auto city2 = get_integer_random(listOfPoints.size());
    if (city1 < city2) {
        for (size_t i = city1; i <= city2; ++i) {
            new_solution[i] = prev_sol[city2 - i + city1];
        }
    } else if (city1 > city2) {
        for (size_t i = city2; i <= city1; ++i) {
            new_solution[i] = prev_sol[city1 - i + city2];
        }
    }
    return new_solution;
}

int main() {
    clock_t _start = clock();
    size_t N = 0;
    std::cin >> N;
    if (N == 0) return 0;
    std::vector<Point> listOfPoints(N);
    for (size_t i = 0; i != N; ++i) {
        std::cin >> listOfPoints[i].x >> listOfPoints[i].y;
    }
    // инициализация отжог
    double temperature = 100000;
    auto temperature_0 = temperature;
    double temperature_save = temperature;
    double alpha = 0.998;
    double stopping_temperature = 1e-8;
    size_t stopping_iteration = 1000000;
    size_t iteration = 0;
    // лучшее решение
    std::vector<size_t> best_solution;
    double best_dist = DBL_MAX;
    // первое решение жадным
    auto cur_solution = get_greedy_sol(listOfPoints);
    auto cur_dist = get_tour_Distance(listOfPoints, cur_solution);
    // соответственно они и будут пока лучшими решениями
    best_solution = cur_solution;
    best_dist = cur_dist;
    // дальше будем пытаться его улучшить
    while (static_cast<double>(clock() - _start) / CLOCKS_PER_SEC < 9.5 && temperature >= stopping_temperature) {
        // ищем новое решение
        auto potential_sol = new_solution(listOfPoints, best_solution);
        auto potential_dist = get_tour_Distance(listOfPoints, potential_sol);
        if (potential_dist < cur_dist) {
            cur_dist = potential_dist;
            cur_solution = potential_sol;
            if (potential_dist < best_dist) {
                best_dist = potential_dist;
                best_solution = potential_sol;
            }
        } else {
            if (get_double_random() < get_probability(cur_dist, potential_dist, temperature)) {
                cur_dist = potential_dist;
                cur_solution = potential_sol;
            }
        }
        temperature = temperature_0 * alpha;
        iteration += 1;
    }
    std::cout << best_dist << "\n";
    for (auto elem : best_solution) std::cout << elem + 1 << " ";
}
