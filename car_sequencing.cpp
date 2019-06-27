//
// Created by guardian on 16.06.19.
//
#include <iostream>
#include <vector>
#include <random>
#include <climits>
#include <cfloat>
#include <ctime>

struct robot {
    size_t window_capacity = 0; // сколько не более чем это число автомобилей могут содержать опцию i
    size_t window_size = 0; // размер окна робота (т.е. сколько автомобилей попадают в окно этого робота)

    default robot() {}
};

size_t get_integer_random(size_t N) {
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_int_distribution<size_t> distr(0, N - 1); // define the range
    return distr(eng);
}

long long
get_full_sequence_error(const std::vector<std::vector<size_t>> &car_type, const std::vector<size_t> &car_sequence,
                        const std::vector<robot> &robots) {
    long long func_value = 0;
    for (size_t i = 0; i != robots.size(); ++i) {// для всех роботов
        size_t window = robots[i].window_size; // смотрим на размер окна i-го робота
        int w_left = -window;
        int w_right = 0;
        while ((w_right < car_sequence.size() + window)) {
            long long counter = 0;
            if (w_left < 0) {
                for (int car = 0; car != w_right + 1; ++car) {
                    counter += car_type[car_sequence[car]][i];
                }
            } else if (w_left == 0) { ;
            } else {
                for (int car = w_left; car < w_right && car < car_sequence.size(); ++car) {
                    counter += car_type[car_sequence[car]][i];
                }
            }
            if (counter > robots[i].window_capacity)
                func_value += (counter - robots[i].window_capacity);
            ++w_right;
            ++w_left;
        }
    }
    return func_value;
}

// подаем ему на вход сначала старую последовательность с указанием индексов, а потом новую
long long
get_part_sequence_error(const std::vector<std::vector<size_t>> &car_type, const std::vector<size_t> &car_sequence,
                        const std::vector<robot> &robots) {
    long long func_value = 0;
    for (size_t i = 0; i != robots.size(); ++i) {// для всех роботов
        size_t window = robots[i].window_size; // смотрим на размер окна i-го робота
        int w_left = -window;
        int w_right = 0;
        while ((w_right < car_sequence.size() + window) && (w_left <= index1 && index1 <= w_right)) {
            if (!(w_left <= index1 && index1 <= w_right)) continue;
            long long counter = 0;
            if (w_left < 0) {
                for (int car = 0; car != w_right + 1; ++car) {
                    counter += car_type[car_sequence[car]][i];
                }
            } else if (w_left == 0) { ;
            } else {
                for (int car = w_left; car < w_right && car < car_sequence.size(); ++car) {
                    counter += car_type[car_sequence[car]][i];
                }
            }
            if (counter > robots[i].window_capacity)
                func_value += (counter - robots[i].window_capacity);
            ++w_right;
            ++w_left;
        }
    }
    return func_value;
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


long long
get_greedy_decision_function(const std::vector<std::vector<size_t>> &car_type, const std::vector<size_t> &car_sequence,
                             size_t cur_car, const std::vector<long long> &counter_in_all_cars,
                             const std::vector<long long> &counter_in_sequence) {
    long long sum = 0;
    for (size_t option = 1; option != car_type[0].size(); ++option) {
        // считаем для последовательности
        sum += ((car_type[cur_car][option] == 0) ^
                (counter_in_all_cars[option] * 1.0 / car_type.size() >
                 counter_in_sequence[option] * 1.0 / car_sequence.size()));
    }
    return sum;
}

void update_counter_in_sequence(const std::vector<std::vector<size_t>> &car_type,
                                std::vector<long long> &counter_in_sequence,
                                size_t new_car) {
    for (size_t option = 1; option != car_type[0].size(); ++option) {
        counter_in_sequence[option] += car_type[new_car][option];
    }
}


size_t get_heuristic_best_car(const std::vector<size_t> &min_cars, const std::vector<std::vector<size_t>> &car_type,
                              const std::vector<robot> &robots, const std::vector<size_t> &car_sequence,
                              const std::vector<long long> &counter_in_all_cars,
                              const std::vector<long long> &counter_in_sequence) {
    // выбираем первую машину среди самых больших по требуемости опций машин
    if (car_sequence.empty()) {
        size_t max_car = 0;
        long long cur_max = 0;
        for (auto car : min_cars) {
            long long cur_counter = 0;
            for (size_t j = 1; j != car_type[car].size(); ++j) {
                cur_counter += car_type[car][j];
            }
            if (cur_counter > cur_max) {
                max_car = car;
                cur_max = cur_counter;
            }
        }
        return max_car;
    }
    // дальше жадным способом дополняем решение
    long long max_err_i = 0;
    size_t max_car_i = 0;
    for (auto car : min_cars) {
        auto cur_err = get_greedy_decision_function(car_type, car_sequence, car, counter_in_all_cars,
                                                    counter_in_sequence);
        if (cur_err > max_err_i) {
            max_err_i = cur_err;
            max_car_i = car;
        }
    }
    return max_car_i;
}

int get_violation(const std::vector<std::vector<size_t>> &car_type, const std::vector<size_t> &car_sequence,
                  const std::vector<robot> &robots, size_t new_car, size_t option) {
    long long counter = 0;
    if (car_sequence.size() + 1 < robots[option].window_size) {
        for (auto car : car_sequence) {
            counter += car_type[car][option];
        }
        counter += car_type[new_car][option];
    } else {
        for (size_t i = car_sequence.size() + 1 - robots[option].window_size; i < car_sequence.size(); ++i) {
            counter += car_type[car_sequence[i]][option];
        }
        counter += car_type[new_car][option];
    }
    if (counter <= robots[option].window_capacity) return 0;
    else return 1;
}

std::vector<size_t> get_greedy_sol(size_t N, const std::vector<std::vector<size_t>> &car_type,
                                   const std::vector<robot> &robots, const std::vector<size_t> &size_car_type) {
    // тут заранее посчитаем штуки для эврестической функции
    // для всех типов машин
    std::vector<long long> counter_in_all_cars(car_type[0].size(), 0);
    for (size_t option = 0; option != car_type[0].size(); ++option) {
        for (auto car : car_type) counter_in_all_cars[option] += car[option];
    }
    // для нашей последовательности
    std::vector<long long> counter_in_sequence(car_type[0].size(), 0);
    //
    // суммируем кол-во машин каждого типа
    std::vector<size_t> num_of_cars(size_car_type.size(), 0);
    std::vector<size_t> car_sequence(0);
    std::vector<size_t> min_cars;
    for (size_t i = 0; i != N; ++i) { // нужно добавить N типов машин к последовательности
        long long min_sum = LLONG_MAX; // минимальная сумма на i-м шаге
        min_cars = {}; // машины с минимальной суммой
        for (size_t car = 0; car != car_type.size(); ++car) { // для каждой машины
            if (num_of_cars[car] >= size_car_type[car]) continue;
            long long car_sum = 0; // считаем сумму
            for (size_t option = 0; option != car_type[0].size(); ++option) { // ошибку violation
                car_sum += car_type[car][option] * get_violation(car_type, car_sequence, robots, car, option);
            }
            if (car_sum < min_sum) {
                min_sum = car_sum;
                min_cars = {car};
            } else if (car_sum == min_sum) {
                min_cars.push_back(car);
            }
        }
        if (min_cars.size() > 1) {
            auto res_car = get_heuristic_best_car(min_cars, car_type, robots,
                                                  car_sequence, counter_in_all_cars,
                                                  counter_in_sequence);
            car_sequence.push_back(res_car);
            update_counter_in_sequence(car_type, counter_in_sequence, res_car);
            ++num_of_cars[res_car];
        } else {
            car_sequence.push_back(min_cars[0]);
            update_counter_in_sequence(car_type, counter_in_sequence, min_cars[0]);
            ++num_of_cars[min_cars[0]];
        }
    }
    return car_sequence;
}


int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    clock_t _start = clock();
    // кол-во машин в линии сборки
    size_t N = 0;
    std::cin >> N;
    // кол-во различных опций
    size_t L = 0;
    std::cin >> L;
    // кол-во различных типо машин
    size_t K = 0;
    std::cin >> K;
    std::vector<robot> robots(L);
    // максимально возможное допустимое число обрабатываемых автомобилей с опцией i на окне для каждого из роботов
    for (size_t i = 0; i != L; ++i) {
        std::cin >> robots[i].window_capacity;
    }
    // размеры окон роботов
    for (size_t i = 0; i != L; ++i) {
        std::cin >> robots[i].window_size;
    }
    std::vector<size_t> size_car_type(K, 0);
    // описание типов автомобилей (номер типа: кол-во изготавливаемых, дальше l чисел из {0, 1} из опций (есть/нет опция)
    std::vector<std::vector<size_t>> car_type(K, std::vector<size_t>(L));
    for (size_t i = 0; i != K; ++i) {
        std::cin >> size_car_type[i]; // кол-во изготовления машин i-го типа
        for (size_t j = 0; j != L; ++j) {
            std::cin >> car_type[i][j];
        }
    }
    // само решение
    // инициализация отжог
    double temperature = std::sqrt(N);
    double temperature_save = temperature;
    double alpha = 0.999985;
    double stopping_temperature = 1e-8;
    size_t stopping_iteration = 200;
    size_t iteration = 0;
//    // лучшее решение
    std::vector<size_t> best_solution;
    // первое решение жадным
    auto cur_solution = get_greedy_sol(N, car_type, robots, size_car_type);
    auto cur_err = get_full_sequence_error(car_type, cur_solution, robots);
    // соответственно они и будут пока лучшими решениями
    best_solution = cur_solution;
    auto best_err = cur_err;
    // дальше будем пытаться его улучшить
    while (static_cast<double>(clock() - _start) / CLOCKS_PER_SEC < 5.6 && temperature >= stopping_temperature) {
//    while (static_cast<double>(clock() - _start) / CLOCKS_PER_SEC < 5.9) {
        // ищем новое решение
        size_t car1 = 0;
        size_t car2 = 0;
        while (car1 >= car2) {
            car1 = get_integer_random(N);
            car2 = get_integer_random(N);
        }
        auto potential_sol = cur_solution;
        std::swap(potential_sol[car1], potential_sol[car2]);
        auto err1 = get_part_sequence_error(car_type, cur_solution, robots, cur_err, car1, car2);
        auto err2 = get_part_sequence_error(car_type, potential_sol, robots, cur_err, car1, car2);
        auto potential_dist = cur_err - err1 + err2;
        if (potential_dist < cur_err) {
            cur_err = potential_dist;
            cur_solution = potential_sol;
            if (potential_dist < best_err) {
                best_err = potential_dist;
                best_solution = potential_sol;
            }
        }
        else {
            if (get_double_random() < get_probability(cur_err, potential_dist, temperature)) {
                cur_err = potential_dist;
                cur_solution = potential_sol;
            }
        }
        temperature *= alpha;
        iteration += 1;
    }
    std::cout << best_err << "\n";
    for (auto elem : best_solution) std::cout << elem << " ";
}