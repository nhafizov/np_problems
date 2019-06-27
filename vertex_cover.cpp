//
// Created by guardian on 20.05.19.
//
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <set>
#include <iterator>

struct Edge {
    size_t from = 0;
    size_t to = 0;
    long long price = 0;

    Edge() {}

    Edge(size_t _from, size_t _to, double _price) : from(_from), to(_to), price(_price) {}
};


int main() {
    // число вершин в первой доле
    size_t N = 0;
    std::cin >> N;
    // число вершин во второй доле
    size_t M = 0;
    std::cin >> M;
    // веса вершин
    std::vector<size_t> vertex_costs(N);
    for (size_t i = 0; i != N; ++i) {
        std::cin >> vertex_costs[i];
    }
    // граф в виде ребер
    std::vector<Edge> graphE;
    for (size_t i = 0; i != M; ++i) {
        size_t from = 0;
        size_t to = 0;
        std::cin >> from >> to;
        graphE.emplace_back(Edge(--from, --to, 0));
    }

    std::set<size_t> vertex_cover;
    long long sum_costs = 0;
    std::vector<size_t> copy_vertex_costs = vertex_costs;
    for (auto &edge : graphE) {
        if (vertex_costs[edge.from] < vertex_costs[edge.to]) {
            if (vertex_costs[edge.from] == 0) continue;
            vertex_costs[edge.to] -= vertex_costs[edge.from];
            vertex_costs[edge.from] -= vertex_costs[edge.from];
            if (vertex_cover.find(edge.from) == vertex_cover.end()) edge.price += copy_vertex_costs[edge.from];
            vertex_cover.insert(edge.from);
        } else {
            if (vertex_costs[edge.to] == 0) continue;
            vertex_costs[edge.to] -= vertex_costs[edge.to];
            vertex_costs[edge.from] -= vertex_costs[edge.to];
            if (vertex_cover.find(edge.to) == vertex_cover.end()) edge.price += copy_vertex_costs[edge.to];
            vertex_cover.insert(edge.to);
        }
    }
    for (auto &edge : graphE) {
        sum_costs += edge.price;
    }
    // сумма весов
    std::cout << vertex_cover.size() << " " << sum_costs << "\n";
    for (auto &elem : vertex_cover) {
        std::cout << elem + 1 << " ";
    }
}
