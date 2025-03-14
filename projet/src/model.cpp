#include "model.hpp"

#include <mpi.h>

#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <cstring>

namespace {
    double pseudo_random(std::size_t index, std::size_t time_step) {
        std::uint_fast32_t xi = std::uint_fast32_t(index * (time_step + 1));
        std::uint_fast32_t r = (48271 * xi) % 2147483647;
        return r / 2147483646.;
    }

    double log_factor(std::uint8_t value) {
        return std::log(1. + value) / std::log(256);
    }
}

Model::Model(double t_length, unsigned t_discretization, std::array<double, 2> t_wind,
             LexicoIndices t_start_fire_position, double t_max_wind)
    : m_length(t_length),
      m_distance(-1),
      m_geometry(t_discretization),
      m_wind(t_wind),
      m_wind_speed(std::sqrt(t_wind[0] * t_wind[0] + t_wind[1] * t_wind[1])),
      m_max_wind(t_max_wind),
      m_vegetation_map(t_discretization * t_discretization, 255u),
      m_fire_map(t_discretization * t_discretization, 0u) {
    if (t_discretization == 0) {
        throw std::range_error("Le nombre de cases par direction doit être plus grand que zéro.");
    }
    m_distance = m_length / double(m_geometry);
    auto index = get_index_from_lexicographic_indices(t_start_fire_position);
    m_fire_front = new std::uint8_t[t_discretization * t_discretization];
    for (int i = 0; i < t_discretization * t_discretization; i++) {
        m_fire_front[i] = 0;
    }
    m_fire_map[index] = 255u;
    m_fire_front[index] = 255u;
    simulation_cells = 1;

    constexpr double alpha0 = 4.52790762e-01;
    constexpr double alpha1 = 9.58264437e-04;
    constexpr double alpha2 = 3.61499382e-05;

    if (m_wind_speed < t_max_wind)
        p1 = alpha0 + alpha1 * m_wind_speed + alpha2 * (m_wind_speed * m_wind_speed);
    else
        p1 = alpha0 + alpha1 * t_max_wind + alpha2 * (t_max_wind * t_max_wind);
    p2 = 0.3;

    if (m_wind[0] > 0) {
        alphaEastWest = std::abs(m_wind[0] / t_max_wind) + 1;
        alphaWestEast = 1. - std::abs(m_wind[0] / t_max_wind);
    } else {
        alphaWestEast = std::abs(m_wind[0] / t_max_wind) + 1;
        alphaEastWest = 1. - std::abs(m_wind[0] / t_max_wind);
    }

    if (m_wind[1] > 0) {
        alphaSouthNorth = std::abs(m_wind[1] / t_max_wind) + 1;
        alphaNorthSouth = 1. - std::abs(m_wind[1] / t_max_wind);
    } else {
        alphaNorthSouth = std::abs(m_wind[1] / t_max_wind) + 1;
        alphaSouthNorth = 1. - std::abs(m_wind[1] / t_max_wind);
    }
}

// --------------------------------------------------------------------------------------------------------------------
bool Model::update(MPI_Comm computing_comm) {
    int rank, size;
    MPI_Comm_rank(computing_comm, &rank);
    MPI_Comm_size(computing_comm, &size);

    int calc_units = size;
    size_t total_size = m_geometry * m_geometry;
    int unit_size = total_size / calc_units;
    size_t start = unit_size * rank;
    size_t end = std::min(start + unit_size, total_size);
    size_t row_sz = m_geometry;
    // std::cout << "Rank " << rank << " start " << start << " end " << end << std::endl;

    std::uint8_t *next_front;
    next_front = new std::uint8_t[m_geometry * m_geometry];
    memset(next_front, 0, m_geometry * m_geometry);

    for (size_t i = start; i < end; i++) {
        next_front[i] = m_fire_front[i];
    }
    for (size_t i = start; i < end; i++) {
        if (m_fire_front[i] == 0) {
            continue;
        }

        // Récupération de la coordonnée lexicographique de la case en feu :
        LexicoIndices coord = get_lexicographic_from_index(i);
        // if(coord.row > 100 || coord.column > 100)
        //     std::cout << coord.row << " " << coord.column << std::endl;
        // Et de la puissance du foyer
        double power = log_factor(m_fire_front[i]);

        // On va tester les cases voisines pour contamination par le feu :
        if (coord.row < m_geometry - 1) {
            double tirage = pseudo_random(i + m_time_step, m_time_step);
            double green_power = m_vegetation_map[i + m_geometry];
            double correction = power * log_factor(green_power);
            if (tirage < alphaSouthNorth * p1 * correction) {
                m_fire_map[i + m_geometry] = 255.;
                next_front[i + m_geometry] = 255;
                // if(rank == 1 && i > 30000){
                //     std::cout << "Fire front " << rank << " " << coord.row << " " << coord.column << std::endl;
                // }
            }
        }

        if (coord.row > 0) {
            double tirage = pseudo_random(i * 13427 + m_time_step, m_time_step);
            double green_power = m_vegetation_map[i - m_geometry];
            double correction = power * log_factor(green_power);
            if (tirage < alphaNorthSouth * p1 * correction) {
                m_fire_map[i - m_geometry] = 255.;
                next_front[i - m_geometry] = 255;
                // if(rank == 1){
                //     std::cout << "Fire front " << rank << " " << coord.row << " " << coord.column << std::endl;
                // }
            }
        }

        if (coord.column < m_geometry - 1) {
            double tirage = pseudo_random(i * 13427 * 13427 + m_time_step, m_time_step);
            double green_power = m_vegetation_map[i + 1];
            double correction = power * log_factor(green_power);
            if (tirage < alphaEastWest * p1 * correction) {
                m_fire_map[i + 1] = 255.;
                next_front[i + 1] = 255;
                // if(rank == 1){
                //     std::cout << "Fire front " << rank << " " << coord.row << " " << coord.column << std::endl;
                // }
            }
        }

        if (coord.column > 0) {
            double tirage = pseudo_random(i * 13427 * 13427 * 13427 + m_time_step, m_time_step);
            double green_power = m_vegetation_map[i - 1];
            double correction = power * log_factor(green_power);
            if (tirage < alphaWestEast * p1 * correction) {
                m_fire_map[i - 1] = 255.;
                next_front[i - 1] = 255;
                // if(rank == 1){
                //     std::cout << "Fire front " << rank << " " << coord.row << " " << coord.column << std::endl;
                // }
            }
        }
        // Si le feu est à son max,
        if (m_fire_front[i] == 255) {  // On regarde si il commence à faiblir pour s'éteindre au bout d'un moment :
            double tirage = pseudo_random(i * 52513 + m_time_step, m_time_step);
            if (tirage < p2) {
                m_fire_map[i] >>= 1;
                next_front[i] >>= 1;
            }
        } else {
            // Foyer en train de s'éteindre.
            m_fire_map[i] >>= 1;
            next_front[i] >>= 1;
        }
    }
    simulation_cells = 0;
    // A chaque itération, la végétation à l'endroit d'un foyer diminue
    for (size_t i = start; i < end; i++) {
        m_fire_front[i] = next_front[i];
        simulation_cells += (m_fire_front[i] != 0);
    }
    if (rank < size - 1)
        for (size_t i = end; i < end + row_sz; i++) 
            m_fire_front[i] = next_front[i];
        
    if (rank > 0)
        for (size_t i = start - row_sz; i < start; i++)
            m_fire_front[i] = next_front[i];
        

    
    for (size_t i = start; i < end; i++) {
        // std::cout << "Fire front " << rank << " " << f.first << " " << (int)f.second << std::endl;
        if (m_fire_front[i] == 0) {
            continue;
        }

        if (m_vegetation_map[i] > 0)
            m_vegetation_map[i] -= 1;
    }
    m_time_step += 1;

    if (rank < size - 1) {
        MPI_Send(m_fire_front + end, row_sz, MPI_UINT8_T, rank + 1, 12, computing_comm);
        MPI_Send(m_vegetation_map.data() + end - row_sz, row_sz, MPI_UINT8_T, rank + 1, 13, computing_comm);
        
        MPI_Recv(m_fire_front + end - row_sz, row_sz, MPI_UINT8_T, rank - 1, 10, computing_comm, MPI_STATUS_IGNORE);
        MPI_Recv(m_vegetation_map.data() + end, row_sz, MPI_UINT8_T, rank + 1, 11, computing_comm, MPI_STATUS_IGNORE);
    }

    if(rank > 0) {
        MPI_Send(m_fire_front + start - row_sz, row_sz, MPI_UINT8_T, rank - 1, 10, computing_comm);
        MPI_Send(m_vegetation_map.data() + start, row_sz, MPI_UINT8_T, rank - 1, 11, computing_comm);

        MPI_Recv(m_fire_front + start, row_sz, MPI_UINT8_T, rank - 1, 12, computing_comm, MPI_STATUS_IGNORE);
        MPI_Recv(m_vegetation_map.data() + start - row_sz, row_sz, MPI_UINT8_T, rank - 1, 13, computing_comm, MPI_STATUS_IGNORE);
    }
    
    return simulation_cells > 0;
}
// ====================================================================================================================
std::size_t
Model::get_index_from_lexicographic_indices(LexicoIndices t_lexico_indices) const {
    return t_lexico_indices.row * this->geometry() + t_lexico_indices.column;
}
// --------------------------------------------------------------------------------------------------------------------
auto Model::get_lexicographic_from_index(std::size_t t_global_index) const -> LexicoIndices {
    LexicoIndices ind_coords;
    ind_coords.row = t_global_index / this->geometry();
    ind_coords.column = t_global_index % this->geometry();
    return ind_coords;
}