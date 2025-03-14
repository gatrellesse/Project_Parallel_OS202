#include <mpi.h>

#include <cassert>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>

#include "display.hpp"
#include "model.hpp"

using namespace std::string_literals;
using namespace std::chrono_literals;

struct ParamsType {
    double length{1.};
    unsigned discretization{20u};
    std::array<double, 2> wind{0., 0.};
    Model::LexicoIndices start{10u, 10u};
};

void analyze_arg(int nargs, char* args[], ParamsType& params) {
    if (nargs == 0) return;
    std::string key(args[0]);
    if (key == "-l"s) {
        if (nargs < 2) {
            std::cerr << "Manque une valeur pour la longueur du terrain !"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
        params.length = std::stoul(args[1]);
        analyze_arg(nargs - 2, &args[2], params);
        return;
    }
    auto pos = key.find("--longueur=");
    if (pos < key.size()) {
        auto subkey = std::string(key, pos + 11);
        params.length = std::stoul(subkey);
        analyze_arg(nargs - 1, &args[1], params);
        return;
    }

    if (key == "-n"s) {
        if (nargs < 2) {
            std::cerr << "Manque une valeur pour le nombre de cases par "
                         "direction pour la discrétisation du terrain !"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
        params.discretization = std::stoul(args[1]);
        analyze_arg(nargs - 2, &args[2], params);
        return;
    }
    pos = key.find("--number_of_cases=");
    if (pos < key.size()) {
        auto subkey = std::string(key, pos + 18);
        params.discretization = std::stoul(subkey);
        analyze_arg(nargs - 1, &args[1], params);
        return;
    }

    if (key == "-w"s) {
        if (nargs < 2) {
            std::cerr
                << "Manque une paire de valeurs pour la direction du vent !"
                << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values = std::string(args[1]);
        params.wind[0] = std::stod(values);
        auto pos = values.find(",");
        if (pos == values.size()) {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule "
                         "pour définir la vitesse"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(values, pos + 1);
        params.wind[1] = std::stod(second_value);
        analyze_arg(nargs - 2, &args[2], params);
        return;
    }
    pos = key.find("--wind=");
    if (pos < key.size()) {
        auto subkey = std::string(key, pos + 7);
        params.wind[0] = std::stoul(subkey);
        auto pos = subkey.find(",");
        if (pos == subkey.size()) {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule "
                         "pour définir la vitesse"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(subkey, pos + 1);
        params.wind[1] = std::stod(second_value);
        analyze_arg(nargs - 1, &args[1], params);
        return;
    }

    if (key == "-s"s) {
        if (nargs < 2) {
            std::cerr << "Manque une paire de valeurs pour la position du "
                         "foyer initial !"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values = std::string(args[1]);
        params.start.column = std::stod(values);
        auto pos = values.find(",");
        if (pos == values.size()) {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule "
                         "pour définir la position du foyer initial"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(values, pos + 1);
        params.start.row = std::stod(second_value);
        analyze_arg(nargs - 2, &args[2], params);
        return;
    }
    pos = key.find("--start=");
    if (pos < key.size()) {
        auto subkey = std::string(key, pos + 8);
        params.start.column = std::stoul(subkey);
        auto pos = subkey.find(",");
        if (pos == subkey.size()) {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule "
                         "pour définir la vitesse"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(subkey, pos + 1);
        params.start.row = std::stod(second_value);
        analyze_arg(nargs - 1, &args[1], params);
        return;
    }
}

ParamsType parse_arguments(int nargs, char* args[]) {
    if (nargs == 0) return {};
    if ((std::string(args[0]) == "--help"s) || (std::string(args[0]) == "-h")) {
        std::cout <<
            R"RAW(Usage : simulation [option(s)]
  Lance la simulation d'incendie en prenant en compte les [option(s)].
  Les options sont :
    -l, --longueur=LONGUEUR     Définit la taille LONGUEUR (réel en km) du carré représentant la carte de la végétation.
    -n, --number_of_cases=N     Nombre n de cases par direction pour la discrétisation
    -w, --wind=VX,VY            Définit le vecteur vitesse du vent (pas de vent par défaut).
    -s, --start=COL,ROW         Définit les indices I,J de la case où commence l'incendie (milieu de la carte par défaut)
)RAW";
        exit(EXIT_SUCCESS);
    }
    ParamsType params;
    analyze_arg(nargs, args, params);
    return params;
}

bool check_params(ParamsType& params) {
    bool flag = true;
    if (params.length <= 0) {
        std::cerr << "[ERREUR FATALE] La longueur du terrain doit être "
                     "positive et non nulle !"
                  << std::endl;
        flag = false;
    }

    if (params.discretization <= 0) {
        std::cerr << "[ERREUR FATALE] Le nombre de cellules par direction doit "
                     "être positive et non nulle !"
                  << std::endl;
        flag = false;
    }

    if ((params.start.row >= params.discretization) ||
        (params.start.column >= params.discretization)) {
        std::cerr << "[ERREUR FATALE] Mauvais indices pour la position "
                     "initiale du foyer"
                  << std::endl;
        flag = false;
    }

    return flag;
}

void display_params(ParamsType const& params) {
    std::cout << "Parametres définis pour la simulation : \n"
              << "\tTaille du terrain : " << params.length << std::endl
              << "\tNombre de cellules par direction : "
              << params.discretization << std::endl
              << "\tVecteur vitesse : [" << params.wind[0] << ", "
              << params.wind[1] << "]" << std::endl
              << "\tPosition initiale du foyer (col, ligne) : "
              << params.start.column << ", " << params.start.row << std::endl;
}

int main(int nargs, char* args[]) {
    auto params = parse_arguments(nargs - 1, &args[1]);
    display_params(params);
    if (!check_params(params)) return EXIT_FAILURE;

    int rank, world;
    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &world);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::shared_ptr<Displayer> displayer = nullptr;
    Model simu(params.length, params.discretization, params.wind, params.start);
    size_t map_sz = simu.fire_map().size();
    simu.reset_time_step();

    if (rank == 0) {
        displayer = Displayer::init_instance(params.discretization,
                                             params.discretization);
    }

    // auto displayer = Displayer::init_instance( params.discretization,
    // params.discretization );

    // if (rank == 1) {
    //     simu = Model(params.length, params.discretization, params.wind, params.start);
    // }
    MPI_Request is_running_request, fire_map_request, vegetal_map_request;
    int32_t isRunning = 1;

    if (rank == 1) {
        MPI_Isend(&isRunning, 1, MPI_INT32_T, 0, 0, MPI_COMM_WORLD, &is_running_request);
        MPI_Isend(simu.fire_map().data(), map_sz, MPI_UINT8_T, 0, 1, MPI_COMM_WORLD, &fire_map_request);
        MPI_Isend(simu.vegetal_map().data(), map_sz, MPI_UINT8_T, 0, 2, MPI_COMM_WORLD, &vegetal_map_request);
    }

    SDL_Event event;
    auto simulation_start = std::chrono::high_resolution_clock::now();
    double avg_update_time = 0.;
    double avg_display_time = 0.;
    int frame_count = 0;

    // MPI_Group computing_group;
    MPI_Comm computing_comm;
    MPI_Comm_split(MPI_COMM_WORLD, rank != 0, rank, &computing_comm);

    int row_rank, row_size;
    MPI_Comm_rank(computing_comm, &row_rank);
    MPI_Comm_size(computing_comm, &row_size);

    for (auto [k, v] : simu.fire_front()) {
        std::cout << k << " " << v << std::endl;
    }

    // MPI_Comm_free(&computing_comm);
    // std::this_thread::sleep_for(1s);

    while (isRunning) {
        if (rank != 0) {
            auto start_time = std::chrono::high_resolution_clock::now();
            // atualizar fantasmas
            simu.update(computing_comm);
            MPI_Request sync;
            if (rank != 1){
                int front_size = simu.fire_front().size();
                // std::cout << "Rank " << rank << " sending " << front_size << std::endl;
                MPI_Send(&front_size, 1, MPI_INT32_T, 1, 12, MPI_COMM_WORLD);
                std::vector<std::int32_t> front(front_size * 2);
                for (int i = 0, j = 0; i < front_size; i++, j += 2) {
                    front[j] = simu.fire_front()[i].first;
                    front[j + 1] = (int)simu.fire_front()[i].second;
                }
                MPI_Send(front.data(), front_size * 2, MPI_INT32_T, 1, 13, MPI_COMM_WORLD);

                // std::cout << "Rank " << rank << " waiting updated front" << std::endl;
                int updated_front_size;
                MPI_Recv(&updated_front_size, 1, MPI_INT32_T, 1, 12, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                std::vector<std::int32_t> fire_front(2 * updated_front_size);
                MPI_Recv(fire_front.data(), updated_front_size * 2, MPI_INT32_T, 1, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                    // std::cout << rank << " received " << k << " " << v << std::endl;
                // std::cout << "received ~" << rank << "~ " << fire_front[0] << " " << fire_front[2] << " " << fire_front[4] << " " << fire_front[6] << std::endl;
                // std::cout << "Rank " << rank << " received " << fire_front.size() << std::endl;
                simu.set_fire_front(fire_front);
            } else {
                for (int i = 2; i < world; i++) {
                    int front_size;
                    MPI_Recv(&front_size, 1, MPI_INT32_T, i, 12, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    std::vector<std::int32_t> fire_front(front_size * 2);
                    MPI_Recv(fire_front.data(), front_size * 2, MPI_INT32_T, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    // std::cout << "Received fire front from " << i << " of size " << front_size << std::endl;
                    // std::cout << "received ~" << i << "~ " << fire_front[0] << " " << fire_front[2] << " " << fire_front[4] << " " << fire_front[6] << std::endl;
                    // std::cout << "-------\n";
                    simu.update_fire_front(fire_front);
                    // for (auto [k, v] : simu.fire_front()) {
                    //     std::cout << k << " " << v << std::endl;
                    // }
                }
                int updated_front_size = simu.fire_front().size();
                auto aux = simu.fire_front();
                std::vector<std::int32_t> upd_front(updated_front_size * 2);
                for (int i = 0, j = 0; i < updated_front_size; i++, j+= 2) {
                    // std::cout << "mounting " << aux[i].first << " " << (int)aux[i].second << std::endl;
                    upd_front[j] = aux[i].first;
                    upd_front[j + 1] = (int)aux[i].second;
                    // std::cout << "mounted " << upd_front[j] << " " << upd_front[j + 1] << std::endl;
                }
                
                // std::cout << "sending upd ~" << rank << "~ " << upd_front[0] << " " << upd_front[2] << " " << upd_front[4] << " " << upd_front[6] << std::endl;
                
                for (int i = 2;i < world; i++) {
                    MPI_Send(&updated_front_size, 1, MPI_INT32_T, i, 12, MPI_COMM_WORLD);
                    MPI_Send(upd_front.data(), updated_front_size * 2, MPI_INT32_T, i, 13, MPI_COMM_WORLD);
                    // std::cout << "Sent updated fire front to " << i << " of size " << updated_front_size << std::endl;
                }
            }

            // std::cout << "Finished sync " << rank << " " << simu.fire_front().size() << std::endl;
            // for (auto [k, v] : simu.fire_front()) {
            //     std::cout << "~" << rank << "~ " << k << " " << v << std::endl;
            // }
            isRunning = simu.fire_front().size() > 0;


            auto end_time = std::chrono::high_resolution_clock::now();
            avg_update_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();

            // std::cout << "Time : "
            //           << std::chrono::duration_cast<std::chrono::milliseconds>(
            //                  end_time - start_time)
            //                  .count()
            //           << " ms" << std::endl;

            // std::cout << isRunning << " " << map_sz << std::endl;
            if (!isRunning) {
                for (int i = 0; i < 1000; i++) {
                    MPI_Isend(&isRunning, 1, MPI_INT32_T, 0, 0, MPI_COMM_WORLD, &is_running_request);
                    MPI_Isend(simu.fire_map().data(), map_sz, MPI_UINT8_T, 0, 1, MPI_COMM_WORLD, &fire_map_request);
                    MPI_Isend(simu.vegetal_map().data(), map_sz, MPI_UINT8_T, 0, 2, MPI_COMM_WORLD, &vegetal_map_request);
                }
                std::cout << "Stopping simulation" << std::endl;
                break;
            }

            if (rank == 1){

                int is_running_received, is_fire_map_received, is_vegetal_map_received;
                MPI_Test(&is_running_request, &is_running_received, MPI_STATUS_IGNORE);
                MPI_Test(&fire_map_request, &is_fire_map_received, MPI_STATUS_IGNORE);
                MPI_Test(&vegetal_map_request, &is_vegetal_map_received, MPI_STATUS_IGNORE);
                
                if (is_running_received)
                MPI_Isend(&isRunning, 1, MPI_INT32_T, 0, 0, MPI_COMM_WORLD, &is_running_request);
                if (is_fire_map_received)
                MPI_Isend(simu.fire_map().data(), map_sz, MPI_UINT8_T, 0, 1, MPI_COMM_WORLD, &fire_map_request);
                if (is_vegetal_map_received)
                MPI_Isend(simu.vegetal_map().data(), map_sz, MPI_UINT8_T, 0, 2, MPI_COMM_WORLD, &vegetal_map_request);
                
                if ((simu.time_step() & 31) == 0) {
                    std::cout << "Time step " << simu.time_step()
                    << "\n===============" << std::endl;
                    std::cout.flush();
                }
                
            }
            
            // std::this_thread::sleep_for(1s);
        }

        if (rank == 0) {
            std::vector<uint8_t> fire_map(map_sz), vegetal_map(map_sz);
            MPI_Recv(&isRunning, 1, MPI_INT32_T, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // std::cout << "Received isRunning : " << isRunning << std::endl;
            if (!isRunning) break;
            MPI_Recv(fire_map.data(), map_sz, MPI_UINT8_T, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // std::cout << "Received fire_map " << fire_map.size() << std::endl;
            MPI_Recv(vegetal_map.data(), map_sz, MPI_UINT8_T, 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // std::cout << "Received vegetal_map " << vegetal_map.size() << std::endl;

            if (fire_map.size() != map_sz) {
                std::cerr << "Fire map size mismatch" << std::endl;
                return EXIT_FAILURE;
            }
            if (vegetal_map.size() != map_sz) {
                std::cerr << "Vegetal map size mismatch" << std::endl;
                return EXIT_FAILURE;
            }

            // std::cout << "Updating" << std::endl;
            auto start_time = std::chrono::high_resolution_clock::now();
            displayer->update(vegetal_map, fire_map);
            frame_count++;
            auto end_time = std::chrono::high_resolution_clock::now();
            avg_display_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
            // std::cout << "Finished updating" << std::endl;
            // displayer->update(simu.vegetal_map(), simu.fire_map());
            // std::this_thread::sleep_for(0.1s);
        }
        if (SDL_PollEvent(&event) && event.type == SDL_QUIT) break;
    }
    if (rank == 1) {
        std::cout << "Average update time: " << avg_update_time / simu.time_step() / 1000 << " ms" << std::endl;
    } 
    if (rank == 0) {
        std::cout << "Average display time: " << avg_display_time / frame_count / 1000 << " ms" << std::endl;
    }

    std::cout << "Simulation finished " << rank << std::endl;
    auto simulation_end = std::chrono::high_resolution_clock::now();
    auto simulation_duration = std::chrono::duration_cast<std::chrono::milliseconds>(simulation_end - simulation_start).count();
    if (rank == 1)
        std::cout << "Simulation time (n = " << params.discretization << "): " << simulation_duration << " ms" << std::endl;

    MPI_Finalize();
    return EXIT_SUCCESS;
}