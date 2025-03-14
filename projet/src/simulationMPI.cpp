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
    std::size_t simu_size = params.discretization * params.discretization;
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
    
    SDL_Event event;
    auto simulation_start = std::chrono::high_resolution_clock::now();
    double avg_update_time = 0.;
    double avg_display_time = 0.;
    int frame_count = 0;

    // MPI_Group computing_group;
    MPI_Comm computing_comm;
    MPI_Comm_split(MPI_COMM_WORLD, rank != 0, rank, &computing_comm);

    int comp_rank, comp_size;
    MPI_Comm_rank(computing_comm, &comp_rank);
    MPI_Comm_size(computing_comm, &comp_size);
    
    int unit_size = simu_size / comp_size;
    size_t start = unit_size * comp_rank;
    size_t row_sz = params.discretization;
    size_t end = std::min(start + unit_size, simu_size);
    
    // if (rank != 0) {
    //     MPI_Isend(simu.fire_map().data() + start, end - start + 1, MPI_UINT8_T, 0, 101 * rank, MPI_COMM_WORLD, &fire_map_request);
    //     MPI_Isend(simu.vegetal_map().data() + start, end - start + 1, MPI_UINT8_T, 0, 102 * rank, MPI_COMM_WORLD, &vegetal_map_request);
    //     std::cout << "Rank " << rank << " sent first maps, tags = " << 101 * rank << " " << 102 * rank << std::endl;
    // }


    while (isRunning) {
        if (rank != 0) {
            auto start_time = std::chrono::high_resolution_clock::now();
            isRunning = simu.update(computing_comm);

            // isRunning = simu.simulation_cells > 0;
            auto end_time = std::chrono::high_resolution_clock::now();
            avg_update_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();

            // std::cout << "Time : "
            //           << std::chrono::duration_cast<std::chrono::milliseconds>(
            //                  end_time - start_time)
            //                  .count()
            //           << " ms" << std::endl;

            // std::cout << isRunning << " " << map_sz << std::endl;
        
            // int is_running_received, is_fire_map_received, is_vegetal_map_received;
            // // MPI_Test(&is_running_request, &is_running_received, MPI_STATUS_IGNORE);
            // MPI_Test(&fire_map_request, &is_fire_map_received, MPI_STATUS_IGNORE);
            // MPI_Test(&vegetal_map_request, &is_vegetal_map_received, MPI_STATUS_IGNORE);

            // if (is_running_received)
                MPI_Send(&isRunning, 1, MPI_INT32_T, 0, 100 * rank, MPI_COMM_WORLD);
            // if (is_fire_map_received){
                std::vector<uint8_t> fire_map(simu.fire_map());
                MPI_Send(fire_map.data() + start, end - start + 1, MPI_UINT8_T, 0, 101 * rank, MPI_COMM_WORLD);
            // }
            // if (is_vegetal_map_received){
                std::vector<uint8_t> vegetal_map(simu.vegetal_map());
                MPI_Send(vegetal_map.data() + start, end - start + 1, MPI_UINT8_T, 0, 102 * rank, MPI_COMM_WORLD);
            // }

            MPI_Recv(&isRunning, 1, MPI_INT32_T, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            if ((simu.time_step() & 31) == 0) {
                std::cout << "Time step " << simu.time_step()
                            << "\n===============" << std::endl;
                std::cout.flush();
            }

            // std::this_thread::sleep_for(1s);
        }

        if (rank == 0) {
            std::vector<uint8_t> fire_map(map_sz), vegetal_map(map_sz);
            int actually_running = 0;
            for(int i = 1; i < world; i++){
                MPI_Recv(&isRunning, 1, MPI_INT32_T, i, 100 * i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                actually_running |= isRunning;
                // std::cout << "Received isRunning : " << isRunning << std::endl;
                // if (!isRunning) break;
                int comp_rank = i - 1, comp_size = world - 1;
                
                int unit_size = simu_size / comp_size;
                size_t start = unit_size * comp_rank;
                size_t row_sz = params.discretization;
                size_t end = std::min(start + unit_size, simu_size);

                // std::cout << "Waiting for maps from " << i << " tags " << 101 * i << " " << 102 * i << std::endl;
                MPI_Recv(fire_map.data() + start, end - start + 1, MPI_UINT8_T, i, 101 * i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // std::cout << i << " Received fire_map (" << start << ", " << end << ")" << std::endl;
                MPI_Recv(vegetal_map.data() + start, end - start + 1, MPI_UINT8_T, i, 102 * i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // std::cout << i << " Received vegetal_map (" << start << ", " << end << ")" << std::endl;
            }
            isRunning = actually_running;
            for (int i = 1; i < world; i++) {
                MPI_Send(&actually_running, 1, MPI_INT32_T, i, 0, MPI_COMM_WORLD);
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