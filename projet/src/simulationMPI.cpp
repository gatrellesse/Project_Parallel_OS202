#include <string>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <thread>
#include <chrono>
#include <stdio.h>
#include <mpi.h>
#include <numeric>  // For std::accumulate

#include "model.hpp"
#include "display.hpp"

using namespace std::string_literals;
using namespace std::chrono_literals;

struct ParamsType
{
    double length{1.};
    unsigned discretization{20u};
    std::array<double,2> wind{0.,0.};
    Model::LexicoIndices start{10u,10u};
};

void analyze_arg( int nargs, char* args[], ParamsType& params )
{
    if (nargs ==0) return;
    std::string key(args[0]);
    if (key == "-l"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une valeur pour la longueur du terrain !" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.length = std::stoul(args[1]);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    auto pos = key.find("--longueur=");
    if (pos < key.size())
    {
        auto subkey = std::string(key,pos+11);
        params.length = std::stoul(subkey);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }

    if (key == "-n"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une valeur pour le nombre de cases par direction pour la discrétisation du terrain !" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.discretization = std::stoul(args[1]);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    pos = key.find("--number_of_cases=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos+18);
        params.discretization = std::stoul(subkey);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }

    if (key == "-w"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une paire de valeurs pour la direction du vent !" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values =std::string(args[1]);
        params.wind[0] = std::stod(values);
        auto pos = values.find(",");
        if (pos == values.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(values, pos+1);
        params.wind[1] = std::stod(second_value);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    pos = key.find("--wind=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos+7);
        params.wind[0] = std::stoul(subkey);
        auto pos = subkey.find(",");
        if (pos == subkey.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(subkey, pos+1);
        params.wind[1] = std::stod(second_value);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }

    if (key == "-s"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une paire de valeurs pour la position du foyer initial !" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values =std::string(args[1]);
        params.start.column = std::stod(values);
        auto pos = values.find(",");
        if (pos == values.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la position du foyer initial" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(values, pos+1);
        params.start.row = std::stod(second_value);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    pos = key.find("--start=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos+8);
        params.start.column = std::stoul(subkey);
        auto pos = subkey.find(",");
        if (pos == subkey.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(subkey, pos+1);
        params.start.row = std::stod(second_value);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }
}

ParamsType parse_arguments( int nargs, char* args[] )
{
    if (nargs == 0) return {};
    if ( (std::string(args[0]) == "--help"s) || (std::string(args[0]) == "-h") )
    {
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

bool check_params(ParamsType& params)
{
    bool flag = true;
    if (params.length <= 0)
    {
        std::cerr << "[ERREUR FATALE] La longueur du terrain doit être positive et non nulle !" << std::endl;
        flag = false;
    }

    if (params.discretization <= 0)
    {
        std::cerr << "[ERREUR FATALE] Le nombre de cellules par direction doit être positive et non nulle !" << std::endl;
        flag = false;
    }

    if ( (params.start.row >= params.discretization) || (params.start.column >= params.discretization) )
    {
        std::cerr << "[ERREUR FATALE] Mauvais indices pour la position initiale du foyer" << std::endl;
        flag = false;
    }
    
    return flag;
}

void display_params(ParamsType const& params)
{
    std::cout << "Parametres définis pour la simulation : \n"
              << "\tTaille du terrain : " << params.length << std::endl 
              << "\tNombre de cellules par direction : " << params.discretization << std::endl 
              << "\tVecteur vitesse : [" << params.wind[0] << ", " << params.wind[1] << "]" << std::endl
              << "\tPosition initiale du foyer (col, ligne) : " << params.start.column << ", " << params.start.row << std::endl;
}

int main( int nargs, char* args[] )
{
    auto params = parse_arguments(nargs-1, &args[1]);
    display_params(params);
    if (!check_params(params)) return EXIT_FAILURE;

    int rank, world;
    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &world);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::shared_ptr<Displayer> displayer = nullptr;
    
    if (rank == 0){
        displayer = Displayer::init_instance( params.discretization, params.discretization );
    }
    // std::unique_ptr<Model> simu = nullptr;
    
    // Vectors to store results for each run
    int num_runs = 100;
    std::vector<double> time_updates(num_runs, 0.0);
    std::vector<double> time_displays(num_runs, 0.0);
    std::vector<double> avg_time_display(num_runs, 0.0);
    std::vector<double> avg_time_update(num_runs, 0.0);
    std::vector<double> global_times(num_runs, 0.0);
    // auto simu = Model( params.length, params.discretization, params.wind,
    //                    params.start);
    for(int run = 0; run < num_runs; ++run){
        SDL_Event event;
        int n_iterations = 0;
        int maps_sizes = params.discretization * params.discretization;
        double time_display = 0;
        double time_update = 0 ;
        double time_global = 0;
        bool isRunning = true;
        std::pair<bool,double> result;
        MPI_Request send_request1, send_request2, send_request3;
        MPI_Request recv_request1, recv_request2, recv_request3;
        std::unique_ptr<Model> simu = nullptr;
        if (rank != 0) {
            simu = std::make_unique<Model>(params.length, params.discretization, params.wind, params.start);
    
        }
        
        while (isRunning)
        {
            auto start_global_time = std::chrono::high_resolution_clock::now();
            n_iterations += 1;
            if (rank != 0) {
                result = simu->update();
                isRunning = result.first;
                MPI_Isend(&isRunning, 1,  MPI_CXX_BOOL, 0, 101, MPI_COMM_WORLD, &send_request1);
                MPI_Isend(simu->vegetal_map().data(), maps_sizes, MPI_BYTE, 0, 102, MPI_COMM_WORLD, &send_request2);
                MPI_Isend(simu->fire_map().data(), maps_sizes, MPI_BYTE, 0, 103, MPI_COMM_WORLD, &send_request3);
                double elapsed_time = result.second;
                time_update += elapsed_time;
                if ((simu->time_step() & 31) == 0){
                    std::cout << "Time step " << simu->time_step() << "\n===============" << std::endl;
                }
            }

            
            if (rank == 0) {
                std::vector<std::uint8_t> m_vegetal_map(maps_sizes);
                std::vector<std::uint8_t> m_fire_map(maps_sizes);
                
                // Receive isRunning (bool) from rank 1
                MPI_Irecv(&isRunning, 1, MPI_CXX_BOOL, 1, 101, MPI_COMM_WORLD, &recv_request1);
                MPI_Irecv(m_vegetal_map.data(), maps_sizes, MPI_BYTE, 1, 102, MPI_COMM_WORLD, &recv_request2);
                MPI_Irecv(m_fire_map.data(), maps_sizes, MPI_BYTE, 1, 103, MPI_COMM_WORLD, &recv_request3) ;
                // Wait for all receives to complete
                MPI_Wait(&recv_request1, MPI_STATUS_IGNORE);
                MPI_Wait(&recv_request2, MPI_STATUS_IGNORE);
                MPI_Wait(&recv_request3, MPI_STATUS_IGNORE);
                auto start_time = std::chrono::high_resolution_clock::now();
                displayer->update(m_vegetal_map, m_fire_map);
                auto end_time = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed_seq = end_time - start_time;
                time_display += elapsed_seq.count();

                elapsed_seq = end_time - start_global_time;
                time_global += elapsed_seq.count();
            }

            if (SDL_PollEvent(&event) && event.type == SDL_QUIT)
                break;
            //std::this_thread::sleep_for(0.1s);
        }
        if (rank == 0){
            // std::cout << "n_iterations " << n_iterations << std::endl;
            // std::cout << "time_affichage " << time_display << std::endl;
            // std::cout << "Avg time affichage " << time_display / n_iterations << " seconds" << std::endl;
            time_displays[run] = time_display;
            avg_time_display[run] = time_display/n_iterations;
            global_times[run] = time_global;
        }
        else{
            // std::cout << "time_update " << time_update << std::endl;
            // std::cout << "Avg time update " << time_update / n_iterations << " seconds" << std::endl;    
            time_updates[run] = time_update;
            avg_time_update[run] = time_update/n_iterations;
        }
    }
    
    double ovr_time_display = std::accumulate(time_displays.begin(), time_displays.end(), 0.0) / num_runs;
    double ovr_time_updates = std::accumulate(time_updates.begin(), time_updates.end(), 0.0) / num_runs;
    double ovr_avg_time_display = std::accumulate(avg_time_display.begin(), avg_time_display.end(), 0.0) / num_runs;
    double ovr_avg_time_updates = std::accumulate(avg_time_update.begin(), avg_time_update.end(), 0.0) / num_runs;
    double ovr_time_global = std::accumulate(global_times.begin(), global_times.end(), 0.0) / num_runs;
    if(rank==0){
        std::cout << "Overall global time " << ovr_time_global << std::endl;
    }
    if(rank != 0){
        std::cout << "time_update " << ovr_time_updates << std::endl;
        std::cout << "Avg time update " << ovr_avg_time_updates << " seconds" << std::endl;
    }
    else{
        std::cout << "time_affichage " << ovr_time_display << std::endl;
        std::cout << "Avg time affichage " << ovr_avg_time_display << " seconds" << std::endl;        
    }
    MPI_Finalize();
    return EXIT_SUCCESS;
}