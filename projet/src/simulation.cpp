#include <string>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <thread>
#include <chrono>
#include <fstream>
#include <string>
#include <algorithm>
#include <omp.h>
#include <iomanip>
#include <mpi.h>

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
    std::string version{"1"};
};

void analyze_arg( int nargs, char* args[], ParamsType& params )
{
    if (nargs == 0) return;
    std::string key(args[0]);
    if (key == "-v"s) {
        if (nargs < 2)
        {
            std::cerr << "Manque une valeur pour la version du code !" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.version.assign(args[1]);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
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

int main( int nargs, char* args[] ) {
    auto params = parse_arguments(nargs-1, &args[1]);
    //display_params(params);
    if (!check_params(params)) return EXIT_FAILURE;

    int rank, world;
    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &world);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::shared_ptr<Displayer> displayer = nullptr;
    if (rank == 0) displayer = Displayer::init_instance(params.discretization, params.discretization);

    uint32_t maps_sizes = params.discretization * params.discretization;
    float time_update = 0;
    float time_for = 0;
    double time_affichage = 0;
    int n_iterations = 0;

    MPI_Request send_request1, send_request2, send_request3;
    //MPI_Request recv_request1, recv_request2, recv_request3;
    bool isRunning = true;

    std::unique_ptr<Model> simu = nullptr;
    if (rank == 1) simu = std::make_unique<Model>(params.length, params.discretization, params.wind, params.start);

    SDL_Event event;
    double global_start = omp_get_wtime();
    while (isRunning)
    {
        n_iterations++;  
        if(rank == 1) {
            std::vector<std::uint8_t> b_vegetal_map = simu->vegetal_map();
            std::vector<std::uint8_t> b_fire_map = simu->fire_map();
            MPI_Isend(&isRunning, 1,  MPI_CXX_BOOL, 0, 101, MPI_COMM_WORLD, &send_request1);
            MPI_Isend(b_vegetal_map.data(), maps_sizes, MPI_UINT8_T, 0, 102, MPI_COMM_WORLD, &send_request2);
            MPI_Isend(b_fire_map.data(), maps_sizes, MPI_UINT8_T, 0, 103, MPI_COMM_WORLD, &send_request3);
            isRunning = simu->update(&time_update, &time_for);
            int ready1 = 0, ready2 = 0, ready3 = 0;
            while (!ready1 || !ready2 || !ready3) {
                MPI_Test(&send_request1, &ready1, MPI_STATUS_IGNORE);
                MPI_Test(&send_request2, &ready2, MPI_STATUS_IGNORE);
                MPI_Test(&send_request3, &ready3, MPI_STATUS_IGNORE);
            }
        }
        if(rank==0){
            std::vector<std::uint8_t> m_vegetal_map(maps_sizes);
            std::vector<std::uint8_t> m_fire_map(maps_sizes);
            // Receive isRunning (bool) from rank 1
            MPI_Recv(&isRunning, 1, MPI_CXX_BOOL, 1, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(m_vegetal_map.data(), maps_sizes, MPI_UINT8_T, 1, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(m_fire_map.data(), maps_sizes, MPI_UINT8_T, 1, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE) ;
            double start_affichage = omp_get_wtime();
            displayer->update(m_vegetal_map, m_fire_map);
            time_affichage += omp_get_wtime() - start_affichage; 
        }
        if (SDL_PollEvent(&event) && event.type == SDL_QUIT) break;
    }
    if(rank == 1) {
        MPI_Isend(&isRunning, 1,  MPI_CXX_BOOL, 0, 101, MPI_COMM_WORLD, &send_request1);
        MPI_Isend(simu->vegetal_map().data(), maps_sizes, MPI_UINT8_T, 0, 102, MPI_COMM_WORLD, &send_request2);
        MPI_Isend(simu->fire_map().data(), maps_sizes, MPI_UINT8_T, 0, 103, MPI_COMM_WORLD, &send_request3);
    }
    if(rank == 0) std::cout << "global_time " << omp_get_wtime() - global_start << std::endl;
    if(rank == 1) {
        std::ofstream out_file;
        out_file.open("example-"+params.version+".txt");
        for( auto keys : simu->keys_by_step()) {
            for( auto element : keys ) out_file << " " << (int) element ;
            out_file << std::endl;
        }
        std::ofstream update_time_file;
        update_time_file.open("update-time-"+params.version+".txt");
        for( auto element : simu->update_time_by_step()) update_time_file << " " << element ;
        update_time_file << std::endl;
    }
    if(rank == 0) {
        std::ofstream display_time_file;
        display_time_file.open("display-time-"+params.version+".txt");
        for( auto element : displayer->display_time_by_step()) display_time_file << " " << element ;
        display_time_file << std::endl;
    }
    if(rank == 0) std::cout << "n_iterations " << n_iterations << std::endl;
    std::cout << std::fixed;
    if(rank == 1) std::cout << std::setprecision(10) << "avg_time_update " << time_update / n_iterations << std::endl;
    if(rank == 0) std::cout << std::setprecision(10) << "avg_time_affichage " << time_affichage / n_iterations << std::endl;
    if(rank == 0) SDL_Quit();
    MPI_Finalize();
    return EXIT_SUCCESS;
}
