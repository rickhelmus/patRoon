#include <thread>
#include "thread_mgr.h"
#include "so_manager.h"
#include "bruker_api.h"

std::unique_ptr<ThreadingManager> ThreadingManager::instance;


ThreadingManager::ThreadingManager() :
n_threads(std::thread::hardware_concurrency()),
threading_type(CONVERTER_THREADING)
{}

ThreadingManager::~ThreadingManager() {}

ThreadingManager& ThreadingManager::get_instance()
{
    if(!instance)
        instance = std::make_unique<DefaultThreadingManager>();
    return *instance.get();
}

void ThreadingManager::set_opentims_threading()
{
    threading_type = OPENTIMS_THREADING;
    signal_threading_changed();
}

void ThreadingManager::set_converter_threading()
{
    threading_type = CONVERTER_THREADING;
    signal_threading_changed();
}

void ThreadingManager::set_shared_threading()
{
    threading_type = SHARED_THREADING;
    signal_threading_changed();
}


void ThreadingManager::set_num_threads(size_t n)
{
    if(n == 0)
        n_threads = std::thread::hardware_concurrency();
    else
        n_threads = n;

    signal_threads_changed();
}

/*
 * DefaultThreadingManager
 */

DefaultThreadingManager::~DefaultThreadingManager() {}

size_t DefaultThreadingManager::get_no_opentims_threads()
{
    return n_threads * io_overhead;
}


/*
 * BrukerThreadingManager
 */

BrukerThreadingManager::BrukerThreadingManager(const ThreadingManager& prev_instance, const std::string& bruker_so_path) :
ThreadingManager(prev_instance),
bruker_lib(bruker_so_path),
tims_set_num_threads(bruker_lib.symbol_lookup<tims_set_num_threads_t>("tims_set_num_threads"))
{
    set_bruker_threads();
}

BrukerThreadingManager::~BrukerThreadingManager() {}

void BrukerThreadingManager::SetupBrukerThreading(const std::string& bruker_so_path)
{
    ThreadingManager::instance = std::make_unique<BrukerThreadingManager>(ThreadingManager::get_instance(), bruker_so_path);
}

void BrukerThreadingManager::set_bruker_threads()
{
    switch(threading_type)
    {
        case OPENTIMS_THREADING:
            tims_set_num_threads(1);
            break;
        case SHARED_THREADING:
            tims_set_num_threads(sqrt(n_threads * io_overhead) + 0.5); // Square root, rounded up
            break;
        case CONVERTER_THREADING:
            tims_set_num_threads(n_threads);
            break;
        default:
            throw std::logic_error("Invalid threading model");
    }
}

void BrukerThreadingManager::signal_threading_changed()
{
    set_bruker_threads();
}

void BrukerThreadingManager::signal_threads_changed()
{
    set_bruker_threads();
}

size_t BrukerThreadingManager::get_no_opentims_threads()
{
    switch(threading_type)
    {
        case OPENTIMS_THREADING:
            return n_threads * io_overhead;
        case SHARED_THREADING:
            return sqrt(n_threads * io_overhead) + 0.5;
        case CONVERTER_THREADING:
            return 1;
        default:
            throw std::logic_error("Invalid threading model");
    }
}
