#pragma once

#include <cmath>

#include "bruker_api.h"
#include "so_manager.h"

enum OpentimsThreadingType
{
    OPENTIMS_THREADING,
    CONVERTER_THREADING,
    SHARED_THREADING
};

class ThreadingManager
{
 protected:
    static std::unique_ptr<ThreadingManager> instance;
    size_t n_threads;
    const double io_overhead = 1.2;
    OpentimsThreadingType threading_type;

    virtual void signal_threading_changed() = 0;
    virtual void signal_threads_changed() = 0;

 public:
    ThreadingManager();
    ThreadingManager(const ThreadingManager& other) = default;
    virtual ~ThreadingManager();

    static ThreadingManager& get_instance();

    void set_opentims_threading();
    void set_converter_threading();
    void set_shared_threading();

    void set_num_threads(size_t n);
    void set_threading();

    virtual size_t get_no_opentims_threads() = 0;
};

class DefaultThreadingManager final : public ThreadingManager
{
    void signal_threading_changed() override {};
    void signal_threads_changed() override {};
 public:
    DefaultThreadingManager() = default;
    DefaultThreadingManager(const DefaultThreadingManager& other) = default;
    virtual ~DefaultThreadingManager();

    size_t get_no_opentims_threads() override;
};

class BrukerThreadingManager final : public ThreadingManager
{
    const LoadedLibraryHandle bruker_lib;
    tims_set_num_threads_t* const tims_set_num_threads;
    void set_bruker_threads();

    void signal_threading_changed() override;
    void signal_threads_changed() override;
 public:
    BrukerThreadingManager() = delete;
    BrukerThreadingManager(const ThreadingManager& prev_instance, const std::string& bruker_so_path);
    virtual ~BrukerThreadingManager();

    static void SetupBrukerThreading(const std::string& so_path);

    size_t get_no_opentims_threads() override;
};
