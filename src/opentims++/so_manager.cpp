#include "so_manager.h"

#if defined(OPENTIMS_UNIX)

#include <dlfcn.h>

// RAII-style wrapper for results of dlopen()
LoadedLibraryHandle::LoadedLibraryHandle(const std::string& path) : os_handle(nullptr)
{
    os_handle = dlopen(path.c_str(), RTLD_NOW);
    if(os_handle == nullptr)
        throw std::runtime_error(std::string("dlopen(") + path + ") failed, reason: " + dlerror());
}

LoadedLibraryHandle::~LoadedLibraryHandle()
{
    if(os_handle != nullptr)
        dlclose(os_handle);
    // Deliberately not handling errors in dlclose() call here.
}

#elif defined(OPENTIMS_WINDOWS)

#include <libloaderapi.h>
#include <errhandlingapi.h>

LoadedLibraryHandle::LoadedLibraryHandle(const std::string& path) : os_handle(nullptr)
{
    os_handle = LoadLibraryExA(path.c_str(), nullptr, 0);
    if(os_handle == nullptr)
        throw std::runtime_error(std::string("LoadLibraryExA(") + path + ") failed, reason: " + std::to_string(GetLastError()));
}

LoadedLibraryHandle::~LoadedLibraryHandle()
{
    if(os_handle != nullptr)
        FreeLibrary(os_handle);
}

#else

LoadedLibraryHandle::LoadedLibraryHandle(const std::string& path) {};
LoadedLibraryHandle::~LoadedLibraryHandle() {};

#endif
