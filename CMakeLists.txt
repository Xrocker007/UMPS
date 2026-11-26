cmake_minimum_required(VERSION 3.21)

project(UMPS LANGUAGES C CXX)

set(VEGA_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/VegaFEM"
    CACHE PATH "Location of the VegaFEM source")

# Build Vega as a library without its example applications.
set(VEGAFEM_BUILD_COPYLEFT OFF CACHE BOOL "" FORCE)
set(VEGAFEM_BUILD_UTILITIES OFF CACHE BOOL "" FORCE)
set(VEGAFEM_BUILD_UTILITIES_GRAPHICS OFF CACHE BOOL "" FORCE)

add_subdirectory("${VEGA_SOURCE_DIR}" "${CMAKE_BINARY_DIR}/vega")

add_executable(UMPS main.cpp plasticForceModel.cpp)
target_compile_features(UMPS PRIVATE cxx_std_17)
target_link_libraries(UMPS PRIVATE vegafem::vegafem)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

target_link_libraries(UMPS PRIVATE OpenGL::GL GLUT::GLUT)

# Copy runtime dependency DLLs beside our executable.
if(WIN32)
    add_custom_command(TARGET UMPS POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_RUNTIME_DLLS:UMPS>
            $<TARGET_FILE_DIR:UMPS>
        COMMAND_EXPAND_LISTS
    )
endif()