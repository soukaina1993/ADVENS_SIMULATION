cmake_minimum_required(VERSION 3.31)
project(fluids)

#set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
#add_compile_options(-std=c++17)


# Chemin vers les binaires et includes de MSYS2
set(CMAKE_PREFIX_PATH "C:/msys64/ucrt64")

# Détection automatique via pkg-config
find_package(PkgConfig REQUIRED)
#pkg_check_modules(PQ REQUIRED libpq)
pkg_check_modules(PQXX REQUIRED libpqxx)
pkg_check_modules(GDAL REQUIRED gdal)

include_directories(${PQXX_INCLUDE_DIRS} ${GDAL_INCLUDE_DIRS})
link_directories(${PQXX_LIBRARY_DIRS} ${GDAL_LIBRARY_DIRS})



#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
set(CMAKE_CSS_STANDARD_LIBRAIRIES "C:/msys64/ucrt64/bin/libstdc++-6.dll")
set(CMAKE_CXX_STANDARD_LIBRARIES "-static-libgcc -static-libstdc++ -lwsock32 -lws2_32 ${CMAKE_CSS_STANDARD_LIBRARIES}")  #--> gives basic_string error, solved
#set(CMAKE_CXX_STANDARD_LIBRARIES "-static-libgcc -lwsock32 -lws2_32 ${CMAKE_CSS_STANDARD_LIBRARIES}")

#----------------------------
# CoolProp
set(COOLPROP_STATIC_LIBRARY true CACHE BOOL "Use CoolProp")
set(FORCE_BITNESS_64 true)

add_subdirectory ("${CMAKE_SOURCE_DIR}/externals/CoolProp" "CoolProp")
file(GLOB_RECURSE COOL_SOURCES externals/CoolProp/include/*.h externals/CoolProp/src/*.cpp)
add_executable(Prop ${COOL_SOURCES} testCoolProp.cpp)
target_link_libraries (Prop CoolProp)
#------------------------------

# BASIC ENERGY SYSTEMS
file(GLOB_RECURSE BASICES_SOURCES BasicEnergySystem/BasicES.cpp BasicEnergySystem/BasicES.h BasicEnergySystem/BasicES.cc)
file(GLOB_RECURSE EXCHANGER_SOURCES BasicEnergySystem/exchanger/*.cpp BasicEnergySystem/exchanger/*.h BasicEnergySystem/exchanger/*.cc)
file(GLOB_RECURSE TURBINE_SOURCES BasicEnergySystem/turbine/*.cpp BasicEnergySystem/turbine/*.h BasicEnergySystem/turbine/*.cc)
file(GLOB_RECURSE PIPE_SOURCES BasicEnergySystem/pipe/*.cpp BasicEnergySystem/pipe/*.h BasicEnergySystem/pipe/*.cc)
file(GLOB_RECURSE HEATPUMP_SOURCES BasicEnergySystem/heatpump/*.cpp BasicEnergySystem/heatpump/*.h BasicEnergySystem/heatpump/*.cc)

#file(GLOB_RECURSE NETWORK_SOURCES network/*.cpp network/*.h network/*.cc)       # source was modified, do not use!!
file(GLOB_RECURSE QNETWORK_SOURCES network/*.cpp network/*.h network/*.cc)
list(FILTER QNETWORK_SOURCES EXCLUDE REGEX "network/TConnec.*")
list(FILTER QNETWORK_SOURCES EXCLUDE REGEX "network/network.*")

# old versions
file(GLOB_RECURSE HYDRAULIC_SOURCES hydraulic/*.h hydraulic/*.cpp hydraulic/*.cc)
file(GLOB_RECURSE EPANET_SOURCES EPANET/*.h EPANET/*.cpp EPANET/*.cc)
file(GLOB_RECURSE INIT_SOURCES central/*.h central/*.cpp central/*.cc)      # removed

# CONSUMER AND ENVIRONMENT
#file(GLOB_RECURSE CONSUMER_SOURCES consumer/*.h consumer/*.cpp consumer/*.cc)   # source was modified, do not use!!
file(GLOB_RECURSE BUILDING_SOURCES consumer/*Building.h consumer/*Building.cpp consumer/District.h consumer/District.cpp consumer/solar.h consumer/solar.cpp)
file(GLOB_RECURSE ENVIRONMENT_SOURCES environment/*.h environment/*.cpp environment/*.cc)

# FLUIDS AND THERMODYNAMICS
file(GLOB_RECURSE FLUIDS_SOURCES fluids/*.h fluids/*.cpp fluids/*.cc)       # own libraries
list(FILTER FLUIDS_SOURCES EXCLUDE REGEX "fluids/ExtCoolProp.*")
file(GLOB_RECURSE FLUIDSCOOL_SOURCES fluids/ExtCoolProp.h fluids/ExtCoolProp.cpp)   # CoolProp library
#file(GLOB_RECURSE FLUIDAUTRE_SOURCES fluidautre/*.h fluidautre/*.cpp fluidautre/*.cc)      # removed
file(GLOB_RECURSE FLUX_SOURCES Flux/*.h Flux/*.cpp Flux/*.cc)           # removed
file(GLOB_RECURSE FLOW_SOURCES Flow/*.h Flow/*.cpp Flow/*.cc)
file(GLOB_RECURSE STATE_SOURCES state/*.cpp state/*.h state/*.cc)

# HELPERS
file(GLOB_RECURSE IOFILE_SOURCES iofile/*.h iofile/*.cpp iofile/*.cc)
file(GLOB_RECURSE MATH_SOURCES math/*.h math/*.cpp math/*.cc)

# GEOGRAPHY
file(GLOB_RECURSE GEO_SOURCES geodata/*.h geodata/*.cpp geodata/*.cc)

# GENETIC ALGORITHM
file(GLOB_RECURSE GA_SOURCES gao/*.cpp gao/*.h gao/*.cc)

# POSTGRESQL DATABASE
file(GLOB_RECURSE DB_SOURCES db_connect/*.cpp db_connect/*.h db_connect/*.cc)


# Old Files, mainly tests
add_executable(Hydraulic ${FLOW_SOURCES} ${HYDRAULIC_SOURCES} ${STATE_SOURCES} ${FLUIDS_SOURCES} ${BASICES_SOURCES} ${IOFILE_SOURCES} ${EPANET_SOURCES} ${MATH_SOURCES} hydraulicmain.cpp)
add_executable(Hood ${FLOW_SOURCES} ${STATE_SOURCES} ${INIT_SOURCES} ${FLUIDS_SOURCES} ${IOFILE_SOURCES} ${ENVIRONMENT_SOURCES} ${BASICES_SOURCES} ${BUILDING_SOURCES} ${MATH_SOURCES} Hoodmain.cpp)
add_executable(Iofile ${IOFILE_SOURCES} ${FLUIDS_SOURCES} ${MATH_SOURCES} iomain.cpp)
add_executable(Exchanger ${STATE_SOURCES} ${FLOW_SOURCES} ${FLUX_SOURCES} ${EXCHANGER_SOURCES} ${IOFILE_SOURCES} ${BASICES_SOURCES} ${FLUIDS_SOURCES}  ${EPANET_SOURCES} ${MATH_SOURCES} ${GA_SOURCES} BasicESmain.cpp)
add_executable(Newcom ${STATE_SOURCES} ${FLOW_SOURCES} ${FLUX_SOURCES} ${GA_SOURCES} ${TURBINE_SOURCES} ${HEATPUMP_SOURCES} ${BASICES_SOURCES} ${INIT_SOURCES} ${FLUIDS_SOURCES} ${IOFILE_SOURCES} ${ENVIRONMENT_SOURCES} ${BUILDING_SOURCES} ${MATH_SOURCES} ${COMPOSANT_SOURCES} ${EPANET_SOURCES} testcomposant.cpp)
add_executable(expmain ${STATE_SOURCES} ${FLOW_SOURCES} ${FLUX_SOURCES} ${TURBINE_SOURCES} ${HEATPUMP_SOURCES} ${BASICES_SOURCES} ${INIT_SOURCES} ${FLUIDS_SOURCES} ${IOFILE_SOURCES} ${ENVIRONMENT_SOURCES} ${BUILDING_SOURCES} ${MATH_SOURCES} ${COMPOSANT_SOURCES} ${EPANET_SOURCES} ${GA_SOURCES} expmain.cpp)
add_executable(Fluid ${STATE_SOURCES} ${FLOW_SOURCES} ${FLUX_SOURCES} ${FLUIDS_SOURCES} ${IOFILE_SOURCES} ${MATH_SOURCES} fluidsMethod.cpp)
add_executable(hxmain ${STATE_SOURCES} ${FLOW_SOURCES} ${FLUX_SOURCES} ${EXCHANGER_SOURCES} ${BASICES_SOURCES} ${INIT_SOURCES} ${FLUIDS_SOURCES} ${IOFILE_SOURCES} ${ENVIRONMENT_SOURCES} ${BUILDING_SOURCES} ${MATH_SOURCES} ${COMPOSANT_SOURCES} ${HYDRAULIC_SOURCES} ${EPANET_SOURCES} ${GA_SOURCES} exchmain.cpp)
add_executable(Geo ${GEO_SOURCES} testGeo.cpp)

# Quantum network solver test file
add_executable(Testnet ${STATE_SOURCES} ${FLOW_SOURCES} ${GEO_SOURCES} ${GA_SOURCES} ${BASICES_SOURCES} ${PIPE_SOURCES} ${BUILDING_SOURCES} ${ENVIRONMENT_SOURCES} ${FLUIDS_SOURCES} ${FLUIDSCOOL_SOURCES} ${IOFILE_SOURCES} ${MATH_SOURCES} ${HEATPUMP_SOURCES} ${QNETWORK_SOURCES} ${DB_SOURCES} testnetwork.cpp)
target_link_libraries (Testnet CoolProp ${PQXX_LIBRARIES} ${GDAL_LIBRARIES})

# Quantum network solver
add_executable(QNet ${STATE_SOURCES} ${FLOW_SOURCES} ${GEO_SOURCES} ${GA_SOURCES} ${BASICES_SOURCES} ${PIPE_SOURCES} ${BUILDING_SOURCES} ${ENVIRONMENT_SOURCES} ${FLUIDS_SOURCES} ${IOFILE_SOURCES} ${MATH_SOURCES} ${HEATPUMP_SOURCES} ${QNETWORK_SOURCES} ${DB_SOURCES} QNetworkmain.cpp)
#add_executable(QuickCopy quickcopy.cpp)
add_executable(PrepareDatabase ${STATE_SOURCES} ${FLOW_SOURCES} ${GEO_SOURCES} ${GA_SOURCES} ${BASICES_SOURCES} ${PIPE_SOURCES} ${BUILDING_SOURCES} ${ENVIRONMENT_SOURCES} ${FLUIDS_SOURCES} ${IOFILE_SOURCES} ${MATH_SOURCES} ${HEATPUMP_SOURCES} ${QNETWORK_SOURCES} ${DB_SOURCES} PrepareDatabase.cpp)

# ADVENS_Fluid executable for Shiny App
add_executable(ADVENS_Fluid ${STATE_SOURCES} ${FLUIDS_SOURCES} ${FLUIDSCOOL_SOURCES} ${IOFILE_SOURCES} ${MATH_SOURCES} fluids.cpp)
target_link_libraries (ADVENS_Fluid CoolProp)

# Scratch Files
add_executable(Scratch C:/Users/cornelia.blanke/AppData/Roaming/JetBrains/CLion2025.2/scratches/scratch.cpp)
add_executable(MinExample C:/Users/cornelia.blanke/AppData/Roaming/JetBrains/CLion2025.2/scratches/MinimalExample.cpp)
add_executable(PostgreSQLExample ${GEO_SOURCES} C:/Users/cornelia.blanke/AppData/Roaming/JetBrains/CLion2025.2/scratches/PostgreSQLExample.cpp)
#add_executable(Hello C:/Users/cornelia.blanke/AppData/Roaming/JetBrains/CLion2024.3/scratches/hello.cpp)
# -----------------------------------------------------------------

## GDAL and PostgreSQL - needs MSYS2
find_package(GDAL CONFIG REQUIRED)
add_executable(Jamil ${STATE_SOURCES} ${FLOW_SOURCES} ${GEO_SOURCES} ${GA_SOURCES} ${BASICES_SOURCES} ${PIPE_SOURCES} ${BUILDING_SOURCES} ${ENVIRONMENT_SOURCES} ${FLUIDS_SOURCES} ${IOFILE_SOURCES} ${MATH_SOURCES} ${HEATPUMP_SOURCES} ${QNETWORK_SOURCES} testJamil.cpp)
target_link_libraries(Geo GDAL::GDAL)
target_link_libraries(MinExample GDAL::GDAL)

find_package(PostgreSQL REQUIRED)
#target_link_libraries(PostgreSQLExample pqxx pq)
target_link_libraries(QNet ${PQXX_LIBRARIES} ${GDAL_LIBRARIES})
target_link_libraries(PrepareDatabase ${PQXX_LIBRARIES} ${GDAL_LIBRARIES})
target_link_libraries(PostgreSQLExample ${PQXX_LIBRARIES} ${GDAL_LIBRARIES})
target_link_libraries(Jamil ${PQXX_LIBRARIES} ${GDAL_LIBRARIES})
# -----------------------------------------------------------------

## needs mingw-w64
string(APPEND CMAKE_CXX_STANDARD_LIBRARIES " -Wl,-Bstatic,--whole-archive -lwinpthread -Wl,--no-whole-archive")
find_package (Python3 COMPONENTS Interpreter Development)
add_subdirectory(extern/pybind11)
pybind11_add_module(Example example.cpp)

add_executable(PyBind example2.cpp)
target_link_libraries(PyBind PRIVATE pybind11::embed)