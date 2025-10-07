cmake_minimum_required(VERSION 3.27.4 FATAL_ERROR)

project(thomson)

find_package(ROOT REQUIRED)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if(${CMAKE_SYSTEM} MATCHES "Linux")
    if (${CMAKE_SYSTEM} MATCHES "generic")
        set(DAS_PATH "/home/user/Documents/VS code/TSignal/build/")
    else()
        set(DAS_PATH /home/korzh/Work/lib/TSignal/lib/)
    endif()
    else()
        message(FATAL_ERROR "Not run Linux")
endif()
    
    
file(GLOB LIBS ${DAS_PATH}/*.so)


include(${ROOT_USE_FILE})
set(GUI_HEADERS ${PROJECT_SOURCE_DIR}/include/ThomsonGUI.h)
set(LINK ${PROJECT_SOURCE_DIR}/LinkDef/LinkDef.h)

set(LIB_NAME ${PROJECT_NAME}GUI)

file(GLOB SRC ${PROJECT_SOURCE_DIR}/src/*.cpp)

root_generate_dictionary(
    G_${LIB_NAME} 
    ${GUI_HEADERS} 
    MODULE ${LIB_NAME} 
    LINKDEF ${LINK}
    OPTIONS "-p"
    NO_INSTALL
)

add_library(${LIB_NAME} SHARED ${SRC} ${PROJECT_BINARY_DIR}/G_${LIB_NAME}.cxx)
target_include_directories(${LIB_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(${LIB_NAME} PUBLIC
    ROOT::Core
    ROOT::Gui
    ROOT::Gpad
    ROOT::RIO
    ROOT::Tree
    ROOT::Hist
    ROOT::Physics
    ${LIBS}
)


set(MAIN_SRC ${PROJECT_SOURCE_DIR}/src/main/main.cpp)
add_executable(${PROJECT_NAME} ${MAIN_SRC})
target_link_libraries(${PROJECT_NAME} PRIVATE ${LIB_NAME})