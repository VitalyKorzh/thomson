cmake_minimum_required(VERSION 3.27.4 FATAL_ERROR)

project(thomson)

find_package(ROOT REQUIRED COMPONENTS RIO Net Gui)
find_package(Optimizer REQUIRED)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -fPIC")

if(${CMAKE_SYSTEM} MATCHES "Linux")
    if (${CMAKE_SYSTEM} MATCHES "generic")
        #set(DAS_PATH "/home/user/Documents/VS code/TSignal/build/")
    else()
        set(DAS_PATH /home/korzh/Work/lib/TSignal/lib/)
    endif()
    else()
        message(FATAL_ERROR "Not run Linux")
endif()
    
    
file(GLOB DAS_LIBS ${DAS_PATH}/*.so)


include(${ROOT_USE_FILE})
set(GUI_HEADERS ${PROJECT_SOURCE_DIR}/include/ThomsonGUI.h)
set(LINK ${PROJECT_SOURCE_DIR}/LinkDef/LinkDef.h)
set(LIB_NAME ${PROJECT_NAME}GUI)
file(GLOB SRC ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB SRC_GUI ${PROJECT_SOURCE_DIR}/src/GUI/*.cpp)
file(GLOB SRC_THOMSON_COUNTER ${PROJECT_SOURCE_DIR}/src/thomsonCounter/*.cpp)

root_generate_dictionary(
    G_${LIB_NAME} 
    ${GUI_HEADERS} 
    MODULE ${LIB_NAME} 
    LINKDEF ${LINK}
    OPTIONS "-p"
)

add_library(${LIB_NAME} SHARED ${SRC_GUI} ${SRC_THOMSON_COUNTER} ${PROJECT_BINARY_DIR}/G_${LIB_NAME}.cxx)
target_include_directories(${LIB_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(${LIB_NAME} PUBLIC
    ROOT::Core
    ROOT::Gui
    ROOT::Gpad
    ROOT::RIO
    ROOT::Tree
    ROOT::Hist
    ROOT::Physics
    Optimizer::Optimizer
    ${ROOT_LIBRARIES}
    ${DAS_LIBS}
)

add_executable(${PROJECT_NAME} ${SRC})
target_link_libraries(${PROJECT_NAME} PRIVATE ${LIB_NAME})