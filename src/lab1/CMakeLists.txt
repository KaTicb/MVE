cmake_minimum_required(VERSION 3.22)

add_executable(lab1)

target_sources(lab1 PRIVATE main.cxx ${CMAKE_CURRENT_LIST_DIR}/../dataMDD.bin)

add_custom_target(
        build_lab1_plots DEPENDS lab1
        COMMAND python ${CMAKE_CURRENT_LIST_DIR}/draw_static_plot.py
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Running Python script to build plots for lab1"
)