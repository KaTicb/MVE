cmake_minimum_required(VERSION 3.22)

add_executable(lab2)

target_sources(lab2 PRIVATE main.cxx ${CMAKE_CURRENT_LIST_DIR}/../dataMDD.bin)

add_custom_target(
        build_lab2_plots DEPENDS lab2
        COMMAND python ${CMAKE_CURRENT_LIST_DIR}/draw_dynamic_plot.py
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Running Python script to build plots for lab2"
)