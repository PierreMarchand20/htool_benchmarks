foreach(compression ${compressions})
    foreach(type ${types})
        foreach(clustering ${clusterings})
            foreach(admissible_condition ${admissible_conditions})
                foreach(kernel ${kernels})

                    add_executable(${type}_${compression}_${clustering}_${admissible_condition}_${kernel} ${type}_${compression}_${clustering}_${admissible_condition}_${kernel}.cpp)
                    target_link_libraries(${type}_${compression}_${clustering}_${admissible_condition}_${kernel} htool_benchmark)

                endforeach()
            endforeach()
        endforeach()
    endforeach()
endforeach()
