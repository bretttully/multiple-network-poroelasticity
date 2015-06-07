file(GLOB src_files
    "*.h"
    "*.cpp"
)
# this is here to ensure that the IDE knows about this file
# and can then run the code formatter
list(APPEND src_files "${src_dir}/../uncrustify.cfg")

add_library(pympet SHARED ${src_files})
target_link_libraries(pympet boost_python)
install(TARGETS pympet DESTINATION lib)
