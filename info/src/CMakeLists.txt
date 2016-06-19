set(SOURCE_FILES "PluginMain.cpp" "PrimitiveGeneratorNode.h" "PrimitiveGeneratorNode.cpp" "PrimitiveGeneratorLoc.h" "PrimitiveGeneratorLoc.h" "PrimitiveGeneratorLoc.cpp" "PrimitiveGeneratorCommand.h" "PrimitiveGeneratorCommand.cpp" "DagUtils.cpp" "PrimitiveGeneratorProfiles.h" "DagUtils.h" "AETemplate.h" "icons.h" "Utils.h"
)

find_package(Maya REQUIRED)

include_directories(${MAYA_INCLUDE_DIR})
link_directories(${MAYA_LIBRARY_DIR})
add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} ${MAYA_LIBRARIES})

MAYA_PLUGIN(${PROJECT_NAME})
