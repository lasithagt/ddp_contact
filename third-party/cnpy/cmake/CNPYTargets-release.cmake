#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "cnpy" for configuration "Release"
set_property(TARGET cnpy APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(cnpy PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libcnpy.so.3.4.1"
  IMPORTED_SONAME_RELEASE "libcnpy.so.3"
  )

list(APPEND _IMPORT_CHECK_TARGETS cnpy )
list(APPEND _IMPORT_CHECK_FILES_FOR_cnpy "${_IMPORT_PREFIX}/lib/libcnpy.so.3.4.1" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
