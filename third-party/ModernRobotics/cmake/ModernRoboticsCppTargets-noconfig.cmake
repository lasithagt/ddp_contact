#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "ModernRoboticsCpp" for configuration ""
set_property(TARGET ModernRoboticsCpp APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(ModernRoboticsCpp PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libModernRoboticsCpp.so.3.4.1"
  IMPORTED_SONAME_NOCONFIG "libModernRoboticsCpp.so.3"
  )

list(APPEND _IMPORT_CHECK_TARGETS ModernRoboticsCpp )
list(APPEND _IMPORT_CHECK_FILES_FOR_ModernRoboticsCpp "${_IMPORT_PREFIX}/lib/libModernRoboticsCpp.so.3.4.1" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
