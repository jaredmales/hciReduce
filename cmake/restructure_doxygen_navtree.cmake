# Reshape Doxygen's generated sidebar into hciReduce's documentation hierarchy.
# Doxygen's layout language cannot nest generated indexes (Todo, files, classes,
# and namespaces) below a documentation group, so adjust the emitted navigation
# data after a successful Doxygen run.

if(NOT DEFINED HCIREDUCE_DOC_HTML_DIR)
    message(FATAL_ERROR "HCIREDUCE_DOC_HTML_DIR is required")
endif()

set(_hcireduce_navtree "${HCIREDUCE_DOC_HTML_DIR}/navtreedata.js")
if(NOT EXISTS "${_hcireduce_navtree}")
    message(FATAL_ERROR "Doxygen did not generate ${_hcireduce_navtree}")
endif()

file(READ "${_hcireduce_navtree}" _hcireduce_navtree_contents)
string(FIND "${_hcireduce_navtree_contents}" "var NAVTREE =" _hcireduce_navtree_start)
string(FIND "${_hcireduce_navtree_contents}" "var NAVTREEINDEX =" _hcireduce_navtree_index_start)
if(_hcireduce_navtree_start LESS 0 OR _hcireduce_navtree_index_start LESS 0)
    message(FATAL_ERROR "Unexpected Doxygen navtreedata.js format")
endif()

string(SUBSTRING "${_hcireduce_navtree_contents}" 0 ${_hcireduce_navtree_start} _hcireduce_navtree_header)
string(SUBSTRING "${_hcireduce_navtree_contents}" ${_hcireduce_navtree_index_start} -1 _hcireduce_navtree_index)

set(_hcireduce_navtree_data [=[var NAVTREE =
[
  [ "hciReduce", "index.html", [
    [ "User's Guide", "introduction.html", [
      [ "Introduction", "introduction.html", null ],
      [ "Common Configuration and Workflow", "group__common__user__guide.html", null ],
      [ "klipReduce", "group__klipreduce__user__guide.html", null ],
      [ "p4Reduce", "group__p4reduce__user__guide.html", [
        [ "Introduction", "group__p4reduce__user__guide.html", null ],
        [ "Pixel-Local Processing", "group__p4__local__user__guide.html", null ],
        [ "Negative-Planet Optimization", "group__p4__optimizer__user__guide.html", null ],
        [ "Frozen-Model PSF Responses", "group__p4__psf__user__guide.html", null ],
        [ "Rotated-Frame Regression (Negative Result)", "group__p4__rotated__user__guide.html", null ]
      ] ],
      [ "hciAnalyze", "group__hcianalyze__user__guide.html", null ]
    ] ],
    [ "Developer's Guide", "group__programming__guide.html", [
      [ "Introduction", "group__programming__guide.html", null ],
      [ "KLIP Algorithm", "group__developer__klip__algorithm.html", null ],
      [ "P4 Algorithm", "group__developer__p4__algorithm.html", null ],
      [ "Library API", "group__programming__library.html", null ],
      [ "Testing and Verification", "group__hcireduce__testing.html", [
        [ "Building Tests", "group__testing__building.html", null ],
        [ "Coverage Report", "group__testing__coverage.html", null ],
        [ "Unit Tests", "group__unit__tests.html", [
          [ "HCIobservation Unit Tests", "group__HCIobservation__unit__tests.html", null ],
          [ "HCI Configuration Type Unit Tests", "group__HCI__unit__tests.html", null ],
          [ "ADIobservation Unit Tests", "group__ADIobservation__unit__tests.html", null ],
          [ "KLIPreduction Unit Tests", "group__KLIPreduction__unit__tests.html", null ],
          [ "P4 PCA Unit Tests", "group__P4PCA__unit__tests.html", null ],
          [ "P4 Pixel Grid Unit Tests", "group__P4PixelGrid__unit__tests.html", null ],
          [ "P4 Reduction Unit Tests", "group__P4Reduction__unit__tests.html", null ],
          [ "klipReduce Application Unit Tests", "group__klipReduce__unit__tests.html", null ],
          [ "p4Reduce Application Unit Tests", "group__p4Reduce__unit__tests.html", null ]
        ] ]
      ] ],
      [ "To-Do", "todo.html", null ],
      [ "Namespaces", "namespaces.html", [
        [ "Namespace List", "namespaces.html", "namespaces_dup" ],
        [ "Namespace Members", "namespacemembers.html", [
          [ "All", "namespacemembers.html", null ],
          [ "Functions", "namespacemembers_func.html", null ],
          [ "Enumerations", "namespacemembers_enum.html", null ]
        ] ]
      ] ],
      [ "Classes", "annotated.html", [
        [ "Class List", "annotated.html", "annotated_dup" ],
        [ "Class Index", "classes.html", null ],
        [ "Class Hierarchy", "hierarchy.html", "hierarchy" ],
        [ "Class Members", "functions.html", [
          [ "All", "functions.html", "functions_dup" ],
          [ "Functions", "functions_func.html", null ],
          [ "Variables", "functions_vars.html", null ],
          [ "Typedefs", "functions_type.html", null ]
        ] ]
      ] ],
      [ "Files", "files.html", [
        [ "File List", "files.html", "files_dup" ],
        [ "File Members", "globals.html", [
          [ "All", "globals.html", null ],
          [ "Functions", "globals_func.html", null ]
        ] ]
      ] ]
    ] ],
    [ "References", "citelist.html", null ]
  ] ]
];

]=])

file(WRITE "${_hcireduce_navtree}" "${_hcireduce_navtree_header}${_hcireduce_navtree_data}${_hcireduce_navtree_index}")

# Doxygen's external index files contain paths through its original root-level
# Topics, To-Do, References, Namespaces, Classes, and Files entries. Rewrite
# those paths to match the tree above. Process the root entries first, then
# remove the old Topics level from the two documentation groups.
file(GLOB _hcireduce_navtree_indexes "${HCIREDUCE_DOC_HTML_DIR}/navtreeindex*.js")
foreach(_hcireduce_navtree_index_file IN LISTS _hcireduce_navtree_indexes)
    file(READ "${_hcireduce_navtree_index_file}" _hcireduce_navtree_index_contents)

    string(REGEX REPLACE "(:)\\[1," "\\1[1,3," _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(:)\\[1\\]" "\\1[1,3]" _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(:)\\[2," "\\1[1,5," _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(:)\\[2\\]" "\\1[1,5]" _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(:)\\[3," "\\1[2," _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(:)\\[3\\]" "\\1[2]" _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(:)\\[4," "\\1[1,6," _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(:)\\[4\\]" "\\1[1,6]" _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(:)\\[5," "\\1[1,7," _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(:)\\[5\\]" "\\1[1,7]" _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(:)\\[6," "\\1[1,8," _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(:)\\[6\\]" "\\1[1,8]" _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")

    string(REGEX REPLACE "(:)\\[0,1," "\\1[1," _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(:)\\[0,1\\]" "\\1[1]" _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(:)\\[0,0," "\\1[0," _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(:)\\[0,0\\]" "\\1[0]" _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")

    # The custom User's Guide has a standalone Introduction before Doxygen's
    # four user-guide groups. Preserve each page's selected node after adding
    # that child; otherwise links select the preceding sibling.
    string(REGEX REPLACE "(\\\"introduction\\.html\\\":)\\[[0-9,]+\\]" "\\1[0,0]"
                         _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(\\\"group__common__user__guide\\.html\\\":)\\[0\\]" "\\1[0,1]"
                         _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(\\\"group__klipreduce__user__guide\\.html\\\":)\\[0,1\\]" "\\1[0,2]"
                         _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    # Select P4's Introduction child for the overview URL. The p4Reduce parent
    # links to that same page, so navigating through either entry leaves the
    # P4 subtree expanded in the sidebar.
    string(REGEX REPLACE "(\\\"group__p4reduce__user__guide\\.html\\\":)\\[0,2\\]" "\\1[0,3,0]"
                         _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")

    string(REGEX REPLACE "(\\\"group__p4__local__user__guide\\.html\\\":)\\[[0-9,]+\\]" "\\1[0,3,1]"
                         _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(\\\"group__p4__optimizer__user__guide\\.html\\\":)\\[[0-9,]+\\]" "\\1[0,3,2]"
                         _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(\\\"group__p4__psf__user__guide\\.html\\\":)\\[[0-9,]+\\]" "\\1[0,3,3]"
                         _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(\\\"group__p4__rotated__user__guide\\.html\\\":)\\[[0-9,]+\\]" "\\1[0,3,4]"
                         _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(\\\"group__hcianalyze__user__guide\\.html[^\\\"]*\\\":)\\[[0-9,]+\\]" "\\1[0,4]"
                         _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")

    # The Developer's Guide parent and Introduction child share one URL. Keep
    # the subtree open and select Introduction, just as for p4Reduce.
    string(REGEX REPLACE "(\\\"group__programming__guide\\.html[^\\\"]*\\\":)\\[[0-9,]+\\]" "\\1[1,0]"
                         _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(\\\"group__developer__klip__algorithm\\.html[^\\\"]*\\\":)\\[[0-9,]+\\]" "\\1[1,1]"
                         _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(\\\"group__developer__p4__algorithm\\.html[^\\\"]*\\\":)\\[[0-9,]+\\]" "\\1[1,2]"
                         _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(\\\"group__programming__library\\.html[^\\\"]*\\\":)\\[[0-9,]+\\]" "\\1[1,3]"
                         _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(\\\"group__hcireduce__testing\\.html[^\\\"]*\\\":)\\[[0-9,]+\\]" "\\1[1,4]"
                         _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(\\\"group__testing__building\\.html[^\\\"]*\\\":)\\[[0-9,]+\\]" "\\1[1,4,0]"
                         _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(\\\"group__testing__coverage\\.html[^\\\"]*\\\":)\\[[0-9,]+\\]" "\\1[1,4,1]"
                         _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(\\\"group__unit__tests\\.html[^\\\"]*\\\":)\\[[0-9,]+\\]" "\\1[1,4,2]"
                         _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")

    set(_hcireduce_unit_test_groups
        HCIobservation
        HCI
        ADIobservation
        KLIPreduction
        P4PCA
        P4PixelGrid
        P4Reduction
        klipReduce
        p4Reduce)
    set(_hcireduce_unit_test_index 0)
    foreach(_hcireduce_unit_test_group IN LISTS _hcireduce_unit_test_groups)
        string(REGEX REPLACE
               "(\\\"group__${_hcireduce_unit_test_group}__unit__tests\\.html[^\\\"]*\\\":)\\[[0-9,]+\\]"
               "\\1[1,4,2,${_hcireduce_unit_test_index}]"
               _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
        math(EXPR _hcireduce_unit_test_index "${_hcireduce_unit_test_index} + 1")
    endforeach()

    # These generated indexes must resolve inside the custom tree; otherwise
    # synchronizing them collapses the Developer's Guide subtree.
    string(REGEX REPLACE "(\\\"todo\\.html\\\":)\\[[0-9,]+\\]" "\\1[1,5]"
                         _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(\\\"files\\.html\\\":)\\[[0-9,]+\\]" "\\1[1,8,0]"
                         _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")
    string(REGEX REPLACE "(\\\"citelist\\.html\\\":)\\[[0-9,]+\\]" "\\1[2]"
                         _hcireduce_navtree_index_contents "${_hcireduce_navtree_index_contents}")

    file(WRITE "${_hcireduce_navtree_index_file}" "${_hcireduce_navtree_index_contents}")
endforeach()

# Doxygen names its generated citation page "Bibliography". hciReduce exposes
# that single page as the higher-level References section.
set(_hcireduce_references_page "${HCIREDUCE_DOC_HTML_DIR}/citelist.html")
if(EXISTS "${_hcireduce_references_page}")
    file(READ "${_hcireduce_references_page}" _hcireduce_references_contents)
    string(REPLACE "<title>hciReduce: Bibliography</title>" "<title>hciReduce: References</title>"
                   _hcireduce_references_contents "${_hcireduce_references_contents}")
    string(REPLACE "<div class=\"title\">Bibliography " "<div class=\"title\">References "
                   _hcireduce_references_contents "${_hcireduce_references_contents}")
    file(WRITE "${_hcireduce_references_page}" "${_hcireduce_references_contents}")
endif()
