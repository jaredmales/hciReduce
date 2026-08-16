# Updates to Documentation

I want to restructure the hciReduce documentation.  

    1. remove mxlib derived grouping: The doxygen documentation currently has mxlib derived groups, e.g. "Image Processing Files", "High Contrast Imaging.  These aren't need here, so we can have a flatter structure.

    2. The TOC goes hciReduce->User's Guide -> hciReduce.  The inner hciReduce shouldnt' be there

    3. We should have a "User's Guide" for klipReduce and p4Reduce, and a "Programming Guide" that documents the library, etc.

    4. Under the "User's Guide", we want to document the interfaces and "how it works" for the two apps for the scientist, not the programmer.  We should have a section for common configuration and functionality, and then sections for the two algorithms.

    5.  The "Programming Guide" should have an overview intro / main page that describes the class hierarchy and flow.