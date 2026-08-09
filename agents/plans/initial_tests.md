# Tests needed for hciReduce

We need to create a test system for hciReduce.  We'll follow the general flow of the mxlib test system, including use of catch2.  See AGENTS.md for some guidance.

In parallel let's start a doxygen system, also following mxlib.  Docs should look just like the mxlib docs in appearance. We should include the same provisions for including test coverage output.

After getting organized and setting up the system, start with HCIobservation.hpp/cpp.  Propose a suite of tests.  We want to cover both logistics (i.e. does the file load the way it's supposed to, etc) and physics (do we apply the filter we think we do).

Take this in two parts: first do the organization, then document the proposed test suite below under "#plan".  Do not alter this prompt.

Note: a code review is in-scope, and this includes any identified issues with mxlib.

# plan
