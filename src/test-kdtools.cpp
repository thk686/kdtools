/*
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 *
 * For your own packages, ensure that your test files are
 * placed within the `src/` folder, and that you include
 * `LinkingTo: testthat` within your DESCRIPTION file.
 */

// All test files should include the <testthat.h>
// header file.
#include <testthat.h>

#include <strdist.h>
using namespace keittlab::strdist;

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
context("C++ string distance") {

  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("string distance is correct") {
    expect_true(levenshtein("", "") == 0);
    expect_true(levenshtein("test", "") == 4);
    expect_true(levenshtein("", "testing") == 7);
    expect_true(levenshtein("test", "testing") == 3);
    expect_true(levenshtein("asdfasdlkjadskfjakdjfadkjfa;kldjfa;kjf;alkjf",
                            "a;klja;lkjapioaewoihja;knaa;kjd;kljrea") == 29);
  }

}
