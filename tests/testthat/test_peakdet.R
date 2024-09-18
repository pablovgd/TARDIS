test_that(".find_peak_border works", {
    ints <- c(3, 4, 5, 20, 10, 8, 3, 4, 5)
    sc <- c(FALSE, diff(diff(ints) > 0) != 0, FALSE)

    res <- .find_peak_border(sc, 4)
    expect_true(is.numeric(res))
    expect_true(length(res) == 2)
    expect_equal(res, c(left = 1, right = 7))

    ## first right border is too close
    ints <- c(3, 4, 5, 20, 10, 12, 20, 40)
    sc <- c(FALSE, diff(diff(ints) > 0) != 0, FALSE)
    res <- .find_peak_border(sc, 4)
    expect_equal(res, c(left = 1, right = 8))
    ## accepting a sign change in the next point.
    res <- .find_peak_border(sc, 4, min_dist = 1)
    expect_equal(res, c(left = 1, right = 5))

    ## firts right border is exactly at the minimum required distance
    ints <- c(3, 4, 5, 20, 10, 8, 20, 40)
    sc <- c(FALSE, diff(diff(ints) > 0) != 0, FALSE)
    res <- .find_peak_border(sc, 4)
    expect_equal(res, c(left = 1, right = 6))

    ## checking left border
    ints <- c(10, 5, 3, 4, 5, 10, 20, 10, 8, 3, 5)
    sc <- c(FALSE, diff(diff(ints) > 0) != 0, FALSE)
    res <- .find_peak_border(sc, 7)
    expect_equal(res, c(left = 3, right = 10))

    ## left sign change at limit
    ints <- c(10, 5, 5, 10, 20, 10, 8, 3, 5)
    sc <- c(FALSE, diff(diff(ints) > 0) != 0, FALSE)
    res <- .find_peak_border(sc, 5)
    expect_equal(res, c(left = 3, right = 8))

    ## left sign change too closelimit
    ints <- c(10, 5, 5, 20, 10, 8, 3, 5)
    sc <- c(FALSE, diff(diff(ints) > 0) != 0, FALSE)
    res <- .find_peak_border(sc, 4)
    expect_equal(res, c(left = 1, right = 7))
})

test_that("find_peak_points works", {
    ## Single peak
    ints <- c(3, 3, 4, 5, 8, 10, 20, 8, 7, 1, 1)
    res <- find_peak_points(seq_along(ints), ints, targetRtime = 6)
    expect_true(is.numeric(res))
    expect_true(length(res) == 3)
    expect_equal(names(res), c("left", "right", "peak_index"))
    expect_equal(res, c(left = 2, right = 11, peak_index = 7))
    
    ## multiple peaks
    ints <- c(4, 3, 5, 8, 7, 4, 1, 2, 1, 3, 1, 6, 10, 12, 11, 9, 7, 5, 1)
    ## first peak
    res <- find_peak_points(seq_along(ints), ints, targetRtime = 6)
    expect_equal(res, c(left = 2, right = 7, peak_index = 4))
    ## second peak
    res <- find_peak_points(seq_along(ints), ints, targetRtime = 15)
    expect_equal(res, c(left = 11, right = 19, peak_index = 14))
    ## 50%
    ints <- c(4, 3, 5, 5.2, 5.5, 4, 1, 2, 1, 3, 1, 6, 10, 12, 11, 9, 7, 5, 1)
    ## asking for first peak, getting second
    res <- find_peak_points(seq_along(ints), ints, targetRtime = 6)
    expect_equal(res, c(left = 11, right = 19, peak_index = 14))
    
    ## NA in input
    ints <- c(4, 3, 5, 8, 7, 4, 1, 2, 1, 3, 1, 6, 10, 12, 11, 9, 7, 5, 1)
    ints[2] <- NA
    res <- find_peak_points(seq_along(ints), ints, targetRtime = 6)
    expect_equal(res, c(left = 1, right = 7, peak_index = 4))
    ## NA close to apex
    ints[2] <- 3
    ints[13] <- NA
    res <- find_peak_points(seq_along(ints), ints, targetRtime = 15)
    expect_equal(res, c(left = 11, right = 19, peak_index = 14))
})

################################################################################
##             ORIGINAL CODE
##
test_that("find_true_occurence works", {
    ## Unit tests are the same as with new code, with the exception that
    ## the closest acceptable minimum is 3 points/index away
    ints <- c(3, 4, 5, 20, 10, 8, 3, 4, 5)
    sc <- c(FALSE, diff(diff(ints) > 0) != 0, FALSE)

    res <- find_true_occurrence(sc, ints, 4)
    expect_true(is.numeric(res))
    expect_true(length(res) == 2)
    expect_equal(res, c(left_true = 1, right_true = 7))

    ## first right border is too close
    ints <- c(3, 4, 5, 20, 10, 12, 20, 40)
    sc <- c(FALSE, diff(diff(ints) > 0) != 0, FALSE)
    res <- find_true_occurrence(sc, ints, 4)
    expect_equal(res, c(left_true = 1, right_true = 8))

    ## first right border is exactly at the minimum required distance
    ints <- c(3, 4, 5, 20, 10, 8, 5, 20, 40)
    sc <- c(FALSE, diff(diff(ints) > 0) != 0, FALSE)
    res <- find_true_occurrence(sc, ints, 4)
    expect_equal(res, c(left_true = 1, right_true = 7))

    ## checking left border
    ints <- c(10, 5, 3, 4, 5, 10, 20, 10, 8, 3, 5)
    sc <- c(FALSE, diff(diff(ints) > 0) != 0, FALSE)
    res <- find_true_occurrence(sc, ints, 7)
    expect_equal(res, c(left_true = 3, right_true = 10))

    ## left sign change at limit
    ints <- c(10, 3, 5, 5, 10, 20, 10, 8, 3, 5)
    sc <- c(FALSE, diff(diff(ints) > 0) != 0, FALSE)
    res <- find_true_occurrence(sc, ints, 5)
    expect_equal(res, c(left_true = 2, right_true = 9))

    ## left sign change too closelimit
    ints <- c(10, 5, 5, 20, 10, 8, 3, 5)
    sc <- c(FALSE, diff(diff(ints) > 0) != 0, FALSE)
    res <- find_true_occurrence(sc, ints, 4)
    expect_equal(res, c(left_true = 1, right_true = 7))
})

test_that("find_peak_points_orig works", {
    ## Single peak
    ints <- c(3, 3, 4, 5, 8, 10, 20, 8, 7, 1, 1)
    res <- find_peak_points_orig(seq_along(ints), ints, 6)
    expect_true(is.numeric(res))
    expect_true(length(res) == 3)
    expect_equal(names(res), c("left_true", "right_true", "peakindex"))
    expect_equal(res, c(left_true = 2, right_true = 11, peakindex = 7))
    res_b <- find_peak_points(seq_along(ints), ints, 6)
    expect_equal(unname(res), unname(res_b))
    
    ## multiple peaks
    ints <- c(4, 3, 4, 5, 8, 7, 4, 1, 2, 1, 3, 1, 6, 10, 12, 11, 9, 7, 5, 1)
    ## first peak
    res <- find_peak_points_orig(seq_along(ints), ints, 6)
    expect_equal(res, c(left_true = 2, right_true = 8, peakindex = 5))
    res_b <- find_peak_points(seq_along(ints), ints, 6)
    expect_equal(unname(res), unname(res_b))
    ## second peak
    res <- find_peak_points_orig(seq_along(ints), ints, 15)
    expect_equal(res, c(left_true = 12, right_true = 20, peakindex = 15))
    res_b <- find_peak_points(seq_along(ints), ints, 15)
    expect_equal(unname(res), unname(res_b))
    ## 50%
    ints <- c(4, 3, 4, 5, 5.2, 5.5, 4, 1, 2, 1, 3, 1, 6, 10, 12, 11, 9, 7, 5, 1)
    ## asking for first peak, getting second
    res <- find_peak_points_orig(seq_along(ints), ints, 6)
    expect_equal(res, c(left_true = 12, right_true = 20, peakindex = 15))
    res_b <- find_peak_points(seq_along(ints), ints, 6)
    expect_equal(unname(res), unname(res_b))
    
    ## NA in input
    ints <- c(4, 3, 5, 8, 7, 4, 1, 2, 1, 3, 1, 6, 10, 12, 11, 9, 7, 5, 1)
    ints[2] <- NA
    res <- find_peak_points_orig(seq_along(ints), ints, 6)
    expect_equal(res, c(left_true = 1, right_true = 7, peakindex = 4))
    ## NA close to apex
    ints[2] <- 3
    ints[13] <- NA
    res <- find_peak_points_orig(seq_along(ints), ints, 15)
    expect_equal(res, c(left_true = 11, right_true = 19, peakindex = 14))
    res_b <- find_peak_points(seq_along(ints), ints, 15)
    expect_equal(unname(res), unname(res_b))
})

