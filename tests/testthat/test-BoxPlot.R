test_that("BoxPlot() returns a ggplot object", {
  data <- Data2[1:20, ] |>
    t() |>
    as.data.frame()
  data <- cbind(ID = rownames(data), data)
  p <- BoxPlot(File = data)
  expect_true(ggplot2::is_ggplot(p))
})
