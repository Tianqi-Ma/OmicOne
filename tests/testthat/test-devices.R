# Graphics-device hygiene (fct_plots.R): a render never leaves a device open,
# never draws on the wrong device, and never fails because the plot box is
# tiny ("figure margins too large") -- the bugs behind "too many open devices".

test_that("with_null_device closes its device and anything opened inside it", {
  grDevices::graphics.off()
  before <- grDevices::dev.list()
  out <- with_null_device({
    grDevices::pdf(NULL)                   # a stray device left open by the code inside
    42
  })
  expect_equal(out, 42)
  expect_equal(grDevices::dev.list(), before)
  # the previously current device stays current
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  d <- grDevices::dev.cur()
  with_null_device(grDevices::pdf(NULL))
  expect_equal(grDevices::dev.cur(), d)
  grDevices::dev.off(d)
})

test_that("close_devices_since keeps the device it is told to keep", {
  grDevices::graphics.off()
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  keep <- grDevices::dev.cur()
  before <- grDevices::dev.list()
  grDevices::pdf(NULL)
  grDevices::pdf(NULL)
  close_devices_since(before, keep = keep)
  expect_equal(grDevices::dev.list(), before)
  expect_equal(grDevices::dev.cur(), keep)
  grDevices::dev.off(keep)
})

test_that("building composed figures with no device open leaves no device behind", {
  grDevices::graphics.off()
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
  g <- compose_grid(list(p, p), 1, 2, legends = list(p))
  expect_s3_class(g, "omicone_grid")
  expect_null(grDevices::dev.list())       # no default device (Rplots.pdf) opened
})

test_that("an error message draws on a tiny device (grid, no base-graphics margins)", {
  f <- tempfile(fileext = ".png")
  grDevices::png(f, width = 12, height = 12)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(draw_plot_message("Plot error:\nsomething"))
})

test_that("a step plot renders without leaking devices, even when its code opens one", {
  grDevices::graphics.off()
  server <- function(input, output, session) {
    render_step_plot(output, input, "p", function() {
      grDevices::pdf(NULL)                 # misbehaving plotting code
      ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
    }, name = "p")
    render_step_plot(output, input, "bad", function() stop("no data in this cohort"), name = "bad")
  }
  shiny::testServer(server, {
    for (i in 1:5) {
      img <- output$p
      expect_true(is.list(img) && nzchar(img$src))
    }
    expect_no_error(output$bad)            # the error is drawn as a message, not raised
  })
  expect_null(grDevices::dev.list())
})
