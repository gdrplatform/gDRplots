context("Test constants")

test_that("calling const works fine", {
  json_path <- system.file(package = "gDRplots", "settings.json")
  s <- gDRutils::get_settings_from_json(json_path = json_path)
  expect_identical(get_metrics_to_transform(), s$METRICS_TRANSFORMED)
  
})
