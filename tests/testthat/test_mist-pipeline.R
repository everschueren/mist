library(mist)

mist_config <- readMistConfig('../../inst/extdata/small/mist_small_test.yml')
runMistPipeline(mist_config)