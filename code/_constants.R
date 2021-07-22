## system info

nodename <- Sys.info()["nodename"]
n_cores <- parallel::detectCores()

## image, design, analysis info

n_vert <- 20484  ## surface hcp mesh
n_trs <- c(
  Axcpt   = 1220,
  # Axcpt_proactive  = 1220,
  # Axcpt_reactive   = 1220,
  Cuedts  = 1300,
  # Cuedts_proactive = 1300,
  # Cuedts_reactive  = 1300,
  Stern   = 1200,
  # Stern_proactive  = 1200,
  # Stern_reactive   = 1200,
  Stroop  = 1080
  # Stroop_proactive = 1080,
  # Stroop_reactive  = 1180
)

tasks <- c("Axcpt", "Cuedts", "Stern", "Stroop")

subjs_ub55 <- read.csv(here::here("..", "ub55", "in", "ub55_subjects.txt"))
subjs_test <- c()
subjs_retest <- c()



dmcc34 <- c(
  22, 77, 78, 86, 87, 91, 93, 99, 101, 103, 105, 107, 110, 127, 130, 139, 140,
  144, 148, 172, 175, 185, 189, 219, 301, 303, 306, 314, 340, 346, 347, 349, 350, 353
)

target_trs <- list(
  Axcpt = 7:9,
  Cuedts = 9:10,
  Stern = 11:12,
  Stroop = 2:4
)

cue_trs <- list(
  Axcpt = 8:9,  ### ????
  Cuedts = 9:10, ### ????
  Stern = 7:8,
  Stroop = 2:4
)