behav_wave12 <- list(
  Axcpt  = data.table::fread(here::here("in", "behav", "behavior-and-events_wave12_Axcpt.csv")),
  Cuedts = data.table::fread(here::here("in", "behav", "behavior-and-events_wave12_Cuedts.csv")),
  Stern  = data.table::fread(here::here("in", "behav", "behavior-and-events_wave12_Stern.csv")),
  Stroop = data.table::fread(here::here("in", "behav", "behavior-and-events_wave12_Stroop.csv"))
)

