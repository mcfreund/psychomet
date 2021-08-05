behav <- list(
  Axcpt = data.table::fread(here::here("..", "ub55", "in", "ub55_axcpt__behav-events.csv")) %>% rename(resp = target.resp),
  Cuedts = data.table::fread(here::here("..", "ub55", "in", "ub55_cuedts_behav-events.csv")),
  Stern = data.table::fread(here::here("..", "ub55", "in", "ub55_stern_behav-events.csv")),
  Stroop = data.table::fread(here::here("..", "ub55", "in", "ub55_stroop_behav-events.csv"))
)

behav$Stern$resp <- ifelse(behav$Stern$resp == 1, 2, ifelse(behav$Stern$resp == 2, 1, NA))
behav$Stern$cresp <- ifelse(behav$Stern$cresp == 1, 2, ifelse(behav$Stern$cresp == 2, 1, NA))
