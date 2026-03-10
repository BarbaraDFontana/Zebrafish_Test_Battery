# -----------------------------------------------------------
# Zebrafish Behavioral Battery Literature Analysis
# Fontana et al.
#
# This script analyzes zebrafish behavioral assays reported
# in studies using the Novel Tank Test (NTT).
#
# Input:
#   metanalise_zebrafish_screening.xlsx
#
# Output:
#   figures/
#     FigureA_MostCommonTests.pdf
#     FigureB_NTTParameters.pdf
#     FigureC_NumberTestsPerStudy.pdf
#     FigureD_TestCoOccurrenceHeatmap.pdf
# -----------------------------------------------------------


################################################
# Libraries
################################################

library(tidyverse)
library(readxl)
library(janitor)
library(stringr)
library(tidytext)
library(scico)


################################################
# Load and clean data
################################################

data_file <- "zebrafish_screening.xlsx"

df <- read_excel(data_file) %>%
  clean_names() %>%
  mutate(across(where(is.character), ~str_trim(.))) %>%
  mutate(across(where(is.character), ~na_if(., ""))) %>%
  mutate(across(where(is.character), ~na_if(., "na"))) %>%
  mutate(across(c(ntt_present),
                ~case_when(
                  str_to_lower(.) %in% c("yes","y") ~ "yes",
                  str_to_lower(.) %in% c("no","n") ~ "no",
                  TRUE ~ NA_character_
                )))

# protect against duplicated studies
df <- df %>%
  distinct(doi, .keep_all = TRUE)


################################################
# A — Most common behavioral tests
################################################

tests <- df %>%
  separate_rows(behavioral_test, sep = ";|,") %>%
  mutate(test = str_to_lower(str_trim(behavioral_test))) %>%
  filter(test != "", !is.na(test)) %>%
  mutate(test = case_when(
    
    str_detect(test,"novel tank|ntt") ~ "novel tank test",
    
    str_detect(test,"light.*dark|ldt|scototaxis") ~
      "light–dark test",
    
    str_detect(test,"open field|\\boft\\b") ~
      "open field test",
    
    str_detect(test,"social preference|social") ~
      "social preference test",
    
    str_detect(test,"shoal|shoaling") ~
      "shoaling test",
    
    str_detect(test,"aggression|mirror|mirror-induced|miat") ~
      "aggression test",
    
    str_detect(test,
               "novel object|\\bnor\\b|novel approach|\\bnat\\b|object approach") ~
      "novel object approach/recognition",
    
    str_detect(test,"y[- ]?maze") ~
      "y-maze test",
    
    str_detect(test,"t[- ]?maze") ~
      "t-maze test",
    
    str_detect(test,"place preference|conditioned|cpp") ~
      "conditioned place preference",
    
    str_detect(test,"predator") ~
      "predator exposure test",
    
    str_detect(test,"color preference") ~
      "color preference test",
    
    TRUE ~ test
  ))


test_counts <- tests %>%
  filter(test != "novel tank test") %>%
  count(test, sort = TRUE) %>%
  slice_max(n, n = 10)

figA <- ggplot(test_counts,
               aes(x = reorder(test,n), y = n, fill = test)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_scico_d(palette = "navia") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  labs(
    title = "Most frequent tests in NTT studies",
    x = "Behavioral assay",
    y = "Number of studies"
  )


################################################
# B — NTT behavioral parameters
################################################

ntt_params <- df %>%
  filter(str_to_lower(ntt_present) == "yes") %>%
  separate_rows(behavioral_parameters_measured, sep=";|,|\\n|\\r") %>%
  mutate(parameter_raw =
           str_to_lower(str_trim(behavioral_parameters_measured))) %>%
  filter(parameter_raw != "", !is.na(parameter_raw)) %>%
  mutate(parameter = case_when(
    
    # vertical zones
    str_detect(parameter_raw,
               "time.*(top|upper)|time.*(bottom|lower)|(top|upper).*time|(bottom|lower).*time|duration.*(top|upper|bottom|lower)") ~
      "time in top/bottom",
    
    # latency
    str_detect(parameter_raw,
               "latency.*(top|upper|zone|bottom|lower)") ~
      "latency to enter zone",
    
    # locomotion distance
    str_detect(parameter_raw,
               "distance travelled|distance traveled|distance moved|total distance|path length|locomotion") ~
      "distance travelled",
    
    # velocity
    str_detect(parameter_raw,
               "velocity|speed|mean speed|average speed|swimming speed") ~
      "velocity",
    
    # angular movement
    str_detect(parameter_raw,
               "angular velocity|meander") ~
      "angular velocity",
    
    # turning
    str_detect(parameter_raw,
               "turn angle|absolute turn angle") ~
      "turn angle",
    
    # zone transitions
    str_detect(parameter_raw,
               "entries|frequency.*zone|number of entries|transitions|crossings") ~
      "zone entries",
    
    # freezing
    str_detect(parameter_raw,
               "freezing|immobility|inactivity") ~
      "freezing/immobility",
    
    # erratic movement
    str_detect(parameter_raw,
               "erratic|dart") ~
      "erratic movement",
    
    # risk behaviour
    str_detect(parameter_raw,
               "risk assessment|hesitation") ~
      "risk assessment",
    
    # exploration
    str_detect(parameter_raw,
               "trajectory|exploration") ~
      "exploratory activity",
    
    # bottom distance
    str_detect(parameter_raw,
               "distance.*(bottom|lower)|(bottom|lower).*distance") ~
      "distance from bottom",
    
    str_detect(parameter_raw,"\\banxiety\\b") ~
      NA_character_,
    
    TRUE ~ "other"
  )) %>%
  drop_na(parameter)


param_counts <- ntt_params %>%
  filter(parameter != "other") %>%
  distinct(doi, parameter) %>%
  count(parameter, sort = TRUE) %>%
  slice_max(n, n = 10)

figB <- ggplot(param_counts,
               aes(x = reorder(parameter,n), y = n, fill = parameter)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), hjust = -0.1, size = 4) +
  coord_flip() +
  scale_fill_scico_d(palette = "navia") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  labs(
    title = "Popular NTT parameters",
    x = "Behavioral parameter",
    y = "Number of studies"
  ) +
  ylim(0, max(param_counts$n) * 1.1)


################################################
# C — Number of tests per study
################################################

df_testcount <- tests %>%
  group_by(doi) %>%
  summarise(n_tests = n_distinct(test)) %>%
  ungroup()

test_summary <- df_testcount %>%
  mutate(n_tests_group =
           ifelse(n_tests >= 4, "4+", as.character(n_tests))) %>%
  count(n_tests_group) %>%
  mutate(
    pct = n / sum(n),
    label = paste0(n_tests_group,"\n",
                   scales::percent(pct, accuracy = 1))
  )

test_summary$n_tests_group <- factor(
  test_summary$n_tests_group,
  levels = c("1","2","3","4+")
)

figC <- ggplot(test_summary,
               aes(x = 2, y = pct, fill = n_tests_group)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5,2.5) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 4) +
  scale_fill_scico_d(palette = "navia") +
  theme_void(base_size = 14) +
  labs(title = "Number of behavioral tests per study")


################################################
# D — Test co-occurrence heatmap
################################################

test_freq <- tests %>%
  count(test, sort = TRUE)

tests_filtered <- tests %>%
  filter(test %in% test_freq$test[test_freq$n >= 5])

tests_matrix <- tests_filtered %>%
  distinct(doi, test) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = test,
              values_from = value,
              values_fill = 0)

mat <- as.matrix(tests_matrix[,-1])

co_occurrence <- t(mat) %*% mat
diag(co_occurrence) <- NA

co_occurrence_df <- as.data.frame(as.table(co_occurrence)) %>%
  drop_na()

figD <- ggplot(co_occurrence_df,
               aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_scico(palette = "navia", direction = -1) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank()) +
  labs(
    title = "Behavioral test co-occurrence",
    fill = "N of studies"
  )


################################################
# Save figures
################################################

fig_dir <- "figures"
if(!dir.exists(fig_dir)) dir.create(fig_dir)

ggsave(file.path(fig_dir,"FigureA_MostCommonTests.pdf"),figA,8,6)
ggsave(file.path(fig_dir,"FigureB_NTTParameters.pdf"),figB,8,6)
ggsave(file.path(fig_dir,"FigureC_NumberTestsPerStudy.pdf"),figC,8,6)
ggsave(file.path(fig_dir,"FigureD_TestCoOccurrenceHeatmap.pdf"),figD,8,6)