# -----------------------------------------------------------
# Zebrafish Behavioral Battery Survey Analysis
# Fontana et al.
#
# This script processes survey responses from zebrafish
# behavioral neuroscience laboratories and generates
# the figures used in the manuscript.
#
# Input:
#   Forms_Data.csv
#
# Output:
#   figures/*.pdf
# -----------------------------------------------------------

################################################
# Libraries
################################################

library(tidyverse)
library(janitor)
library(stringr)
library(tidytext)
library(forcats)
library(scales)
library(scico)

################################################
# Paths
################################################

project_dir <- "Path with data"
data_file   <- file.path(project_dir, "Forms_Data.csv")

raw <- readr::read_csv(data_file, show_col_types = FALSE)

################################################
# Clean column names
################################################

df <- raw %>%
  clean_names()

################################################
# Helper functions
################################################

# Convert free-text times to minutes
to_minutes <- function(x) {
  x <- tolower(x)
  
  ifelse(is.na(x) | x == "", NA_real_,
         purrr::map_dbl(x, ~{
           s <- .x
           
           s <- str_replace_all(
             s,
             "(\\d+(?:\\.\\d+)?)\\s*h(?:our)?s?",
             function(m) as.numeric(str_extract(m, "\\d+(?:\\.\\d+)?")) * 60
           )
           
           mins <- as.numeric(str_extract_all(s, "\\d+(?:\\.\\d+)?")[[1]])
           
           if (length(mins) == 0) return(NA_real_)
           mean(mins)
         })
  )
}

# Normalize separators
norm_sep <- function(x){
  x %>%
    str_replace_all("→|->|—|>|–|,|;|/|\\+", " -> ") %>%
    str_replace_all("\\s+", " ") %>%
    str_replace_all("\\s*->\\s*", " -> ") %>%
    str_trim()
}

################################################
# Extract survey fields
################################################

df2 <- df %>%
  mutate(
    
    multi_same_day = str_to_lower(coalesce(
      df[[grep("^do_you_perform_more_than_one_behavioral_test", names(df))]],
      NA_character_
    )),
    
    randomize_order = str_to_lower(coalesce(
      df[[grep("^do_you_randomize_the_order_of_tests", names(df))]],
      NA_character_
    )),
    
    n_tests_per_subject = as.numeric(readr::parse_number(
      df[[grep("^how_many_tests_do_you_typically_run", names(df))]]
    )),
    
    order_raw = df[[grep("^if_yes_in_which_order", names(df))]],
    
    intertest_minutes = to_minutes(
      df[[grep("^what_is_the_average_inter_test_interval", names(df))]]
    ),
    
    acclimation_minutes = to_minutes(
      df[[grep("^what_is_the_average_acclimation_time", names(df))]]
    )
  )

################################################
# A — Top tests mentioned in survey
################################################

rules <- tibble::tribble(
  ~group,        ~pattern,
  "NTT",         "(?i)novel\\s*ty?\\s*tank|\\bntt\\b",
  "LDT",         "(?i)light\\s*-?\\s*dark|\\bldt\\b",
  "Shoaling Test","(?i)shoal",
  "Social Preference","(?i)social",
  "Aggression","(?i)aggression|mirror",
  "OFT","(?i)open\\s*field|\\boft\\b"
)

map_to_group <- function(token){
  if (is.na(token) || token == "") return(NA_character_)
  m <- rules$group[str_detect(token, rules$pattern)]
  if (length(m)==0) NA_character_ else m[[1]]
}

orders_long <- df2 %>%
  transmute(order_raw) %>%
  filter(!is.na(order_raw)) %>%
  mutate(seq_norm = norm_sep(order_raw)) %>%
  separate_rows(seq_norm, sep="->") %>%
  mutate(
    token = str_to_lower(str_squish(seq_norm)),
    group = purrr::map_chr(token, map_to_group)
  )

test_counts <- orders_long %>%
  filter(!is.na(group)) %>%
  count(group, sort=TRUE)

p_tests_grouped <- test_counts %>%
  ggplot(aes(x=fct_reorder(group,n), y=n, fill=group))+
  geom_col(width=0.7)+
  coord_flip()+
  scale_fill_scico_d(palette="navia")+
  labs(
    title="Top behavioral tests mentioned in the survey",
    x=NULL,
    y="Mentions"
  )+
  theme_classic(base_size=12)+
  theme(legend.position="none")

################################################
# B — Number of tests per subject
################################################

collapse4plus <- TRUE

tests_tbl <- df2 %>%
  transmute(k=n_tests_per_subject) %>%
  filter(!is.na(k), k>0) %>%
  mutate(k_lab=if_else(k>=4,"4+",as.character(k))) %>%
  count(k_lab,name="n") %>%
  mutate(
    pct=n/sum(n),
    lbl=paste0(k_lab,"\n",percent(pct,accuracy=1))
  )

p_tests_donut <- ggplot(tests_tbl,aes(x=2,y=pct,fill=k_lab))+
  geom_col(width=1)+
  coord_polar(theta="y")+
  xlim(0.5,2.5)+
  geom_text(aes(label=lbl),
            position=position_stack(vjust=0.5),
            size=4)+
  scale_fill_scico_d(palette="navia",name="Tests/fish")+
  labs(title="How many tests per subject (fish)?")+
  theme_void(base_size=12)

################################################
# C — Reasoning keywords
################################################

reasoning <- df2 %>%
  transmute(text=df[[grep("^what_is_your_reasoning",names(df))]]) %>%
  filter(!is.na(text)) %>%
  mutate(text=str_to_lower(text))

data("stop_words")

extra_stop <- tibble(
  word = c(
    "test","tests","tank","fish",
    "behavior","behaviour","social",
    "field","sense","reason","makes",
    "initial","due","component",
    "combined","assessment","aggression"
  )
)

word_freq <- reasoning %>%
  unnest_tokens(word,text) %>%
  anti_join(stop_words,by="word") %>%
  anti_join(extra_stop,by="word") %>%
  count(word,sort=TRUE) %>%
  filter(n>1)

top_reasoning <- ggplot(slice_head(word_freq,n=25),
                        aes(x=fct_reorder(word,n),y=n,fill=word))+
  geom_col(width=0.7)+
  coord_flip()+
  scale_fill_scico_d(palette="navia")+
  labs(
    title="Top words in reasoning behind order",
    x=NULL,
    y="Count"
  )+
  theme_classic(base_size=12)+
  theme(legend.position="none")

################################################
# D — Multiple tests same day
################################################

multi_summary <- df2 %>%
  transmute(answer=case_when(
    n_tests_per_subject>=2~"Yes",
    n_tests_per_subject==1~"No"
  ))%>%
  filter(!is.na(answer))%>%
  count(answer)%>%
  mutate(
    pct=n/sum(n),
    lbl=paste0(answer,"\n",percent(pct))
  )

p_donut <- ggplot(multi_summary,aes(x=2,y=pct,fill=answer))+
  geom_col(width=1)+
  coord_polar(theta="y")+
  xlim(0.5,2.5)+
  geom_text(aes(label=lbl),
            position=position_stack(vjust=0.5),
            size=4)+
  scale_fill_manual(values=c("Yes"="#4A8C7F","No"="#001A33"))+
  labs(title="Multiple behavioral tests on the same day")+
  theme_void()

################################################
# Save figures
################################################

fig_dir <- "figures"
dir.create(fig_dir,showWarnings=FALSE)

ggsave(
  filename = file.path(fig_dir, "tests_grouped.pdf"),
  plot = p_tests_grouped,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(fig_dir, "tests_distribution_donut.pdf"),
  plot = p_tests_donut,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(fig_dir, "tests_reasoning.pdf"),
  plot = top_reasoning,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(fig_dir, "multi_tests_donut.pdf"),
  plot = p_donut,
  width = 8,
  height = 6
)