library(tidyverse)

all_designs <- read_csv("pdc_sim_results.csv") %>% 
  mutate_at(., .vars = c("N", "Scenario", "Design"), 
            ~factor(., ordered = TRUE))

analysis_sum <- all_designs %>% 
  group_by(N, Design, Scenario) %>% 
  summarise(Mv = mean(MeanVal), 
         MeanPerc = mean(PercOracle), 
         Above80 = mean(PercOracle > .8), 
         Above90 = mean(PercOracle > .9)) %>% 
  arrange(N, Scenario, Design)


cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#ggplot(data = sim_summaries, aes(x = PercOracle, color = N)) + 
#  geom_density() + 
#  scale_colour_manual(values = cbp1) + 
#  labs(x = "Percentage of Oracle Value")

cutoffHelper <- function(indf, cutoff.val, group.vars){
  to_return <- indf %>% group_by(N, Scenario, Design) %>% 
    summarise(Power = mean(PercOracle > cutoff.val),
              .groups = "drop_last") %>% 
    mutate(Cutoff = cutoff.val)
  
  return(to_return)
}

perc_oracle_power <- map_dfr(c(seq(0, .5, by = .1), seq(0.5, 1, by = .005)), ~cutoffHelper(indf = all_designs, cutoff.val = .))

plot <- ggplot(data = perc_oracle_power, aes(x = Cutoff, y = Power, color = Design)) + 
  geom_hline(aes(yintercept = 0.9), color = "black") +
  geom_vline(aes(xintercept = 0.9), color = "black") +
  geom_line(size = 0.8) + 
  # xlim(0, 1)+
  # ylim(0,1)+
  scale_x_continuous(breaks = c(0, .2, 0.4, .6, .8,  1)) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_colour_manual(values = cbp1) + 
  labs(x = "Percentage of Oracle Value Cutoff",
       y = "Probability of Attaining Value Above Cutoff") + 
  facet_grid(cols = vars(Scenario),
             rows = vars(N),
             labeller = label_both) + 
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"), axis.text = element_text(size = 6))

# ggsave(filename = "powerPlot.png", plot = plot, device = "png", width = 7, height = 7)

