# TODO: Refactor this hard-coded code into functions

if(FALSE){
  #' Temporary workaround - change the percentage of optimal to only be stage two versus baseline
  
  RenameOracleSummary <- function(oracle.summary){
    oracle_renamed <- oracle.summary %>% 
      rename(DifOracleYSum = DifOracle,
             PercOraclYSum = PercOracle) %>% 
      mutate(DifOracle = DifOracleStg2,
             PercOracle = PercOracleStg2)
    
    return(oracle_renamed)
  }
  
  oracle_summary <- RenameOracleSummary(oracle_summary)
  
  all_designs <- oracle_summary %>% 
    mutate_at(., .vars = c("N", "Scenario"), 
              ~factor(., ordered = TRUE))
  cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  #ggplot(data = sim_summaries, aes(x = PercOracle, color = N)) + 
  #  geom_density() + 
  #  scale_colour_manual(values = cbp1) + 
  #  labs(x = "Percentage of Oracle Value")
  
  cutoffHelper <- function(indf, cutoff.val, group.vars){
    to_return <- indf %>% group_by(N, Scenario) %>% 
      summarise(Power = mean(PercOracle >= cutoff.val),
                .groups = "drop_last") %>% 
      mutate(Cutoff = cutoff.val)
    
    return(to_return)
  }
  
  cutoffHelperDif <- function(indf, cutoff.val, group.vars){
    to_return <- indf %>% group_by(N, Scenario) %>% 
      summarise(Power = mean(DifOracle <= cutoff.val),
                .groups = "drop_last") %>% 
      mutate(Cutoff = cutoff.val)
    
    return(to_return)
  }
  
  perc_oracle_power <- map_dfr(c(seq(0, .5, by = .1), seq(0.5, 1, by = .0025)), ~cutoffHelper(indf = all_designs, cutoff.val = .))
  
  makePlotForScenario <- function(df){
    out <- ggplot(data = df, aes(x = Cutoff, y = Power, color = N)) + 
      # geom_hline(aes(yintercept = 0.9), color = "black") +
      # geom_vline(aes(xintercept = 0.9), color = "black") +
      geom_line(size = 0.8) + 
      # xlim(0, 1)+
      # ylim(0,1)+
      scale_x_continuous(breaks = seq(.7, 1, by = .1),
                         minor_breaks	= seq(.75, .95, by = .1),
                         limits = c(0.7, 1)) +
      scale_y_continuous(breaks = seq(0.5, 1, by = .1),
                         minor_breaks	= NULL, #seq(.05, .95, by = .1),
                         limits = c(0.5, 1)) +
      scale_colour_manual(values = cbp1) + 
      labs(x = "Percentage of Optimal Value Cutoff",
           y = "Probability of Attaining Value Above Cutoff") + 
      theme_minimal() #+
    #theme(panel.grid.minor = element_blank())
    
    return(out)
  }
  
  makePlotForScenarioV2 <- function(df){
    out <- ggplot(data = df, aes(x = Cutoff, y = Power, color = N)) + 
      # geom_hline(aes(yintercept = 0.9), color = "black") +
      # geom_vline(aes(xintercept = 0.9), color = "black") +
      geom_line(size = 1.2) + 
      # xlim(0, 1)+
      # ylim(0,1)+
      scale_x_continuous(breaks = seq(0.6, 1, by = .1),
                         minor_breaks	= seq(.65, .95, by = .1),
                         limits = c(.6, 1)) +
      scale_y_continuous(breaks = seq(0, 1, by = .1),
                         minor_breaks	= NULL, #seq(.05, .95, by = .1),
                         limits = c(0, 1)) +
      scale_colour_manual(values = cbp1) + 
      labs(x = "Percentage of Optimal Value Cutoff",
           y = "Probability of Attaining Value Above Cutoff") +  
      theme_minimal() #+
    #theme(panel.grid.minor = element_blank())
    
    return(out)
  }
  
  makePlotForScenarioV3 <- function(df, n.replicates){
    out <- ggplot(data = df, aes(x = Cutoff, y = Power, color = N)) + 
      # geom_hline(aes(yintercept = 0.9), color = "black") +
      # geom_vline(aes(xintercept = 0.9), color = "black") +
      geom_line(size = 1.2) + 
      # xlim(0, 1)+
      # ylim(0,1)+
      scale_x_reverse(breaks = seq(.8, 1, by = .1),
                      minor_breaks	= seq(.85, .95, by = .1),
                      limits = c(1, .8)) +
      scale_y_continuous(breaks = seq(0, 1, by = .1),
                         minor_breaks	= NULL, #seq(.05, .95, by = .1),
                         limits = c(0, 1)) +
      scale_colour_manual(values = cbp1) + 
      labs(title = "Probability of Attaining Percentage of Optimal Value",
           subtitle = paste0("Simulation replicates per sample size: ", n.replicates),
           x = "Percentage of Optimal Value Cutoff",
           y = "Probability of Attaining Value Above Cutoff") +  
      theme_minimal() #+
    #theme(panel.grid.minor = element_blank())
    
    return(out)
  }
  
  dif_oracle_power <- map_dfr(c(seq(0, .25, by = .0025)), ~cutoffHelperDif(indf = all_designs, cutoff.val = .))%>% 
    mutate(Cutoff = -Cutoff)
  
  makePlotForScenarioDif <- function(df){
    out <- ggplot(data = df, aes(x = Cutoff, y = Power, color = N)) + 
      # geom_hline(aes(yintercept = 0.9), color = "black") +
      # geom_vline(aes(xintercept = 0.9), color = "black") +
      geom_line(size = 1.2) + 
      # xlim(0, 1)+
      # ylim(0,1)+
      scale_x_continuous(breaks = seq(-0.25, 0, by = .05),
                         minor_breaks	= seq(-.225, -.025, by = .05),
                         limits = c(-.25, 0)) +
      scale_y_continuous(breaks = seq(0, 1, by = .1),
                         minor_breaks	= NULL, #seq(.05, .95, by = .1),
                         limits = c(0, 1)) +
      scale_colour_manual(values = cbp1) + 
      labs(x = "Difference from Optimal Value",
           y = "Probability of Attaining Value Above Cutoff") + 
      theme_minimal() #+
    #theme(panel.grid.minor = element_blank())
    
    return(out)
  }
  
  makePlotForScenarioDifV2 <- function(df){
    out <- ggplot(data = df, aes(x = Cutoff, y = Power, color = N)) + 
      # geom_hline(aes(yintercept = 0.9), color = "black") +
      # geom_vline(aes(xintercept = 0.9), color = "black") +
      geom_line(size = 1.2) + 
      # xlim(0, 1)+
      # ylim(0,1)+
      scale_x_continuous(breaks = seq(-0.15, 0, by = .05),
                         minor_breaks	= seq(-.125, -.025, by = .05),
                         limits = c(-.15, 0)) +
      scale_y_continuous(breaks = seq(0, 1, by = .1),
                         minor_breaks	= NULL, #seq(.05, .95, by = .1),
                         limits = c(0, 1)) +
      scale_colour_manual(values = cbp1) + 
      labs(x = "Difference from Optimal Value",
           y = "Probability of Attaining Value Above Cutoff") + 
      theme_minimal() #+
    #theme(panel.grid.minor = element_blank())
    
    return(out)
  }
  
  plotDif <- makePlotForScenarioDifV2(dif_oracle_power)
  
  plotSc1 <- makePlotForScenarioV2(filter(perc_oracle_power, Cutoff >= 0.6))
  plotSc1 <- makePlotForScenarioV3(filter(perc_oracle_power, Cutoff >= 0.8), 3500)
  
  plotDif
  plotSc1
  
  #ggsave(plot = plotSc1, filename = "PercOraclePlotFourTrtSgSc1AltRespDist__630.png", height = 4, width = 5.5, units = "in", device = "png")
  #ggsave(plot = plotDif, filename = "DifOraclePlotFourTrtSgSc2AltRespDist_630.png", height = 4, width = 5.5, units = "in", device = "png")
  
}
