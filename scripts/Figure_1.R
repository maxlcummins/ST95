library(ggtree)
library(treemapify)
library(grid)
library(gridExtra)
library(cowplot)

#Check to see if our dataframe exists prior to running 
if(!exists("Metadata")){
        source("scripts/ST95_all_process.R")
}

#Remove all objects from our environment except the one we need
rm(list=setdiff(ls(), "Metadata"))

#Define colours for sources
source_cols <- c(
        "brown"	,	#	Bovine	#	Source		
        "gold2"	,	#	Canine	#	Source		
        "springgreen4"	,#	Environment	#	Source		
        "cyan2"	,	#	Human	#	Source		
        "black"	,	#	Other	#	Source		
        "mediumorchid1", #	Poultry	#	Source
        "grey"
)

#Define sources
source_vars <- c(
        "Bovine",
        "Canine",
        "Environment",
        "Human",
        "Other",
        "Poultry",
        "Other")

#Bind source name and colour
names(source_cols) <- source_vars

#As above but for continents
continent_cols <- c("#DDCC77", #sand
                    "#EE3377",  #magenta
                    "#009988",   #teal
                    "#4dac26", #dark green,
                    "grey"
)
continent_names <- c("North America","Europe", "Oceania", "Asia", "Unknown")

names(continent_cols) <- continent_names


#Bind continent and source colours/names
var_col <- c(continent_cols, source_cols)

#Reassign the metadata dataframe to something shorter
df <- Metadata

#Replace NAs with Unknown in columns about source, year, continent
df[is.na(df$Revised_Source_Niche),"Revised_Source_Niche"] <- "Unknown"
df$Collection_Year <- as.character(df$Collection_Year)
df[is.na(df$Collection_Year),"Collection_Year"] <- "Unknown"
df[is.na(df$Continent),"Continent"] <- "Unknown"

#Manually define order of bars for bar graph (Based on descending count)
fig1_orders <- c("Human", "Poultry", "Environment", "Bovine", "Canine", "Other", "Unknown")
fig1_text <- sort(table(Metadata$Revised_Source_Niche), decreasing = TRUE) %>% as.data.frame()


#Define dynamic maximum limit of y axis based on highest count of variable type between Continent and Revised_Source_Niche in Metadata
maxlimit <- max(c(table(df$Revised_Source_Niche), table(Metadata$Continent)))

#Plot counts of a given source
fig1a <- ggplot(df,
                aes(x = Revised_Source_Niche, fill = Revised_Source_Niche)) +
        geom_bar(color = "black") +
        theme_minimal() +
        scale_y_continuous(expand = c(0, 0), limits = c(0, maxlimit+25)) +
        scale_x_discrete(limits = fig1_orders) +
        scale_fill_manual(
                name = "A - Isolation Source     ",
                aesthetics = c("colour", "fill"),
                values = var_col,
                na.value = "white"
        ) +
        xlab(label = "Isolation Source") +
        ylab(label = "Genomes (count)")

###############################################################################################

#Extract the data from figure 1a plot to help with plotting the number of strains above a given bar
fig1a_data <- ggplot_build(fig1a)
fig1a_data <- as.data.frame(fig1a_data$data)
fig1a <- fig1a +  annotate("text",
                 x = fig1a_data$x,
                 y = fig1a_data$y + 15,
                 label = fig1a_data$count,
                 size = 3)

#####################################################################################

#Manually define order of bars for bar graph (Based on descending count)
fig2_orders <- c("North America", "Europe", "Oceania", "Asia", "Unknown")
fig2_text <- sort(table(df$Continent), decreasing = TRUE) %>% as.data.frame()

#Plot count of strains form each continent
fig1b <- ggplot(df,
       aes(x = Continent, fill = Continent)) +
        geom_bar(color = "black") +
        theme_minimal() +
        scale_y_continuous(expand = c(0,0), limits = c(0, maxlimit+25)) +
        scale_x_discrete(limits = fig2_orders) +
        scale_fill_manual(name = "B - Continent of Origin",
                          aesthetics = c("colour", "fill"),
                          values = var_col,
                          na.value = "white") +
        xlab(label = "Continent of Origin") +
        ylab(label = "Genomes (count)")


###############################################################################################

#Extract the data from figure 1b plot to help with plotting the number of strains above a given bar
fig1b_data <- ggplot_build(fig1b)
fig1b_data <- as.data.frame(fig1b_data$data)
fig1b <- fig1b +  annotate("text",
                  x = fig1b_data$x,
                  y = fig1b_data$y + 15,
                  label = fig1b_data$count,
                  size = 3)

#####################################################################################

#Plot count of strains form each year
fig1c <- ggplot(df,
       aes(x = Collection_Year, fill = Collection_Year)) +
        geom_bar(color = "black") +
        theme_minimal() +
        scale_y_continuous(expand = c(0,0)) +
        scale_fill_manual(name = "Year",
                          values = c(rep("white", 40), "grey"),
                          na.value = "white") +
        theme(legend.position = 'none', axis.text.x = element_text(hjust=1, angle = 90)) +
        xlab(label = "Collection Year") +
        ylab(label = "Genomes (count)")
        

#Extract the legends from figures 1a and 1b
legend1a <- get_legend(fig1a)
legend1b <- get_legend(fig1b)

#Create a blank plot for using to space out subplots later
spacer <- ggplot() + theme_minimal()

#Plot the two legends together
#Nrow should be two but was set to three for easy padding
legends <- plot_grid(spacer, legend1a, spacer, legend1b, spacer, ncol=1, nrow=5, rel_heights = c(0.1, 2/5, 0, 2/5, 2/5))

#Re-assign figures 1a and 1b without their legends
fig1a <- fig1a + theme(legend.position = 'none')
fig1b <- fig1b + theme(legend.position = 'none')

#Plot figure 1a and b together
plots <- grid.arrange(fig1a, fig1b, ncol=2)

#Plot figure 1a and 1b
fig1ab <- plot_grid(spacer, plots, ncol=1, rel_heights = c(0.8/5, 4.2/5)) 

#Plot figure 1ab and 1c together
fig1_noleg <- plot_grid(fig1ab, fig1c, ncol=1)

#Plot fig1 (fig1abc) with its legend
fig1 <- plot_grid(fig1_noleg, legends, ncol=2, rel_widths = c(4/5, 1/5))

#Add in A, B and C labels for figures 1a, 1b and 1c.
fig1 <- fig1 + 
        annotate(geom = "text", x = 0.018, y = 0.95, label = "A",
                 parse = TRUE,
                 size = 12) + 
        annotate(geom = "text", x = 0.425, y = 0.95, label = "B",
                 parse = TRUE,
                 size = 12)+ 
        annotate(geom = "text", x = 0.018, y = 0.5, label = "C",
                 parse = TRUE,
                 size = 12)

#Visualise Figure 1
fig1

