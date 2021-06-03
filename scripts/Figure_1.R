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


source_cols <- c(
        "brown"	,	#	Bovine	#	Source		
        "gold2"	,	#	Canine	#	Source		
        "springgreen4"	,#	Environment	#	Source		
        "cyan2"	,	#	Human	#	Source		
        "black"	,	#	Other	#	Source		
        "mediumorchid1", #	Poultry	#	Source
        "grey"
)

source_vars <- c(
        "Bovine",
        "Canine",
        "Environment",
        "Human",
        "Other",
        "Poultry",
        "Other")

names(source_cols) <- source_vars

continent_cols <- c("#DDCC77", #sand
                    "#EE3377",  #magenta
                    "#009988",   #teal
                    "#4dac26", #dark green,
                    "grey"
)

continent_names <- c("North America","Europe", "Oceania", "Asia", "Unknown")

names(continent_cols) <- continent_names


var_col <- c(continent_cols, source_cols)


df <- Metadata

df[is.na(df$Revised_Source_Niche),"Revised_Source_Niche"] <- "Unknown"
df$Collection_Year <- as.character(df$Collection_Year)
df[is.na(df$Collection_Year),"Collection_Year"] <- "Unknown"
df[is.na(df$Continent),"Continent"] <- "Unknown"


fig1_orders <- c("Human", "Poultry", "Environment", "Bovine", "Canine", "Other", "Unknown")

fig1_text <- sort(table(Metadata$Revised_Source_Niche), decreasing = TRUE) %>% as.data.frame()


#Define dynamic maximum limit of y axis based on highest count of variable type between Continent and Revised_Source_Niche in Metadata
maxlimit <- max(c(table(df$Revised_Source_Niche), table(Metadata$Continent)))

fig1a <- ggplot(df,
                aes(x = Revised_Source_Niche, fill = Revised_Source_Niche)) +
        geom_bar(color = "black") +
        theme_minimal() +
        scale_y_continuous(expand = c(0, 0), limits = c(0, maxlimit+25)) +
        scale_x_discrete(limits = fig1_orders) +
        #geom_text(
        #        aes(x = fig1_text$Var1,
        #            y = fig1_text$Freq,
        #            label = fig1_text$Freq
        #            ),
        #        position = position_dodge(width = 1),
        #        size = 3
        #) +
        scale_fill_manual(
                name = "A - Isolation Source     ",
                aesthetics = c("colour", "fill"),
                values = var_col,
                na.value = "white"
        ) +
        xlab(label = "Isolation Source") +
        ylab(label = "Genomes (count)")

###############################################################################################

numb<- ggplot_build(fig1a)
hi<- as.data.frame(numb$data)
fig1a <- fig1a +  annotate("text",
                 x = hi$x,
                 y = hi$y + 15,
                 label = hi$count,
                 size = 3)

#####################################################################################


fig2_orders <- c("North America", "Europe", "Oceania", "Asia", "Unknown")

fig2_text <- sort(table(df$Continent), decreasing = TRUE) %>% as.data.frame()

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

numb<- ggplot_build(fig1b)
hi<- as.data.frame(numb$data)
fig1b <- fig1b +  annotate("text",
                  x = hi$x,
                  y = hi$y + 15,
                  label = hi$count,
                  size = 3)

#####################################################################################


fig1c <- ggplot(df,
       aes(x = Collection_Year, fill = Collection_Year)) +
        geom_bar(color = "black") +
        theme_minimal() +
        scale_y_continuous(expand = c(0,0)) +
        #geom_text(aes(label = paste0(percentage,"% (", counts, "/668)")), position = position_stack(vjust=0.5), size = 3) +
        #scale_x_discrete(limits = fig2_orders) +
        scale_fill_manual(name = "Year",
                          #aesthetics = c("colour", "fill"),
                          values = c(rep("white", 40), "grey"),
                          na.value = "white") +
        theme(legend.position = 'none', axis.text.x = element_text(hjust=1, angle = 90)) +
        xlab(label = "Collection Year") +
        ylab(label = "Genomes (count)")
        


legend1a <- get_legend(fig1a)
legend1b <- get_legend(fig1b)
spacer <- ggplot() + theme_minimal()

#nrow should be two but was set to three for easy padding
legends <- plot_grid(spacer, legend1a, spacer, legend1b, spacer, ncol=1, nrow=5, rel_heights = c(0.1, 2/5, 0, 2/5, 2/5))

fig1a <- fig1a + theme(legend.position = 'none')

fig1b <- fig1b + theme(legend.position = 'none')

plots <- grid.arrange(fig1a, fig1b, ncol=2)

fig1ab <- plot_grid(plots, legends, ncol=2, rel_widths = c(4/5, 1/5)) 

fig1ab <- plot_grid(spacer, plots, ncol=1, rel_heights = c(0.8/5, 4.2/5)) 


fig1_noleg <- plot_grid(fig1ab, fig1c, ncol=1)
                  #, rel_widths = c(4/5, 1/5)) 

fig1 <- plot_grid(fig1_noleg, legends, ncol=2, rel_widths = c(4/5, 1/5))

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

fig1

