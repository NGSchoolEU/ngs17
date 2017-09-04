# Best 250 series
# http://www.imdb.com/chart/toptv/
## (1) read the new data with archivist

library(archivist)
series2017 <- aread("mi2-warsaw/RLadies/arepo/45aa16dc4dbf0d87e3e40eb9dc9d18ae")

## (2) or read the old data with Pogromcy Danych

library(PogromcyDanych)
serialeIMDB

## (3) or scrap the data from IMDB database

library(rvest)
library(dplyr)

# read links and titles 
page <- read_html("http://www.imdb.com/chart/toptv/")
series <- html_nodes(page, ".titleColumn a")
titles <- html_text(series)
links <- html_attr(series, "href")
codes <- sapply(strsplit(links, split = "/"), `[`, 3)

# read details for series
allSeries <- lapply(seq_along(codes), function(i) {
  tab <- read_html(paste0("http://www.imdb.com/title/",codes[i],"/epdate?ref_=ttep_ql_4")) %>%
    html_node("table") %>%
    html_table()
  data.frame(Serie = titles[i], tab[,1:4], Season = gsub(tab[,1], pattern="\\..*", replacement=""))
})

# put all together
series2017 <- do.call(rbind, allSeries)
series2017$UserVotes <- as.numeric(gsub(series2017$UserVotes, pattern = "[^0-9]", replacement = ""))
series2017 %>% 
  group_by(Serie) %>%
  mutate(id = seq_along(Serie)) %>%
  ungroup() -> series2017

##
## Time for some datavis
##
selected <- c("Gra o tron", "Breaking Bad", 
              "Sherlock", "Westworld")

dat <- series2017 %>%
  filter(Serie %in% selected)

ggplot(dat, aes(id, UserRating)) + 
  geom_point(aes(color=Season, size=UserVotes)) +
  geom_smooth(se=FALSE, color="grey30") + 
  facet_grid(Serie~.) +
  theme_light() + theme(legend.position="none") +
  scale_color_brewer(palette = 1, type = "qual") +
  ggtitle("User ratings for selected TV series") + 
  xlab("Episode No")


#
# Let's practice
#

mySer <- filter(dat, Serie == "Breaking Bad")
head(mySer)

ggplot(mySer, aes(x = id, y = UserRating)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  geom_smooth(method = "lm", se=FALSE, color="red") 
  
ggplot(mySer, aes(x = id, 
                      y = UserRating,
                      size=UserVotes)) +
  geom_point(aes(color = Season)) +
  geom_smooth() +
  geom_smooth(aes(color = Season)) +
  theme_bw()
   
ggplot(mySer, aes(x = Season, 
                      fill = UserRating > 8.5)) +
  geom_bar(position = "dodge")


ggplot(mySer, aes(x = Season, 
                      y = UserRating)) +
  geom_point(position = "jitter") 

