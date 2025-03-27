
ELF_supers_coords_clean
coords <- get_coords(exp_metaData = PCC_exp_metaData, target_inds = sample(names(PCC_exp_metaData), 100),  
                     easting_col = "DD.easting", 
                     northing_col = "DD.northing", 
                     age_col = "age.class", 
                     year_col = "year")

coords_clean <- subset(coords, subset = !coords[,"UTM.easting"] == ".")
coords_clean[,2:3] <- sapply(coords_clean[,c("UTM.easting", "UTM.northing")], as.numeric)
coords_clean <-  subset(coords_clean, subset = !is.na(coords_clean[,"UTM.easting"]))
coords <- coords_clean
# Define CRS using the sf package
crs.geo <- st_crs("+proj=utm +zone=12 +datum=WGS84")
crs.geo <- st_crs("+proj=longlat +datum=WGS84")


# Make site into a spatial feature using sf package; assign its projection
site <- coords %>%    
  st_as_sf(coords = c("UTM.easting", "UTM.northing"),
           crs = crs.geo)

site.sp<- as(site, "Spatial")

NLCD <- get_nlcd(
  template = site.sp,
  year = 2019,
  label = "test", 
  force.redo = T
)

# site.lim <- ELF_supers_coords_clean %>%    
#   st_as_sf(coords = c("UTM.easting", "UTM.northing"),
#            crs = crs.geo)

tm_shape(NLCD) +
  tm_raster(alpha = 0.75, legend.show = FALSE) +
  tm_shape(site) + tm_bubbles(palette = colors,
                              col = "ID",
                              size = 0.25, 
                              legend.col.show = FALSE)
  


NLCD <- get_nlcd(
  template = site.sp,
  year = 2016,
  label = "test")
  
# Define crs using the raster pacakge
crs.geo.r <- CRS("+proj=utm +zone=12 +datum=WGS84") # set CRS

NED <- projectRaster(NLCD, crs = crs.geo.r)

# Create a dataframe out of the NED raster. This is a necessary step for plotting rasters with ggplot2. 
NED_df1 = as.data.frame(NLCD, xy = TRUE)
head(NED_df1)
NED_df1 <- NED_df1 %>%
  drop_na()
# Rename the elevation column
NED_df <- NED_df1 %>%
  rename(elev_NED = NLCD.Land.Cover.Class)
# Plot the NED with the points overlaid using ggplot. See how the elevation in the points aligns with the NED
NED_elev<-ggplot() +
  geom_raster(data = NED_df, 
              aes(x = x, y = y, fill = elev_NED)) # x and y are defined as the column names, "x" and "y"; "elev_NED" is the elevation column

# Now add on the vegetation plots and their associated elevation from the elevdat data, and change the color of the raster.
p <- ggplot() +
  geom_raster(data = NED_df, aes(x = x, y = y, fill = elev_NED)) + 
  geom_sf(data = site, aes(color = ID)) 

p+ geom_point(aes(x = long, y = lat, color = ID), 
                            data = dat.df, size = 2.5) + 
  geom_path(aes(x = long, y = lat, color = ID), 
            data = dat.df) +
  theme(legend.position = "none") 
  #labs(title=plot.title) #+ 
  #scale_color_manual(values=met.brewer(palette, n = 18))
#geom_sf()
#north(anndat.dat) +


