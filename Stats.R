library(jsonlite)
library(ggplot2)
library(ggformula)
library(scales)
library(lubridate)
library(dplyr)
library(nortest)
library(kSamples)
library(psych)
library("gplots")
library(plotly)
library(geosphere)
library(cowplot)
library(gridExtra)
library(grid)

receiving_station ="IU2KQB"   # Receiving station ID (uploader callsign) as registered on SondeHub
sonda_serial <- "W4140233" #V3330922"#V4730127 #"W1950650"# V4730127" # V4730885"     # Target sonde
data_save_path <-"c:\\tools\\"
saveallimages <-"Y"
uploader_lat <- 0    # Do not valorize
uploader_lon <- 0     # Do not valorize

uploader_alt <- 75
R <- 6371000  

#  lat/lon/alt to (X, Y, Z)
geo_to_cartesian <- function(lat, lon, alt) {
  lat_rad <- lat * pi / 180
  lon_rad <- lon * pi / 180
  x <- (R + alt) * cos(lat_rad) * cos(lon_rad)
  y <- (R + alt) * cos(lat_rad) * sin(lon_rad)
  z <- (R + alt) * sin(lat_rad)
  return(c(x, y, z))
}



# bearing 
calculate_bearing <- function(lat1, lon1, lat2, lon2) {
  lat1_rad <- lat1 * pi / 180
  lon1_rad <- lon1 * pi / 180
  lat2_rad <- lat2 * pi / 180
  lon2_rad <- lon2 * pi / 180
  
  dLon <- lon2_rad - lon1_rad
  y <- sin(dLon) * cos(lat2_rad)
  x <- cos(lat1_rad) * sin(lat2_rad) - sin(lat1_rad) * cos(lat2_rad) * cos(dLon)
  bearing_rad <- atan2(y, x)
  bearing_deg <- (bearing_rad * 180 / pi + 360) %% 360  # 
  return(bearing_deg)
}

#Load data from http://api.v2.sondehub.org/sonde/V4730886
#"V4650263" #V4730127"
#df <- fromJSON("C:\\Users\\cromagnoli\\Downloads\\response_1742754443626.json")
#df <- fromJSON("C:\\Users\\cromagnoli\\Downloads\\V4730885.json")
#df <- fromJSON("C:\\Users\\cromagnoli\\Downloads\\DFM-23009528.json")
#df <-fromJSON(paste0("http://api.v2.sondehub.org/sonde/",sonda_serial))
#df <- fromJSON(paste0("C:\\Users\\cromagnoli\\Downloads\\",sonda_serial,".json"))

# Read from SondeHub
url <- paste0("http://api.v2.sondehub.org/sonde/",sonda_serial)
destfile <- paste0(data_save_path,sonda_serial)

# Download
download.file(url, destfile, mode = "wb")

# Read and (decomp) JSON
df <- fromJSON(gzfile(destfile))

# Housekeeping del dato (Se serve)
# df$humidity <- as.numeric(df$humidity)
# df$lat <- as.numeric(df$lat)
# df$lon <- as.numeric(df$lon)
# df$alt <- as.numeric(df$alt)
# df$vel_h <- as.numeric(df$vel_h)
# df$frame <- as.integer(df$frame)
# df$temp <- as.numeric(df$temp)
# df$datetime <- as.POSIXct(df$datetime, format="%Y-%m-%dT%H:%M:%OSZ", tz="UTC")
# df$time_uploaded <- as.POSIXct(df$time_uploaded, format="%Y-%m-%dT%H:%M:%OSZ", tz="UTC")
# df$timestamp <- as.POSIXct(df$`@timestamp`, format="%Y-%m-%dT%H:%M:%OSZ", tz="UTC")
summary(df)

dfloc <- first(df[df$uploader_callsign ==receiving_station,])
coords <- strsplit(first(dfloc$uploader_position), ",")[[1]]
uploader_lat <- as.numeric(coords[1])
uploader_lon <- as.numeric(coords[2])

#  Uploader to XYZ 3D
uploader_xyz <- geo_to_cartesian(uploader_lat, uploader_lon, uploader_alt)
head(df$datetime)

# Distance, elevation and bearing (to be optimized)
df <- df %>%
  rowwise() %>%
  mutate(
    sonda_xyz = list(geo_to_cartesian(lat, lon, alt)),  # Coord 3D
    distance_3d = sqrt(sum((unlist(sonda_xyz) - uploader_xyz)^2)) / 1000,  # Distance in km
    height_diff = alt - uploader_alt,  # Sonde to  uploader altitude
    elevation_angle = atan(height_diff / sqrt(sum((unlist(sonda_xyz)[1:2] - uploader_xyz[1:2])^2))) * (180 / pi),  # Elevation angle
    bearing = calculate_bearing(uploader_lat, uploader_lon, lat, lon)  #  bearing
  )

dfmio <- df[df$uploader_callsign ==receiving_station,]
df_unique <- df[!duplicated(df[, c("bearing", "distance_3d","elevation_angle")]), ]
df$bearing
################## Simple plot
p_bsnr_solo<-ggplot(dfmio, aes(x = bearing, y = distance_3d, color = snr)) +
  #geom_point(aes(size = ifelse(uploader_callsign ==receiving_station, 5, 2))) +
  geom_point(size=1.1,alpha = 0.15) +
  scale_color_gradient(low = "blue", high = "red") +  
  scale_x_continuous(expand = expansion(0, 0), 
                     limits = c(0, 360),
                     breaks = 0:3 * 90) +
  scale_y_continuous(expand = expansion(c(0, 0.05)),limits = c(0, max(dfmio$distance_3d))) +
  theme_minimal() +coord_radial(r.axis.inside = TRUE) +
  ggtitle(paste(sonda_serial, "- Bearing vs SNR"))+
  labs(x = "Bearing", y="Distance (km)")


p_esnr_solo<-ggplot(dfmio, aes(x = 90 -elevation_angle, y = distance_3d, color = snr)) +
  scale_x_continuous(expand = expansion(0, 0), 
                     limits = c(0, 90),
                     breaks = 0:1 * 90,
                     label = c("90", "0")) +
  geom_point(size=1.1,alpha = 0.15) +
  scale_color_gradient(low = "blue", high = "red") + 
  scale_y_continuous(expand = expansion(c(0, 0.05)),limits = c(0, max(dfmio$distance_3d))) +  
  coord_radial(start = 0, end = 0.5 * pi) +
  theme_minimal()+
  ggtitle(paste(sonda_serial,"- Elevation vs. SNR"))+
  labs(x = "Elevation", y="Distance (km)")

# Comparative plot

p_btp<-ggplot(df, aes(x = bearing, y = distance_3d)) +
  geom_point(size=0.1,alpha = 0.15) +
  scale_x_continuous(expand = expansion(0, 0), 
                     limits = c(0, 360),
                     breaks = 0:3 * 90) +
  geom_point(size=0.8,alpha = 0.15) +
  #scale_color_gradient(low = "blue", high = "red") + 
  scale_y_continuous(expand = expansion(c(0, 0.05)),limits = c(0, max(df$distance_3d))) +
  theme_minimal() +coord_radial(r.axis.inside = TRUE) +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),axis.title.y = element_text(vjust = 1, margin = margin(r = 8)))+
  ggtitle("Bearing full path")+
  labs(x = "Bearing", y="Distance (km)")

p_bsnr<-ggplot(dfmio, aes(x = bearing, y = distance_3d, color = snr)) +
  #geom_point(aes(size = ifelse(uploader_callsign ==receiving_station, 5, 2))) +
  geom_point(size=1.1,alpha = 0.15) +
  scale_color_gradient(low = "blue", high = "red") +  
  scale_x_continuous(expand = expansion(0, 0), 
                     limits = c(0, 360),
                     breaks = 0:3 * 90) +
  scale_y_continuous(expand = expansion(c(0, 0.05)),limits = c(0, max(df$distance_3d))) +
  theme_minimal() +coord_radial(r_axis_inside = TRUE) +
  theme(ggplot2::unit(0.8, "lines"),plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),axis.title.y = element_text(vjust = 1, margin = margin(r = 8)))+
  guides(color = guide_colorbar(barwidth = 0.7))+
  ggtitle("Bearing vs SNR")+
  labs(x = "Bearing", y="Distance (km)")

p_esnr<-ggplot(dfmio, aes(x = 90 -elevation_angle, y = distance_3d, color = snr)) +
  scale_x_continuous(expand = expansion(0, 0), 
                     limits = c(0, 90),
                     breaks = 0:1 * 90,
                     label = c("90", "0")) +
  geom_point(size=1.1,alpha = 0.15) +
  scale_color_gradient(low = "blue", high = "red") + 
  scale_y_continuous(expand = expansion(c(0, 0.05)),limits = c(0, max(df$distance_3d))) +  
  coord_radial(start = 0, end = 0.5 * pi) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),axis.title.y = element_text(vjust = 1, margin = margin(r = 8)))+
  ggtitle("Elevation vs. SNR")+
  guides(color = guide_colorbar(barwidth = 0.7))+
  labs(x = "Elevation", y="Distance (km)")

p_etp <- ggplot(df, aes(x = 90 -elevation_angle, y = distance_3d)) +
  scale_x_continuous(expand = expansion(0, 0), 
                     limits = c(0, 90),
                     breaks = 0:1 * 90,
                     label = c("90", "0")) +
  geom_point(size=0.8,alpha = 0.15) +
  #scale_color_gradient(low = "blue", high = "red") + 
  scale_y_continuous(expand = expansion(c(0, 0.05)),limits = c(0, max(df$distance_3d))) +  
  coord_radial(start = 0, end = 0.5 * pi) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),axis.title.y = element_text(vjust = 1, margin = margin(r = 8)))+
  ggtitle("Elevation full path")+
  labs(x = "Elevation", y="Distance (km)")


grid.arrange(p_bsnr,p_btp,p_esnr, p_etp, ncol = 2,top=paste("Flight of sonde",sonda_serial,"seen by",receiving_station))

grid.lines(x = c(0.5, 0.5), y = c(0.05, 0.95), 
           gp = gpar(col = "black", lwd = 1, lty = "dotted", alpha = 0.5))
grid.lines(x = c(0.05, 0.95), y = c(0.48, 0.48), 
           gp = gpar(col = "black", lwd = 1, lty = "dotted", alpha = 0.5))


# save plot
if (saveallimages=="Y") {
  png(paste0(data_save_path,"p_comparative_",sonda_serial,".png"), width = 6, height = 6, units = "in", res = 300)
  
  grid.newpage()  
  multi_plot <- grid.arrange(p_bsnr, p_btp, p_esnr, p_etp, ncol = 2, 
                             top = paste("Flight of sonde", sonda_serial, "seen by", receiving_station))
  
  grid.lines(x = c(0.5, 0.5), y = c(0.05, 0.95), 
             gp = gpar(col = "black", lwd = 1, lty = "dotted", alpha = 0.5))
  grid.lines(x = c(0.05, 0.95), y = c(0.48, 0.48), 
             gp = gpar(col = "black", lwd = 1, lty = "dotted", alpha = 0.5))
  
  dev.off()
  ggsave(paste0(data_save_path,"p_bsnr_solo_",sonda_serial,".png"), plot = p_bsnr_solo, width = 6, height = 6, dpi = 300, bg="white")
  ggsave(paste0(data_save_path,"p_esnr_solo_",sonda_serial,".png"), plot = p_esnr_solo, width = 6, height = 6, dpi = 300, bg="white")
  
}

# Interactive plot

fig <- plot_ly(dfmio, 
               r = ~distance_3d,  
               theta = ~elevation_angle, 
               color = ~snr, 
               colors = c('blue', 'red'), 
               marker = list(size = 10),
               text = ~paste(
                 "Bearing: ", round(bearing, 2), "°\n",
                 "Elevazione: ", round(elevation_angle, 2), "°\n",
                 "Altezza: ", alt, " m\n",
                 "Distanza 3D: ", round(distance_3d, 2), " km\n",
                 "SNR: ", snr, " dB"),
               hoverinfo = 'text',
               type = 'scatterpolar', 
               mode = 'markers') 


fig <- fig %>%
  add_trace(r = c(0, max(df$distance_3d)), theta = c(0, 0), 
            type = 'scatterpolar', mode = 'lines',
            line = list(color = "black", width = 2),
            showlegend = FALSE,
            inherit = FALSE) %>%   # 👈 Evita di ereditare dati non necessari
  add_annotations(x = 0, y = max(df$distance_3d) + 0.5, text = "N", 
                  showarrow = FALSE, font = list(size = 16, color = "black")) %>%
  layout(
    title = "Grafico Polare Elevazione -Distanza 3D-SNR",
    polar = list(
      radialaxis = list(title = "Distanza in linea d'aria (km)"),
      angularaxis = list(
        title = "Bearing (°)", 
        tickvals = seq(0, 360, by = 30)  
      )
    )
  )

fig <- plot_ly(dfmio, 
               r = ~distance_3d,  
               theta = ~bearing, 
               color = ~snr, 
               colors = c('blue', 'red'), 
               marker = list(size = 10),
               text = ~paste(
                 "Bearing: ", round(bearing, 2), "°\n",
                 "Elevazione: ", round(elevation_angle, 2), "°\n",
                 "Altezza: ", alt, " m\n",
                 "Distanza 3D: ", round(distance_3d, 2), " km\n",
                 "SNR: ", snr, " dB"),
               hoverinfo = 'text',
               type = 'scatterpolar', 
               mode = 'markers') 


fig <- fig %>%
  add_trace(r = c(0, max(df$distance_3d)), theta = c(0, 0), 
            type = 'scatterpolar', mode = 'lines',
            line = list(color = "black", width = 2),
            showlegend = FALSE,
            inherit = FALSE) %>%   # 👈 Evita di ereditare dati non necessari
  add_annotations(x = 0, y = max(df$distance_3d) + 0.5, text = "N", 
                  showarrow = FALSE, font = list(size = 16, color = "black")) %>%
  layout(
    title = "Grafico Polare Bearing -Distanza 3D-SNR",
    polar = list(
      radialaxis = list(title = "Distanza in linea d'aria (km)"),
      angularaxis = list(
        title = "Bearing (°)", 
        tickvals = seq(0, 360, by = 30)  
      )
    )
  )

fig
########################
# Just test do not use #
########################
ggplot(dfmio, aes(x = bearing, y = distanza_3d, color = snr)) +
  geom_point(size = 3) +
  scale_color_viridis_c() +
  coord_polar(theta = "x") +  # Trasforma l'asse X in un sistema polare
  labs(x = "Longitudine", y = "Altitudine", color = "Latitudine") +
  theme_minimal()
ggplot(df, aes(x = heading, y = snr)) +
  geom_point(size = 1, color = "blue") +  # Punti blu
  coord_polar(theta = "x", start = 0) +   # Sistema polare
  scale_x_continuous(breaks = c(0, 90, 180, 270),  # Imposta Nord, Est, Sud, Ovest
                     labels = c("N", "E", "S", "W"), 
                     limits = c(0, 360)) +
  labs(x = "Bearing", y = "Distanza") + 
  theme_minimal()


bearing(c(uploader_lon, uploader_lat), c(df$lon, df$lat))
df <-df %>%
  mutate(bearing = bearing(c(uploader_lon, uploader_lat), cbind(lon, lat)))


df <- df %>%
  rowwise() %>%  # Per applicare una funzione per ogni riga
  mutate(
    bearing = bearing(c(uploader_lon, uploader_lat), c(lon, lat)),  # Calcola il bearing
    distance_horizontal = distVincentySphere(c(uploader_lon, uploader_lat), c(lon, lat)) / 1000,
    height_diff = alt - uploader_height,  # Differenza di altitudine in metri
    distance_3d = sqrt(distance_horizontal^2 + height_diff^2) / 1000  
    
  )
df <- df %>%
  mutate(
    elevation_angle = atan((alt - uploader_alt) / (distance_km * 1000)) * (180 / pi),  # Convertiamo l'angolo in gradi
    elevation_angle_3d = atan(height_diff / distance_horizontal) * (180 / pi)
  )
df$bearing <- df
ggplot(df, aes(x = bearing, y = distance_km)) +
  geom_point(size = 1, color = "blue") +  # Punto che rappresenta il bearing
  coord_polar(theta = "x", start = 0) +   # Sistema polare
  scale_x_continuous(
    breaks = c(0, 90, 180, 270), 
    labels = c("N", "E", "S", "W"), 
    limits = c(0, 360)
  ) +
  labs(title = "Bearing e Distance",
       x = "Bearing (°)", y = "Distance(km)") + 
  theme_minimal()


dfmio <- df[df$uploader_callsign ==receiving_station,]
ggplot(dfmio, aes(x = bearing, y = distance_km, color = snr)) +
  geom_point(size = 1) +  # Aggiungi i punti con dimensione
  coord_polar(start = 0) +  # Imposta il grafico in coordinate polari
  scale_color_gradient(low = "blue", high = "red") +  # Imposta la scala di colore per SNR
  labs(
    title = "Grafico Polare: Bearing vs Distanza con Segnale/Rumore",
    x = "Bearing (gradi)",
    y = "Distanza (km)",
    color = "SNR"
  ) +
  theme_minimal()  # Tema grafico minimale
ggplot(dfmio, aes(x = distance_km, y = alt, color = snr)) +
  geom_point(size =1 ) +  # Aggiungi i punti con dimensione
  scale_color_gradient(low = "blue", high = "red") +  # Imposta la scala di colore per SNR
  labs(
    title = "Grafico: Distanza vs Altezza con SNR",
    x = "Distanza (km)",
    y = "Altezza (m)",
    color = "SNR"
  ) +
  theme_minimal()  # Tema grafico minimale

min(df$distance_km)
ggplot(dfmio, aes(x = elevation_angle, y = distance_km, color = snr)) +
  geom_point(size = 4) +  # Aggiungi i punti
  coord_polar(start = 0) + # Imposta il grafico polare
  scale_color_gradient(low = "blue", high = "red") +  # Gradiente di colori per SNR
  xlim(0, 90) +
  theme_minimal() +
  labs(
    title = "Grafico Polare di Elevazione, Distanza e Segnale Rumore",
    x = "Angolo di Elevazione (gradi)",
    y = "Distanza (km)",
    color = "Segnale/Rumore (SNR)"
  )
ggplot(dfmio, aes(x = elevation_angle, y = distance_km, fill = snr)) +
  geom_bar(stat = "identity", width = 1) +  # Usa geom_bar per creare il grafico
  coord_polar(start = 0) +  # Imposta il grafico polare
  scale_fill_gradient(low = "blue", high = "red") +  # Gradiente di colori per SNR
  theme_minimal() +
  labs(
    title = "Grafico Polare di Elevazione, Distanza e Segnale Rumore",
    x = "Angolo di Elevazione (gradi)",
    y = "Distanza (km)",
    fill = "Segnale/Rumore (SNR)"
  )
fig <- plot_ly(dfmio, 
               r = ~distance_km, 
               theta = ~elevation_angle, 
               color = ~snr, 
               colors = c('blue', 'red'), 
               marker = list(size = 10),
               type = 'scatterpolar', 
               mode = 'markers') %>%
  layout(
    title = "Grafico Polare di Elevazione, Distanza e Segnale Rumore",
    polar = list(
      radialaxis = list(title = "Distanza (km)"),
      angularaxis = list(title = "Angolo di Elevazione (gradi)", tickvals = seq(0, 90, by = 10))
    ),
    coloraxis = list(colorbar = list(title = "Segnale/Rumore (SNR)"))
  )
fig <- plot_ly(dfmio, 
               r = ~distance_km, 
               theta = ~elevation_angle, 
               color = ~snr, 
               colors = c('blue', 'red'), 
               marker = list(size = 10),
               text = ~paste("Altezza: ", alt, " m\n",
                             "Elevazione: ", round(elevation_angle, 2), "°\n",
                             "Distanza: ", round(distance_km, 2), " km"),  # Aggiungi elevazione e distanza
               hoverinfo = 'text',  # Mostra solo il testo al passaggio del mouse
               type = 'scatterpolar', 
               mode = 'markers') %>%
  layout(
    title = "Grafico Polare di Elevazione, Distanza e Segnale Rumore",
    polar = list(
      radialaxis = list(title = "Distanza (km)"),
      angularaxis = list(title = "Angolo di Elevazione (gradi)", tickvals = seq(0, 90, by = 10))
    ),
    coloraxis = list(colorbar = list(title = "Segnale/Rumore (SNR)"))
  )



#############################################################

R <- 6371000  

# Funzione per convertire lat/lon/alt in coordinate cartesiane (X, Y, Z)
geo_to_cartesian <- function(lat, lon, alt) {
  lat_rad <- lat * pi / 180
  lon_rad <- lon * pi / 180
  x <- (R + alt) * cos(lat_rad) * cos(lon_rad)
  y <- (R + alt) * cos(lat_rad) * sin(lon_rad)
  z <- (R + alt) * sin(lat_rad)
  return(c(x, y, z))
}


# Convertiamo Villanterio in coordinate 3D
villanterio_xyz <- geo_to_cartesian(uploader_lat, uploader_lon, uploader_alt)

# Funzione per calcolare il bearing tra due punti geografici
calculate_bearing <- function(lat1, lon1, lat2, lon2) {
  lat1_rad <- lat1 * pi / 180
  lon1_rad <- lon1 * pi / 180
  lat2_rad <- lat2 * pi / 180
  lon2_rad <- lon2 * pi / 180
  
  dLon <- lon2_rad - lon1_rad
  y <- sin(dLon) * cos(lat2_rad)
  x <- cos(lat1_rad) * sin(lat2_rad) - sin(lat1_rad) * cos(lat2_rad) * cos(dLon)
  bearing_rad <- atan2(y, x)
  bearing_deg <- (bearing_rad * 180 / pi + 360) %% 360  # Convertiamo in gradi e normalizziamo tra 0-360
  return(bearing_deg)
}

df <- df %>%
  rowwise() %>%
  mutate(
    sonda_xyz = list(geo_to_cartesian(lat, lon, alt)),  # Convertiamo in coordinate 3D
    distance_3d = sqrt(sum((unlist(sonda_xyz) - villanterio_xyz)^2)) / 1000,  # Distanza in km
    height_diff = alt - uploader_alt,  # Differenza di altitudine
    elevation_angle = atan(height_diff / sqrt(sum((unlist(sonda_xyz)[1:2] - villanterio_xyz[1:2])^2))) * (180 / pi),  # Angolo di elevazione
    bearing = calculate_bearing(uploader_lat, uploader_lon, lat, lon)  # Calcoliamo il bearing
  )
#%>%
#  filter(elevation_angle >= 0 & elevation_angle <= 90)  # Filtriamo per evitare valori non validi
min(df$elevation_angle)

dfmio <- df[df$uploader_callsign ==receiving_station,]
# Creiamo il grafico polare con il bearing, distanza 3D e segnale rumore
fig <- plot_ly(dfmio, 
               r = ~distance_3d,  
               theta = ~bearing, 
               color = ~snr, 
               colors = c('blue', 'red'), 
               marker = list(size = 5),
               text = ~paste("Altezza: ", alt, " m\n",
                             "Elevazione: ", round(elevation_angle, 2), "°\n",
                             "Distanza 3D: ", round(distance_3d, 2), " km\n",
                             "Bearing: ", round(bearing, 2), "°"),
               hoverinfo = 'text',
               type = 'scatterpolar', 
               mode = 'markers') %>%
  layout(
    title = "Grafico Polare Bearing-Distanza 3D-SNR",
    polar = list(
      radialaxis = list(title = "Distanza in linea d'aria (km)"),
      angularaxis = list(
        title = "Bearing (°)", 
        tickvals = seq(0, 360, by = 30)  # Tick ogni 30 gradi
      )
    )
  )
plot_ly(dfmio, 
        r = ~distance_3d,  
        theta = ~elevation_angle, 
        color = ~snr, 
        colors = c('blue', 'red'), 
        marker = list(size = 5),
        text = ~paste("Altezza: ", alt, " m\n",
                      "Elevazione: ", round(elevation_angle, 2), "°\n",
                      "Distanza 3D: ", round(distance_3d, 2), " km\n",
                      "Bearing: ", round(bearing, 2), "°"),
        hoverinfo = 'text',
        type = 'scatterpolar', 
        mode = 'markers') %>%
  layout(
    title = "Grafico Polare Bearing-Distanza 3D-SNR",
    polar = list(
      radialaxis = list(title = "Distanza in linea d'aria (km)"),
      angularaxis = list(
        title = "Bearing (°)", 
        tickvals = seq(0, 360, by = 30)  # Tick ogni 30 gradi
      )
    )
  )

fig <- plot_ly(dfmio, 
               r = ~distance_3d,  
               theta = ~bearing, 
               color = ~snr, 
               colors = c('blue', 'red'), 
               marker = list(size = 10),
               text = ~paste("Altezza: ", alt, " m\n",
                             "Elevazione: ", round(elevation_angle, 2), "°\n",
                             "Distanza 3D: ", round(distance_3d, 2), " km\n",
                             "Bearing: ", round(bearing, 2), "°\n",
                             "SNR: ", snr, " dB"),
               hoverinfo = 'text',
               type = 'scatterpolar', 
               mode = 'markers') 

fig <- fig %>%
  add_trace(r = c(0, max(df$distance_3d)), theta = c(0, 0), 
            type = 'scatterpolar', mode = 'lines',
            line = list(color = "black", width = 2),
            showlegend = FALSE,
            inherit = FALSE) %>%   # 👈 Evita di ereditare dati non necessari
  add_annotations(x = 0, y = max(df$distance_3d) + 0.5, text = "N", 
                  showarrow = FALSE, font = list(size = 16, color = "black")) %>%
  layout(
    title = "Grafico Polare Elevazione-Distanza 3D-SNR",
    polar = list(
      radialaxis = list(title = "Distanza in linea d'aria (km)"),
      angularaxis = list(
        title = "Bearing (°)", 
        tickvals = seq(0, 360, by = 30)  
      )
    )
  )

fig
fig

p<-ggplot(dfmio, aes(x = bearing, y = distance_3d, color = snr)) +
  geom_point(alpha = 0.15) +
  scale_x_continuous(expand = expansion(0, 0), 
                     limits = c(0, 360),
                     breaks = 0:3 * 90) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  theme_minimal() +coord_radial() +
  ggtitle("Bearing vs SNR")+
  labs(x = "Bearing", y="Distance (km)")

p<-ggplot(dfmio, aes(x = 90 -elevation_angle, y = distance_3d, color = snr)) +
  geom_point(alpha = 0.15) +
  scale_x_continuous(expand = expansion(0, 0), 
                     limits = c(0, 90),
                     breaks = 0:3 * 90) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +  
  coord_radial(start = 0, end = 0.5 * pi) +
  theme_minimal()+
  theme(axis.text.x =  element_blank() )+
  ggtitle("Elevation vs. SNR")+
  labs(x = "Elevation", y="Distance (km)")

+  # Rimuove le linee di griglia
  quarter <- p+ annotate(geom="text", x = 90, y = 40, label="Annotation",color="Red")


quarter <- p + coord_radial(start = -pi/2, end = pi/2) + 
  ggtitle("0 to +90 degrees") 
+
  theme_minimal()
quarter <- p + coord_radial(start = 0, end = 0.5 * pi,reverse ="tetha") 
quarter <- p + coord_radial(start = 0, end = 0.5 * pi)+
  ggtitle("0 to +0.5π")
p + coord_radial(inner.radius = 1/5)

p<-ggplot(dfmio, aes(x = bearing, y = distance_3d)) +
  stat_density(geom = "area", fill = "grey50") +
  scale_x_continuous(expand = expansion(0, 0), 
                     limits = c(0, 360),
                     breaks = 0:3 * 90) +
  scale_y_continuous(expand = expansion(c(0, 0.05)))

p + coord_radial()
