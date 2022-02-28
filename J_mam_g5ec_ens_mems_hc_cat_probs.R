# R code to read NMME ensemble members stored in a netcdf file downloaded from
#  IRI data library and calculate ensemble probabilities
# lats, lons, years set by expert script
#
# probabilities stored in cpt format for imput into cpt probability forecast verification (PFV)
#

# install.packages(c("ncdf4"))

library(ncdf4)

tercile = TRUE # quintile if FALSE
glosea5 = TRUE # ECMWF if FALSE

if (glosea5 & tercile) probs_file='g5_ppn_ke_jan_MAM9316_etprobs.txt'
if (glosea5 & (tercile == FALSe)) probs_file='g5_ppn_ke_jan_MAM9316_eqprobs.txt'
if ((glosea5 == FALSE) & tercile) probs_file='ec_ppn_ke_jan_MAM9316_etprobs.txt'
if ((glosea5 == FALSE) & (tercile == FALSE)) probs_file='ec_ppn_ke_jan_MAM9316_eqprobs.txt'

tmths= "-01-01"

if (glosea5) hc_name <- "g5_ppn_ke_jan_MAM9316_mems"
if (glosea5 == FALSE) hc_name <- "ec_ppn_ke_jan_MAM9316_mems"
hc_fname <- paste(hc_name, ".nc", sep = "")


# open  NetCDF file containing NMME ensemble members
ncin <- nc_open(hc_fname)

print(ncin)

lons <- ncvar_get(ncin, "X")
nlon <- dim(lons)
lats <- ncvar_get(ncin, "Y")
nlat <- dim(lats)
yrs <- ncvar_get(ncin, "S")
ny <- dim(yrs)
mems <- ncvar_get(ncin,"M")
nem <- dim(mems)

if (glosea5) prec0 <- ncvar_get(nchc, "precipitation_flux")
if (glosea5 == FALSE) prec0 <- ncvar_get(nchc, "prcp")
probs <- array( , dim=c(nlon,nlat,ny))

p1 <- array(prec0, , dim=c(nlon,nlat,nem,ny))
probs <- array(0,, dim=c(nlon,nlat,3,ny))

for (lon1 in 1:nlon) { 
   for (lat1 in  1:nlat) {
     for (yr in 1:ny)  {
        p2=p1[lon1,lat1,,-yr]
        p2j=p1[lon1,lat1,,yr]
        if (tercile) tce = quantile(p2, c(0.333,0.667))
        if (tercile == FALSE) tce = quantile(p2, c(0.2,0.8))
        p3= (p2j > tce[1]) + (p2j > tce[2])
        
        probs[lon1,lat1,1,yr] = 100*sum(p3 == 0)/nem
        probs[lon1,lat1,3,yr] = 100*sum(p3 == 2)/nem
        probs[lon1,lat1,2,yr]= 100.0-(probs[lon1,lat1,1,yr]+probs[lon1,lat1,3,yr])  
     } 
   }
}

# SAVE PROBS in CPT format

lines = c()
lines = c(lines, paste("xmlns:cpt=http://iri.columbia.edu/CPT/v10/"))
lines = c(lines, paste("cpt:ncats=3"))
 
if (tercile) h1 =  "cpt:field=precipitation, cpt:C=1, cpt:clim_prob=0.333333333333, cpt:T="
if (tercile == FALSE) h1 =  "cpt:field=precipitation, cpt:C=1, cpt:clim_prob=0.2, cpt:T="
h1 = paste0(h1, toString(as.integer(yrs[1]/12+1960)))
h1 = paste0(h1, tmths)
h1 = paste(h1, ", cpt:nrow=")
h1 = paste0(h1, toString(nlat))
h1 = paste(h1, ", cpt:ncol=")
h1 = paste0(h1, toString(nlon))
h1 = paste(h1, ", cpt:row=Y, cpt:col=X, cpt:units=probability (%), cpt:missing=-1.00000000000")

h2 = paste(lons,collapse=" ")

for (yr in 1:ny) {
   if (yr == 1) {lines = c(lines, h1)
   } else {
      h4= "cpt:C=1, cpt:clim_prob=0.33333, cpt:T="
      h4 = paste0(h4, toString(as.integer(yrs[yr]/12+1960)))
      h4 = paste0(h4, tmths)
      lines= c(lines, h4)
   }
   lines =c(lines, h2) 
  
   for (lat1 in 1:nlat) {
      h3= paste(probs[,lat1,1,yr],collapse=" ")
      h3= paste(lats[lat1], h3)
      lines = c(lines,h3)
   }

   if (tercile) lines=c(lines, "cpt:C=2, cpt:clim_prob=0.333333333334")
   if (tercile == FALSE) lines=c(lines, "cpt:C=2, cpt:clim_prob=0.6")
   lines =c(lines, h2) 
   for (lat1 in 1:nlat) {
      h3= paste(probs[,lat1,2,yr],collapse=" ") 
      h3= paste(lats[lat1], h3)
      lines = c(lines,h3)
    }
   
   if (tercile) lines=c(lines, "cpt:C=3, cpt:clim_prob=0.333333333333")
   if (tercile == FALSE) lines=c(lines, "cpt:C=3, cpt:clim_prob=0.2")
   lines =c(lines, h2) 
   for (lat1 in 1:nlat) {
     h3= paste(probs[,lat1,3,yr],collapse=" ")
     h3= paste(lats[lat1], h3)
     lines = c(lines,h3)
   }
   
   
}

# open file, write lines, close file
out = file(probs_file)
writeLines(lines, out)
close(out)


