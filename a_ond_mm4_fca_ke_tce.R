######
# Combine ncal sets of calibrated GCM forecasts

i_files=c("g5_ppn_ke_aug_OND2020_ctprobs.txt",
          "ec_ppn_ke_aug_OND2020_ctprobs.txt",
          "g5_ppn_ke_aug_OND2020_etprobs.txt",
          "ec_ppn_ke_aug_OND2020_etprobs.txt")
          

num_models = 4.0

mmw=c(1.0,1.0,1.0,1.0)

calw=0.5

o_file = "mm4_ppn_ke_aug_OND2020_tprobs.txt" 

# IMPORTANT 3 INPUT FILES MUST BE ON SAME GRID

year1=2020     # first year of forecast dataset and analysis
year2=2020    # last year of forecast dataset and analysis
ny=year2-year1+1

fc_file = i_files[1] 

# READ HEADER FROM FIRST FILE AND EXTRACT NO Longs and lats

stopifnot(file.exists(fc_file))
lines_fc = readLines(fc_file)

date_lines_fc = grep(pattern='[0-9]{4}-[0-9]{2}/[0-9]{2}', lines_fc)

line3=strsplit(lines_fc[date_lines_fc[1]],split=",")

vc5 = line3[[1]][5]
vc6 = line3[[1]][6]
#
# vc5 is " cpt:nrow=60" from which we want to extract the number 60. We use
# strsplit to split the string along the '=', take the second element, and
# convert to numeric. We do the same to extract the number of columns from vc6
irows = as.numeric(strsplit(vc5, '=')[[1]][2])
icols = as.numeric(strsplit(vc6, '=')[[1]][2])

p_above = array(0, dim=c(irows,icols,ny))
p_below = array(0, dim=c(irows,icols,ny))

# Loop through all the years , store probs and obs category in same size numeric vectors

for (year in year1:year2) {
   print(year)  
  
   for (mdl in 1:num_models) {
     
      fc_file = i_files[mdl]
      if (mdl>2) mw=(1.0-calw)*2 else mw=calw*2
      mw=mw*mmw[mdl]
      
      blw = read.table(fc_file, sep="",  skip=(year-year1)*(irows*3+6)+4, nrows=irows, header=FALSE )
      b_probs0 = as.matrix(blw)
      b_probs = b_probs0[, -1]
      lats= b_probs0[,1]
      
      abv = read.table(fc_file, sep="",  skip=(year-year1)*(irows*3+6)+irows*2+8, nrows=irows, header=FALSE )
      a_probs = as.matrix(abv)
      a_probs = a_probs[, -1]
      
      
  # REPEAT EXTRACTION WITH ABOVE NORMAL PROBS (c=3)   
  
      w = (b_probs > -0.05) & (a_probs > -0.05) 
      
      p_below[,,year-year1+1] = p_below[,,year-year1+1] + b_probs*mw
      p_above[,,year-year1+1] = p_above[,,year-year1+1] + a_probs*mw
   }  
}

p_below = p_below/ num_models
p_above = p_above/ num_models

# Save combined file  
  
l_out = lines_fc[1:2]
l=3

for (year in year1:year2) {
  
   l_out = c(l_out,lines_fc[l] )
   l_out = c(l_out,lines_fc[l+1] )
   l=l+2
   for (lat1 in 1:irows) {
       h3= paste(p_below[lat1,,year-year1+1],collapse=" ")
       h3= paste(lats[lat1], h3)
       l_out = c(l_out,h3)
   }
   l=l+ irows
   
   l_out = c(l_out,lines_fc[l] )
   l_out = c(l_out,lines_fc[l+1] )
   l=l+2
   for (lat1 in 1:irows) {
     px = p_below[lat1,,year-year1+1]+p_above[lat1,,year-year1+1]
     p_normal=(100-px)*(px>0) - 1.0*(px<0)
     h3= paste(p_normal,collapse=" ")
     h3= paste(lats[lat1], h3)
     l_out = c(l_out,h3)
   }
   l=l+ irows
   
   l_out = c(l_out,lines_fc[l] )
   l_out = c(l_out,lines_fc[l+1] )
   l=l+2
   for (lat1 in 1:irows) {
     h3= paste(p_above[lat1,,year-year1+1],collapse=" ")
     h3= paste(lats[lat1], h3)
     l_out = c(l_out,h3)
   }
   l=l+ irows
   
    
}

out = file(o_file)
writeLines(l_out,out)
close(out)

