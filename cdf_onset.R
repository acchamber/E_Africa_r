infile = 'chirps_daily_waf_jjaS.tsv'
outfile = 'chirps_onset10_waf_jjas.tsv'

# onset is when rainfall reaches onset_threshold x seasonal mean

onset_threshold = 0.1 

# comment the following line to allow overwriting the output file
stopifnot(file.exists(outfile) == FALSE)


###################################################################
# read data line by line
###################################################################

stopifnot(file.exists(infile))
lines = readLines(infile)


###################################################################
# the following transforms the data into a list of matrices
###################################################################

# find lines that start with a date (2 numbers - space - 3 letters - space - 4
# numbers)
date_lines = grep(pattern='^[0-9]{2} [A-Za-z]{3} [0-9]{4}', lines)

# sanity check: are all data blocks of equal size?
block_len = date_lines[2] - date_lines[1]
stopifnot(all(diff(date_lines) == block_len))

# loop over blocks 
rain = list()
for (date_line in date_lines) {

  # extract current block and split by tabs, trim whitespace
  the_lines = lines[date_line:(date_line+block_len-1)]
  tab = strsplit(the_lines, '\t')
  tab = t(sapply(tab, trimws))

  # parse date, lats, lons, and data
  the_date = as.Date(tab[1,1], format='%d %B %Y')
  the_lons = as.numeric(tab[1, -1])
  the_lats = as.numeric(tab[-1, 1])
  the_data = apply(tab[-1, -1], 2, as.numeric)
  dimnames(the_data) = list(lat = the_lats, lon = the_lons)

  # replace NAs
  the_data[ the_data == -999 ] = NA

  # append to list
  rain[[as.character(the_date)]] = the_data

}


###################################################################
# calculate seasonal totals (per year), and mean total
###################################################################

years = format(as.Date(names(rain)), format='%Y')
uniq_years = unique(years)
totals = list()
for (yr in uniq_years) {
  inds = which(years == yr)
  arr = simplify2array(rain[inds])
  totals[[yr]] = apply(arr, c('lat', 'lon'), sum)
}

mean_total = apply(simplify2array(totals), c('lat', 'lon'), mean)


###################################################################
# loop over years and gridpoints and check on which day total rainfall exceeds
# 10% of mean 
###################################################################

years = format(as.Date(names(rain)), format='%Y')
uniq_years = unique(years)
onset = list()
for (yr in uniq_years) {
  inds = which(years == yr)
  arr = simplify2array(rain[inds])
  
  the_onset = NA * arr[, , 1]

  for (ii in 1:nrow(the_onset)) {
    for (jj in 1:ncol(the_onset)) {
      vec = arr[ii, jj, ]
      exceed_10pct = which(cumsum(vec) > onset_threshold * mean_total[ii, jj])
      the_onset[ii, jj] = exceed_10pct[1]
    }
  }

  # convert NA's to -999
  the_onset[ is.na(the_onset) ] = -999
  
  onset[[yr]] = the_onset
}


###################################################################
# write onset data to file 
###################################################################

# use date of the first entry in the input file for the beginning of the
# seasons
day = format(as.Date(names(rain)[1]), format='%d %b') # "01 Mar"
dates = paste(day, names(onset))

lines = c()
for (ii in 1:length(onset)) {

  the_onset = onset[[ii]]

  # first line: date, longitudes
  lines = c(lines, paste(c(dates[ii], colnames(the_onset)), 
                         collapse=' \t '))

  # the other lines: latitude, data
  data = cbind(rownames(the_onset), the_onset)
  lines = c(lines, apply(data, 1, paste, collapse=' \t '))

}
lines = c(lines, '') # finish with empty line

# open file, write lines, close file
out = file(outfile)
writeLines(lines, out)
close(out)

