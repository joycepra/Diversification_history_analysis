###################################################################################
######################## Processing Sequence Data #################################

#Visualize and filter loci to trim sites with high variation suggestive of sequencing 
#and/or alignment errors and remove suspicious clusters of possible paralogs based
#on a maximum pairwise sequence divergence 

#Based on J. P. Huang scripts (Huang, 2016)


###################################################################################

#(1) Counting variables sites from pyrad .loci output

data.loci <- scan('BD1.loci', what = 'character', sep = '\n'); #change name of the .loci file

break.lines <- grep('//', data.loci);

var.vec <- NULL;

for(i in 1:length(break.lines)){
	
	test <- data.loci[break.lines[i]]
	temp <- unlist(strsplit(test,''))
	n.var.sites <- sum(temp == '-') + sum(temp == '*')
	var.vec <- c(var.vec,n.var.sites)
	
}


#variable site position

var.site.vec <- NULL;

for(i in 1:length(break.lines)){
	
	test <- data.loci[break.lines[i]];
	temp <- unlist(strsplit(test,''));
	s1 <- which(temp == '-');
	sM <- which(temp == '*');
	s1ed <- s1 - 16; # change to the number of space&character before the first DNA character
	sMed <- sM - 16; #change  to the number of space&character before the first DNA character
	var.site.vec <- c(var.site.vec, s1ed, sMed);
	
}

par(mfrow = c(2,1));
hist(var.vec, breaks = c(seq(-1, 100, by = 1)), xlim = c(-1,100), main = 'Number of segregating sites per locus', xlab = 'Number of segregating sites');
abline(v = quantile(var.vec, 0.025), lwd = 1.5, lty = 2, col = 'red');
abline(v = quantile(var.vec, 0.5), lwd = 2, col = 'red');
abline(v = quantile(var.vec, 0.975), lwd = 1.5, lty = 2, col = 'red');
hist(var.site.vec, xlim = c(-1,150), breaks = c(-1:150), xlab = 'Position along the loci', main = 'The position of segregating sites');
abline(h = 10100, col = 'red', lwd = 2);
abline(v = 110, col = 'red', lwd = 2);
abline(v = 100, col = 'red', lty = 2, lwd = 1.5);
abline(v = 120, col = 'red', lty = 2, lwd = 1.5);

###################################################################################

# (2) Chop off the 1st site and the sites after ### in all loci

data.loci <- scan('BD1.loci', what = 'character', sep = '\n');
break.lines <- grep('//', data.loci);

for(i in 1:length(data.loci)){   #the length of the .loci file
	temp <- data.loci[i];
	s.pattern <- '(.{110}).+'; # based on the histogram built above
	r.pattern <- '\\1';
	t.txt <- sub(s.pattern, r.pattern, temp);
	data.loci[i] <- t.txt;
}
for(iter in 1:length(break.lines)){
	line.id <- break.lines[iter];
	data.loci[line.id] <- paste(data.loci[line.id], '|', sep = '');
}

write.table(data.loci, file = 'BD1_110.loci', quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE); 


###################################################################################

# (3) Plot data after chopping


stat.file <- read.csv('BD1_110.loci');

par(mfrow = c(2,2));
hist(var.vec, breaks = c(seq(-1, 100, by = 1)), xlim = c(-1,40), main = 'Number of segregating sites per locus', xlab = 'Number of segregating sites');
abline(v = quantile(var.vec, 0.025), lwd = 1.5, lty = 2, col = 'red');
abline(v = quantile(var.vec, 0.5), lwd = 2, col = 'red');
abline(v = quantile(var.vec, 0.975), lwd = 1.5, lty = 2, col = 'red');
hist(var.site.vec, xlim = c(-1,120), breaks = c(0:120), xlab = 'Position along the loci', main = 'The position of segregating sites');

#then the summaries dis of theta & pd
#par(mfrow = c(2,1));
hist(stat.file[,1], main = 'Number of taxa per locus', xlab = 'Number of taxa', xlim = c(2,21), breaks = c(-1:20));
hist(stat.file[,2]*100, main = 'Maximum p-distance', xlab = '% pairwise divergence (interspecific)', xlim = c(-0.5,15.5), breaks = c(seq(-0.5, 15.5, by = 1)));


###################################################################################

# (4) Changing names of the samples if needed

dat <- scan('BD1_110.loci', what = 'character', sep = '\n');

file.list <- list.files(, pattern = '.csv'); # need a csv file with all the samples names

for(i in 1:length(file.list)){
  temp.tax <- scan(file.list[i], what = 'character', sep ='\n');
  for(iter in 1:length(temp.tax)){
    s.pattern <- temp.tax[iter];
    temp.p <- paste(i+10, '_', sep = '');
    rep.pattern <- paste(temp.p, s.pattern, sep = '');
    dat <- gsub(s.pattern, rep.pattern, dat);
  }
}

s.pattern <- '//';
rep.pattern <- '//   ';
dat<-gsub(s.pattern, rep.pattern, dat);

write.table(dat, file = 'loci_name_changed.loci', sep = ',', quote = FALSE, row.names = FALSE, col.names = FALSE)#save it for future use


###################################################################################


#(5) editing .loci file based on per site theta within taxon & also % pairwise divergence between taxa
#Additionally, remove loci that have less than 4 taxa & those that are not variable

#install.packages("pegas")

library(pegas)

chopped.loci <- scan('loci_name_changed.loci', what = 'character', sep = '\n');#change file name
break.lines <- grep('//', chopped.loci);
index.lines <- c(0, break.lines);

loci.seq.length <- 110;#how long are the aligned sequences
#max.theta.intra <- 0.02;#dont need this here
max.div.inter <- 0.178;#based on data available in D´Elía et al. (2015)
null.data <- rep(0, length(break.lines)*2);#for creating a data matrix
stat.mat <- matrix(null.data, nrow = length(break.lines), ncol = 2);#data matrix for sum stats
colnames(stat.mat) <- c('N_taxa', 'max_p_dis');

for(i in 1:length(break.lines)){
	start.line <- index.lines[i] + 1;
	end.line <- break.lines[i] - 1;
	temp.alignment <- chopped.loci[start.line:end.line];
	
	s.pattern <- '>(.+)_.+[[:space:]]+(.+)';#
	sp.vec <- gsub(s.pattern, '\\1', temp.alignment);
	seq.vec <- gsub(s.pattern, '\\2', temp.alignment);
	uni.sp.vec <- unique(sp.vec);
	stat.mat[i,1] <- length(uni.sp.vec);	

#	theta.vec <- NULL;#calculate theta first
#	for(iter in 1:length(uni.sp.vec)){
#		temp.sp <- uni.sp.vec[iter];
#		temp.sample.index <- grep(temp.sp, sp.vec);
#		if(length(temp.sample.index) > 1){#at least 2 individual in a taxon
#			var.sites <- NULL;
#			for(ii in 1:(length(temp.sample.index)-1)){
#				tt.seq1 <- unlist(strsplit(seq.vec[temp.sample.index[ii]], ''));
#				tt.seq2 <- unlist(strsplit(seq.vec[temp.sample.index[ii + 1]], ''));
#				temp.var.sites <- which(tt.seq1 != tt.seq2 & tt.seq1 != 'N' & tt.seq1 != '-' & tt.seq2 != 'N' & tt.seq2 != '-');
#				var.sites <- c(var.sites, temp.var.sites);
#			}
#			n.segregating.sites <- length(unique(var.sites));
#			temp.theta <- theta.s(n.segregating.sites, length(temp.sample.index))/loci.seq.length;
#			theta.vec <- c(theta.vec, temp.theta);
#		}else{#only one ind for a given taxon, just fill in 0
#			theta.vec <- c(theta.vec, 0);
#		}
#	}
#	max.theta <- max(theta.vec);
#	stat.mat[i, 2] <- max.theta;
		
	p.dis.vec <- NULL;#p.dis across samples from the same or diff taxa
	for(iii in 1:(length(seq.vec)-1)){
		ttt.seq1 <- unlist(strsplit(seq.vec[iii], ''));
		for(itt in (iii+1):length(seq.vec)){
			ttt.seq2 <- unlist(strsplit(seq.vec[itt], ''));
			dis.sites <- which(ttt.seq1 != ttt.seq2 & ttt.seq1 != 'N' & ttt.seq2 != 'N' & ttt.seq1 != '-' & ttt.seq2 != '-');
			p.dis <- length(dis.sites)/loci.seq.length;
			p.dis.vec <- c(p.dis.vec, p.dis);
			#cat(iii, '\t', itt, '\n');
		}
	}
	max.p.dis <- max(p.dis.vec);
	stat.mat[i, 2] <- max.p.dis;					
	cat(i, '\n');	
}

#write.table(stat.mat, file = 'stat_chopped_loci.csv', sep = ',', quote = FALSE, row.names = FALSE);#save it for future use

less.4.taxa <- which(stat.mat[,1] < 4);#those loci that have less than 4 taxa;
#bad.theta <- which(stat.mat[,2] > max.theta.intra & stat.mat[,2] != 'NA');
bad.pd <- which(stat.mat[,2] > max.div.inter & stat.mat[,2] != 'NA');
non.var <- which(stat.mat[,2] == 0);#those that have no variable sites between taxa (so also within taxon);
unwanted.loci <- unique(c(less.4.taxa, bad.pd, non.var));

new.stat.mat <- stat.mat[-unwanted.loci,];
#write.table(new.stat.mat, file = 'afterchopped_EDloci_stat.csv', sep = ',', quote = FALSE, row.names = FALSE);#save it for future use

unwanted.lines <- NULL
for(ite in 1:length(unwanted.loci)){
	loci.id <- unwanted.loci[ite];
	start.line <- index.lines[loci.id] + 1;
	end.line <- break.lines[loci.id];
	range.lines <- start.line:end.line;
	unwanted.lines <- c(unwanted.lines, range.lines);
	cat(ite, '\n')
}

new.loci <- chopped.loci[-unwanted.lines];
write.table(new.loci, file = 'BD1_chopped&edited.loci', quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE);#change file name to your preference
