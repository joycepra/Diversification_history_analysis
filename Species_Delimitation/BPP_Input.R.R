###############################################################################
######################## SPECIES DELIMITATION #################################

#Convert pyrad .loci output to the format used in BPP and iBPP software 
#Based on J. P. Huang scripts 



#Split samples by locus and create the input file for iBPP

data.loci <- scan('BD1_chopped&edited.loci', what = 'character', sep = '\n');#used the edited .loci file
species.vec <- c('>11', '>12', '>13', '>14','>15','>16','>17','>18','>19', '>20', '>21', '>25', '>26', '>27', '>28', 
                 '>29', '>30','>31', '>32', '>33', '>34', '>35', '>36', '>37', '>38', '>39', '>40','>41', '>42', 
                 '>43', '>44', '>45', '>46', '>47', '>49', '>50','>51', '>53', '>54', '>55', '>56', '>57', '>58', '>59', 
                 '>60','>63','>64', '>65', '>66', '>67', '>68', '>69', '>70', '>71', '>72','>73', '>74', '>75', '>76', 
                 '>77', '>78');#enter here the taxa


break.lines <- grep('//', data.loci);
index.lines <- c(0, break.lines);
loci.count <- 0;#number of loci counter. will print the number of loci at the end

for(ii in 1:87424){#the number of loci in the .loci file. Can be obtained by typing length(break.lines).
  
  num.1 <- index.lines[ii] + 1;
  num.2 <- index.lines[ii + 1] - 1;
  
  temp.text <- data.loci[c(num.1:num.2)];
  
  num.sp <- 0;
  ind.vec <- NULL;
  for(iter in 1:length(species.vec)){
    temp.ind.vec <- grep(species.vec[iter], temp.text);
    ind.vec <- c(ind.vec, temp.ind.vec);
    num.sp <- num.sp + (length(temp.ind.vec) > 0);
    
    if(num.sp == length(species.vec)){
      loci.count <- loci.count + 1;
      
      subset.temp <- temp.text[ind.vec];
      
      string.ed <- subset.temp[1];
      find.p <- '>11_[[:alnum:]]+[[:space:]]+'#for calc length of the loci; data here are all 110#change the regex according to how you name your samples
      sub.pa <- '';
      temp.dna.seq <- sub(find.p,sub.pa,string.ed);
      dna.seq <- strsplit(temp.dna.seq, '');
      loci.length <- length(dna.seq[[1]])
      n.ind <- length(subset.temp);
      first.line <- paste(n.ind, loci.length, sep = ' ');#first line for bpp each locus
      
      new.string <- NULL;
      a11.vec <- grep(species.vec[1], subset.temp);
      count <- 1;
      for(iii in a11.vec){
        s.pattern <- '>11_[[:alnum:]]+';#regex for the 1st species
        bpp.tail <- paste('^', count, sep = '');
        rep.pattern <- paste('	H', bpp.tail, sep = '');#how you want the sample to be named in bpp file
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a12.vec <- grep(species.vec[2], subset.temp);
      count <- 1;
      for(iii in a12.vec){
        s.pattern <- '>12_[[:alnum:]]+';#regex for the 2nd species
        bpp.tail <- paste('^', count+1, sep = '');#count + 1 because there are 1 sample of vec _11
        rep.pattern <- paste('	I', bpp.tail, sep = '');#how you want to name the 2nd species
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }  
      a13.vec <- grep(species.vec[3], subset.temp);
      count <- 1;
      for(iii in a13.vec){
        s.pattern <- '>13_[[:alnum:]]+';#regex for the 2nd species
        bpp.tail <- paste('^', count+2, sep = '');#count + 1 because there are 1 sample vec _11 and 1 sample vec _12
        rep.pattern <- paste('	I', bpp.tail, sep = '');#how you want to name the 2nd species
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }   
      a14.vec <- grep(species.vec[4], subset.temp);
      count <- 1;
      for(iii in a14.vec){
        s.pattern <- '>14_[[:alnum:]]+';#regex for the 2nd species
        bpp.tail <- paste('^', count+3, sep = '');
        rep.pattern <- paste('	I', bpp.tail, sep = '');#how you want to name the 2nd species
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }  
      a15.vec <- grep(species.vec[5], subset.temp);
      count <- 1;
      for(iii in a15.vec){
        s.pattern <- '>15_[[:alnum:]]+';#regex for the 2nd species
        bpp.tail <- paste('^', count+4, sep = '');
        rep.pattern <- paste('	I', bpp.tail, sep = '');#how you want to name the 2nd species
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      } 
      a16.vec <- grep(species.vec[6], subset.temp);
      count <- 1;
      for(iii in a16.vec){
        s.pattern <- '>16_[[:alnum:]]+';#regex for the 2nd species
        bpp.tail <- paste('^', count+5, sep = '');
        rep.pattern <- paste('	I', bpp.tail, sep = '');#how you want to name the 2nd species
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }   
      a17.vec <- grep(species.vec[7], subset.temp);
      count <- 1;
      for(iii in a17.vec){
        s.pattern <- '>17_[[:alnum:]]+';#regex for the 2nd species
        bpp.tail <- paste('^', count+6, sep = '');
        rep.pattern <- paste('	I', bpp.tail, sep = '');#how you want to name the 2nd species
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }   
      a18.vec <- grep(species.vec[8], subset.temp);
      count <- 1;
      for(iii in a18.vec){
        s.pattern <- '>18_[[:alnum:]]+';#regex for the 2nd species
        bpp.tail <- paste('^', count+7, sep = '');
        rep.pattern <- paste('	I', bpp.tail, sep = '');#how you want to name the 2nd species
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      } 
      a19.vec <- grep(species.vec[9], subset.temp);
      count <- 1;
      for(iii in a19.vec){
        s.pattern <- '>19_[[:alnum:]]+';#regex for the 1st species
        bpp.tail <- paste('^', count+8, sep = '');
        rep.pattern <- paste('	I', bpp.tail, sep = '');#how you want the sample to be named in bpp file
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a20.vec <- grep(species.vec[10], subset.temp);
      count <- 1;
      for(iii in a20.vec){
        s.pattern <- '>20_[[:alnum:]]+';#regex for the 2nd species
        bpp.tail <- paste('^', count+9, sep = '');
        rep.pattern <- paste('	I', bpp.tail, sep = '');#how you want to name the 2nd species
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }			
      a21.vec <- grep(species.vec[11], subset.temp);
      count <- 1;
      for(iii in a21.vec){
        s.pattern <- '>21_[[:alnum:]]+';#regex for the 3rd sp
        bpp.tail <- paste('^', count+10, sep = '');
        rep.pattern <- paste('	I', bpp.tail, sep = '');#names for the 3rd sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
    }				
      a25.vec <- grep(species.vec[12], subset.temp);
      count <- 1;
      for(iii in a25.vec){
        s.pattern <- '>25_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+11, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	B', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }			
      a26.vec <- grep(species.vec[13], subset.temp);
      count <- 1;
      for(iii in a26.vec){
        s.pattern <- '>26_[[:alnum:]]+';#regex for the 2nd species
        bpp.tail <- paste('^', count+12, sep = '');#
        rep.pattern <- paste('	D', bpp.tail, sep = '');#how you want to name the 2nd species
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }			
      a27.vec <- grep(species.vec[14], subset.temp);
      count <- 1;
      for(iii in a27.vec){
        s.pattern <- '>27_[[:alnum:]]+';#regex for the 3rd sp
        bpp.tail <- paste('^', count+13, sep = '');
        rep.pattern <- paste('	E', bpp.tail, sep = '');#names for the 3rd sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }				
      a28.vec <- grep(species.vec[15], subset.temp);
      count <- 1;
      for(iii in a28.vec){
        s.pattern <- '>28_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+14, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	C', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }		
      a29.vec <- grep(species.vec[16], subset.temp);
      count <- 1;
      for(iii in a29.vec){
        s.pattern <- '>29_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+15, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	C', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }		
      a30.vec <- grep(species.vec[17], subset.temp);
      count <- 1;
      for(iii in a30.vec){
        s.pattern <- '>30_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+18, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	C', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }			
      a31.vec <- grep(species.vec[18], subset.temp);
      count <- 1;
      for(iii in a31.vec){
        s.pattern <- '>31_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+19, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	G', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }			
      a32.vec <- grep(species.vec[19], subset.temp);
      count <- 1;
      for(iii in a32.vec){
        s.pattern <- '>32_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+21, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	G', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
       }		
      a33.vec <- grep(species.vec[20], subset.temp);
      count <- 1;
      for(iii in a33.vec){
        s.pattern <- '>33_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+23, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	G', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
       }		
      a34.vec <- grep(species.vec[21], subset.temp);
      count <- 1;
      for(iii in a34.vec){
        s.pattern <- '>34_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+24, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	G', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }		
      a35.vec <- grep(species.vec[22], subset.temp);
      count <- 1;
      for(iii in a35.vec){
        s.pattern <- '>35_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+25, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	G', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }	
      a36.vec <- grep(species.vec[23], subset.temp);
      count <- 1;
      for(iii in a36.vec){
        s.pattern <- '>36_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+27, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	G', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }	
      a37.vec <- grep(species.vec[24], subset.temp);
      count <- 1;
      for(iii in a37.vec){
        s.pattern <- '>37_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+33, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	G', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }	
      a38.vec <- grep(species.vec[25], subset.temp);
      count <- 1;
      for(iii in a38.vec){
        s.pattern <- '>38_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+34, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	G', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }	
      a39.vec <- grep(species.vec[26], subset.temp);
      count <- 1;
      for(iii in a39.vec){
        s.pattern <- '>39_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+35, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	G', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }	
      a40.vec <- grep(species.vec[27], subset.temp);
      count <- 1;
      for(iii in a40.vec){
        s.pattern <- '>40_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+37, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	G', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }	
      a41.vec <- grep(species.vec[28], subset.temp);
      count <- 1;
      for(iii in a41.vec){
        s.pattern <- '>41_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+38, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	G', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a42.vec <- grep(species.vec[29], subset.temp);
      count <- 1;
      for(iii in a42.vec){
        s.pattern <- '>42_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+42, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	G', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a43.vec <- grep(species.vec[30], subset.temp);
      count <- 1;
      for(iii in a43.vec){
        s.pattern <- '>43_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+45, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	G', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a44.vec <- grep(species.vec[31], subset.temp);
      count <- 1;
      for(iii in a44.vec){
        s.pattern <- '>44_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+46, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	G', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }	
      a45.vec <- grep(species.vec[32], subset.temp);
      count <- 1;
      for(iii in a45.vec){
        s.pattern <- '>45_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+47, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	G', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a46.vec <- grep(species.vec[33], subset.temp);
      count <- 1;
      for(iii in a46.vec){
        s.pattern <- '>46_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+48, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	G', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a47.vec <- grep(species.vec[34], subset.temp);
      count <- 1;
      for(iii in a47.vec){
        s.pattern <- '>47_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+49, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	G', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a49.vec <- grep(species.vec[35], subset.temp);
      count <- 1;
      for(iii in a49.vec){
        s.pattern <- '>49_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+53, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	I', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a50.vec <- grep(species.vec[36], subset.temp);
      count <- 1;
      for(iii in a50.vec){
        s.pattern <- '>50_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+54, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	I', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
       a51.vec <- grep(species.vec[37], subset.temp);
      count <- 1;
      for(iii in a51.vec){
        s.pattern <- '>51_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+57, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	I', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a53.vec <- grep(species.vec[38], subset.temp);
      count <- 1;
      for(iii in a53.vec){
        s.pattern <- '>53_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+59, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	I', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a54.vec <- grep(species.vec[39], subset.temp);
      count <- 1;
      for(iii in a54.vec){
        s.pattern <- '>54_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+60, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	I', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a55.vec <- grep(species.vec[40], subset.temp);
      count <- 1;
      for(iii in a55.vec){
        s.pattern <- '>55_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+61, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	I', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a56.vec <- grep(species.vec[41], subset.temp);
      count <- 1;
      for(iii in a56.vec){
        s.pattern <- '>56_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+62, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	I', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a57.vec <- grep(species.vec[42], subset.temp);
      count <- 1;
      for(iii in a57.vec){
        s.pattern <- '>57_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+63, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	I', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a58.vec <- grep(species.vec[43], subset.temp);
      count <- 1;
      for(iii in a58.vec){
        s.pattern <- '>58_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+64, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	I', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a59.vec <- grep(species.vec[44], subset.temp);
      count <- 1;
      for(iii in a59.vec){
        s.pattern <- '>59_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+65, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	I', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a60.vec <- grep(species.vec[45], subset.temp);
      count <- 1;
      for(iii in a60.vec){
        s.pattern <- '>60_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+66, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	F', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a63.vec <- grep(species.vec[46], subset.temp);
      count <- 1;
      for(iii in a63.vec){
        s.pattern <- '>63_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+67, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	A', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a64.vec <- grep(species.vec[47], subset.temp);
      count <- 1;
      for(iii in a64.vec){
        s.pattern <- '>64_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+68, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	A', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a65.vec <- grep(species.vec[48], subset.temp);
      count <- 1;
      for(iii in a65.vec){
        s.pattern <- '>65_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+69, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	A', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a66.vec <- grep(species.vec[49], subset.temp);
      count <- 1;
      for(iii in a66.vec){
        s.pattern <- '>66_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+70, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	A', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a67.vec <- grep(species.vec[50], subset.temp);
      count <- 1;
      for(iii in a67.vec){
        s.pattern <- '>67_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+71, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	A', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a68.vec <- grep(species.vec[51], subset.temp);
      count <- 1;
      for(iii in a68.vec){
        s.pattern <- '>68_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+72, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	A', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a69.vec <- grep(species.vec[52], subset.temp);
      count <- 1;
      for(iii in a69.vec){
        s.pattern <- '>69_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+75, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	A', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a70.vec <- grep(species.vec[53], subset.temp);
      count <- 1;
      for(iii in a70.vec){
        s.pattern <- '>70_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+76, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	A', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a71.vec <- grep(species.vec[54], subset.temp);
      count <- 1;
      for(iii in a71.vec){
        s.pattern <- '>71_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+77, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	A', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a72.vec <- grep(species.vec[55], subset.temp);
      count <- 1;
      for(iii in a72.vec){
        s.pattern <- '>72_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+78, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	A', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a73.vec <- grep(species.vec[56], subset.temp);
      count <- 1;
      for(iii in a73.vec){
        s.pattern <- '>73_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+79, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	A', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a74.vec <- grep(species.vec[57], subset.temp);
      count <- 1;
      for(iii in a74.vec){
        s.pattern <- '>74_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+80, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	A', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a75.vec <- grep(species.vec[58], subset.temp);
      count <- 1;
      for(iii in a75.vec){
        s.pattern <- '>75_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+81, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	A', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a76.vec <- grep(species.vec[59], subset.temp);
      count <- 1;
      for(iii in a76.vec){
        s.pattern <- '>76_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+82, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	A', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a77.vec <- grep(species.vec[60], subset.temp);
      count <- 1;
      for(iii in a77.vec){
        s.pattern <- '>77_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+83, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	A', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      a78.vec <- grep(species.vec[61], subset.temp);
      count <- 1;
      for(iii in a78.vec){
        s.pattern <- '>78_[[:alnum:]]+';#regex for the 4th
        bpp.tail <- paste('^', count+84, sep = '');#only one sample of Gmor
        rep.pattern <- paste('	A', bpp.tail, sep = '');#names for the 4th sp
        new.string <- c(new.string, sub(s.pattern, rep.pattern, subset.temp[iii]));
        count <- count + 1;
      }
      temp.bpp.file <- c(first.line, new.string, '\n');
      write.table(temp.bpp.file, append = TRUE, file='BPP.txt', quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE);#change the file name here
      
    }
  }
  
}
cat('a total of', loci.count, 'loci');#let you know how many loci are there saved in the bpp file
